#!/usr/bin/env python3
"""
Route A: fold Q_rad/M_rad into interval residual + Rump-with-radius.

Given mid/rad CR matrices (Q±R_Q, M±R_M) from Taylor Arb assembly:

  1. Spectral radius bound ρ_R ≥ ||R||_2 (via Frobenius or power iteration).
  2. Rump PSD on A_mid with extra = ρ_R  ⇒  every matrix in A_mid ± R is PSD
     when Rump succeeds with that extra (conservative).
  3. Interval residual enclosure for (Q-λM)v with radii.

Language: strengthens Route A; not certified N(λ). Hard map unchanged.

Usage:
  python route_A_rump_radius.py
  python route_A_rump_radius.py --N1 4 --N2 2 --N3 2 --taylor-p 4
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.linalg import eigh, norm

_HERE = Path(__file__).resolve().parent
_IE = _HERE.parent / "independent_exclusion"
sys.path.insert(0, str(_IE))
sys.path.insert(0, str(_HERE))

from m3_certify import rump_psd_certificate  # noqa: E402
from route_A_rump import assemble_taylor_mid_rad  # noqa: E402
from route_A_counting import DEFAULT_Y, diagnostic_glb_lower  # noqa: E402
from cr_prototype import KAPPA1  # noqa: E402


def spectral_radius_bound(R: np.ndarray, method: str = "fro") -> float:
    """Upper bound on ||R||_2."""
    R = np.asarray(R, dtype=float)
    if method == "fro":
        return float(norm(R, ord="fro"))
    # power iteration on R^T R
    n = R.shape[0]
    v = np.ones(n) / math.sqrt(n)
    for _ in range(40):
        w = R.T @ (R @ v)
        nw = norm(w)
        if nw < 1e-300:
            return 0.0
        v = w / nw
    return float(math.sqrt(max(v @ (R.T @ (R @ v)), 0.0)))


def rump_psd_with_radius(
    Amid: np.ndarray,
    Arad: np.ndarray,
    extra_floor: float = 0.0,
) -> Dict[str, Any]:
    """
    Certify λ_min(Amid) > ρ + extra_floor with ρ ≥ ||Arad||_2.

    If successful, every real symmetric matrix in the interval hull
    [Amid − Arad, Amid + Arad] (entrywise) has λ_min > extra_floor,
    under the domination ||Δ||_2 ≤ ||Arad||_2 (entrywise |Δ| ≤ Arad).
    """
    Amid = 0.5 * (Amid + Amid.T)
    Arad = np.abs(Arad)
    rho_fro = spectral_radius_bound(Arad, "fro")
    rho_2 = spectral_radius_bound(Arad, "2")
    rho = max(rho_fro, rho_2)  # use max of bounds; fro is valid ≥ ||·||_2
    # Actually ||·||_2 ≤ ||·||_F so fro is safe upper bound; prefer fro for safety
    rho = rho_fro
    extra = rho + extra_floor
    ok, c = rump_psd_certificate(Amid, extra=extra)
    try:
        lam_min = float(eigh(Amid, eigvals_only=True, subset_by_index=[0, 0])[0])
    except Exception:
        lam_min = float("nan")
    return dict(
        rump_ok=bool(ok),
        rump_c=float(c),
        rho_F=rho_fro,
        rho_2_est=rho_2,
        extra=float(extra),
        float_lam_min=lam_min,
        margin=float(lam_min - extra) if math.isfinite(lam_min) else float("nan"),
    )


def interval_residual(
    Qmid: np.ndarray,
    Qrad: np.ndarray,
    Mmid: np.ndarray,
    Mrad: np.ndarray,
    v: np.ndarray,
    lam: float,
) -> Dict[str, float]:
    v = np.asarray(v, dtype=float).ravel()
    av = np.abs(v)
    r_mid = Qmid @ v - lam * (Mmid @ v)
    r_rad = np.abs(Qrad) @ av + abs(lam) * (np.abs(Mrad) @ av)
    return dict(
        mid_norm=float(norm(r_mid)),
        rad_norm=float(norm(r_rad)),
        encl_norm=float(norm(np.abs(r_mid) + r_rad)),
        max_component_encl=float(np.max(np.abs(r_mid) + r_rad)),
    )


def run(
    N1: int = 4,
    N2: int = 2,
    N3: int = 2,
    Y: float = DEFAULT_Y,
    taylor_p: int = 4,
) -> Dict[str, Any]:
    t0 = time.time()
    print(
        f"=== Route A Rump-with-radius  {N1}x{N2}x{N3} Y={Y} p={taylor_p} ==="
    )
    pack = assemble_taylor_mid_rad(N1, N2, N3, Y=Y, taylor_p=taylor_p)
    Q, Rq = pack["Q_mid"], pack["Q_rad"]
    M, Rm = pack["M_mid"], pack["M_rad"]
    top = pack["top"]

    free = np.ones(Q.shape[0], dtype=bool)
    free[top] = False
    idx = np.where(free)[0]
    Qf, Mf = Q[np.ix_(idx, idx)], M[np.ix_(idx, idx)]
    Rqf, Rmf = Rq[np.ix_(idx, idx)], Rm[np.ix_(idx, idx)]

    # float spectrum
    wN = eigh(Q, M, eigvals_only=True, subset_by_index=[0, min(5, Q.shape[0] - 1)])
    wD = eigh(Qf, Mf, eigvals_only=True, subset_by_index=[0, min(4, idx.size - 1)])
    first_pos = next((float(x) for x in wN if x > 1e-6), float("nan"))
    first_D = float(wD[0])

    # Rump-with-radius on Dirichlet pencils Qf - t Mf
    ts = [0.0, 1.0, 0.5 * first_D, 0.9 * first_D]
    rump_rows = []
    for t in ts:
        Amid = Qf - t * Mf
        Arad = Rqf + abs(t) * Rmf  # |Δ(Q - tM)| ≤ R_Q + |t| R_M
        rec = rump_psd_with_radius(Amid, Arad)
        rec["t"] = float(t)
        rec["label"] = "dirichlet"
        rump_rows.append(rec)
        print(
            f"  Dir t={t:.4f}  rump_ok={rec['rump_ok']}  "
            f"ρ_F={rec['rho_F']:.3e}  λmin={rec['float_lam_min']:.3e}  "
            f"margin={rec['margin']:.3e}",
            flush=True,
        )

    # pure mid Rump (extra=0) for comparison
    mid_only = rump_psd_with_radius(Qf, np.zeros_like(Qf))
    mid_only["t"] = 0.0
    mid_only["label"] = "dirichlet_mid_only"

    # residual for first Dirichlet eigenvector
    _, vecs = eigh(Qf, Mf, subset_by_index=[0, 0])
    v = vecs[:, 0]
    resid = interval_residual(Qf, Rqf, Mf, Rmf, v, first_D)

    out = dict(
        language=(
            "Rump PSD with extra ≥ ||Q_rad|+|t|M_rad||_F so the whole "
            "interval matrix hull is PSD when successful. Taylor Arb mid/rad. "
            "Not certified counting."
        ),
        mesh=dict(N1=N1, N2=N2, N3=N3, Y=Y),
        taylor_p=taylor_p,
        n_dofs=pack["n_dofs"],
        n_dirichlet=int(idx.size),
        worst_I1_rel=pack["worst_I1_rel"],
        worst_I3_rel=pack["worst_I3_rel"],
        Q_rad_fro=float(norm(Rq, ord="fro")),
        M_rad_fro=float(norm(Rm, ord="fro")),
        lam_N=[float(x) for x in wN],
        lam_D=[float(x) for x in wD],
        first_positive_N=first_pos,
        first_D=first_D,
        glb_sketch=diagnostic_glb_lower(first_pos, pack["h_max_upper"]),
        rump_with_radius=rump_rows,
        rump_mid_only_t0=mid_only,
        rump_radius_any=any(r["rump_ok"] for r in rump_rows),
        interval_residual_first_D=resid,
        KAPPA1=float(KAPPA1),
        seconds=time.time() - t0,
        status="YELLOW",
        hard_map_note="counting_certified remains false",
    )
    print(
        f"interval residual first D: encl={resid['encl_norm']:.3e}  "
        f"rad={resid['rad_norm']:.3e}"
    )
    print(
        f"rump_radius any PASS? {out['rump_radius_any']}  "
        f"mid-only t0 PASS? {mid_only['rump_ok']}"
    )
    print(f"({out['seconds']:.2f}s)")
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--N1", type=int, default=4)
    p.add_argument("--N2", type=int, default=2)
    p.add_argument("--N3", type=int, default=2)
    p.add_argument("--Y", type=float, default=DEFAULT_Y)
    p.add_argument("--taylor-p", type=int, default=4)
    p.add_argument("--json-out", type=str, default="")
    args = p.parse_args(argv)

    out = run(args.N1, args.N2, args.N3, Y=args.Y, taylor_p=args.taylor_p)
    path = (
        Path(args.json_out)
        if args.json_out
        else _HERE / "route_A_rump_radius_result.json"
    )
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {path}")
    # infrastructure: assembly ok + mid-only Rump on Dirichlet
    ok = out["rump_mid_only_t0"]["rump_ok"]
    print("PASS infrastructure" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
