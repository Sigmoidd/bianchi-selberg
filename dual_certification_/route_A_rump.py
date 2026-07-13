#!/usr/bin/env python3
"""
Route A: Taylor-moment Arb weights + Rump PSD probe on Picard core.

1. Per-tet ∫ y^{-1}, ∫ y^{-3} via Taylor expansion around mean height
   (same spirit as m3_certify.ell_entries / weighted moments).
2. Assemble dense Q,M mid+rad on level-1 CR mesh.
3. Rump BIT 46 PSD certificate on shifted pencils:
     A_t = Q_mid - t * M_mid   (try t below first positive float eigenvalue)
   Success ⇒ λ_min(Q - t M) > 0 in the float sense of Rump on the mid matrix
   (radii not yet folded into a full interval eigensolver).

Language: strengthens Route A Arb path; NOT certified counting N(λ)=[0,0].
Hard map unchanged.

Usage:
  python route_A_rump.py
  python route_A_rump.py --N1 4 --N2 2 --N3 2 --taylor-p 4
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
import time
from math import factorial
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.linalg import eigh

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
_IE = _ROOT / "independent_exclusion"
sys.path.insert(0, str(_IE))
sys.path.insert(0, str(_HERE))

import flint  # noqa: E402
from flint import arb  # noqa: E402

from cr_prototype import build_mesh, geometry, assemble, KAPPA1  # noqa: E402
from m3_certify import (  # noqa: E402
    tet_arb_data,
    mid_rad,
    upper,
    amax,
    amin,
    rump_psd_certificate,
)
from route_A_counting import DEFAULT_Y, diagnostic_glb_lower  # noqa: E402

flint.ctx.prec = 128


def taylor_integral_y_power(
    yv: List[Any],
    vol: Any,
    power: int,
    p: int = 4,
) -> Any:
    """
    Arb enclosure of ∫_T y^{-power} dx via Taylor of g(y)=y^{-power} at ybar.

    Uses multinomial expansion of (y-ybar)^k moments of barycentric
    coordinates on the tet (exact for polynomials ≤3 via 6V formulas).
    Remainder: |R_p| ≤ C_p ymin^{-power-p-1} Dmax^{p+1} * vol.
    """
    ybar = sum(yv, arb(0)) / 4
    d = [y - ybar for y in yv]
    # g^{(k)}(ybar)/k! for g = y^{-s}, s=power
    # g^{(k)} = (-1)^k (s)(s+1)...(s+k-1) y^{-s-k}
    out = arb(0)
    s = power
    for k in range(p + 1):
        if k == 0:
            gk = ybar ** (-s)
        else:
            rising = arb(1)
            for j in range(k):
                rising *= arb(s + j)
            gk = arb((-1) ** k) * rising / factorial(k) * ybar ** (-s - k)
        # ∫ (y-ybar)^k = 6V * sum_{|α|=k} (prod d_i^{α_i}/α_i!) * k! / (k+3)!
        # more carefully: E[λ^α] for Dirichlet(1,1,1,1)
        mom = arb(0)
        for alpha in itertools.product(range(k + 1), repeat=4):
            if sum(alpha) != k:
                continue
            mono = arb(factorial(k))
            for i in range(4):
                mono *= (d[i] ** alpha[i]) / factorial(alpha[i])
            # ∫ λ^α = 6V * prod(α_i!) / (k+3)!
            mom += mono * arb(factorial(k)) / factorial(k + 3)
        out += gk * 6 * vol * mom

    ymin = yv[0]
    for y in yv[1:]:
        ymin = amin(ymin, y)
    ymax = yv[0]
    for y in yv[1:]:
        ymax = amax(ymax, y)
    Dmax = amax(ymax - ybar, ybar - ymin)
    # remainder bound: |g^{(p+1)}(ξ)|/(p+1)! * D^{p+1}
    # |g^{(p+1)}| ≤ (s)(s+1)...(s+p) ymin^{-s-p-1}
    rising = arb(1)
    for j in range(p + 1):
        rising *= arb(s + j)
    rem_coef = rising / factorial(p + 1) * ymin ** (-s - p - 1)
    rem = rem_coef * (Dmax ** (p + 1)) * vol
    ball = arb(0).union(rem).union(-rem)
    return out + ball


def assemble_taylor_mid_rad(
    N1: int,
    N2: int,
    N3: int,
    Y: float = DEFAULT_Y,
    taylor_p: int = 4,
) -> Dict[str, Any]:
    mesh = build_mesh(N1, N2, N3, Y, curved=True)
    geo = geometry(mesh)
    X, tets = mesh["X"], mesh["tets"]
    nt = len(tets)
    nfr = geo["nfaces"]
    tf = geo["tet_faces"]
    Mloc = np.full((4, 4), -1 / 20.0) + np.eye(4) * (9 / 20.0)

    Q_mid = np.zeros((nfr, nfr))
    Q_rad = np.zeros((nfr, nfr))
    M_mid = np.zeros((nfr, nfr))
    M_rad = np.zeros((nfr, nfr))

    worst_I1_rel = 0.0
    worst_I3_rel = 0.0
    h_max_u = 0.0
    n_ok = 0

    for e in range(nt):
        P = X[list(tets[e])]
        vol_a, grads, hT = tet_arb_data(P)
        yv = [arb(float(P[i, 2])) for i in range(4)]
        I1 = taylor_integral_y_power(yv, vol_a, power=1, p=taylor_p)
        I3 = taylor_integral_y_power(yv, vol_a, power=3, p=taylor_p)
        I1m, I1r = mid_rad(I1)
        I3m, I3r = mid_rad(I3)
        if abs(I1m) > 0:
            worst_I1_rel = max(worst_I1_rel, I1r / abs(I1m))
        if abs(I3m) > 0:
            worst_I3_rel = max(worst_I3_rel, I3r / abs(I3m))
        h_max_u = max(h_max_u, upper(hT))
        n_ok += 1 if bool(vol_a > 0) else 0

        gphi = -3.0 * np.array(
            [[float(grads[a][k].mid()) for k in range(3)] for a in range(4)]
        )
        GGT = gphi @ gphi.T
        Se_mid = I1m * GGT
        Se_rad = I1r * np.abs(GGT)
        Me_mid = I3m * Mloc
        Me_rad = I3r * np.abs(Mloc)
        faces = tf[e]
        for a in range(4):
            for b in range(4):
                ia, ib = int(faces[a]), int(faces[b])
                Q_mid[ia, ib] += Se_mid[a, b]
                Q_rad[ia, ib] += Se_rad[a, b]
                M_mid[ia, ib] += Me_mid[a, b]
                M_rad[ia, ib] += Me_rad[a, b]

    _Q, _M, _a, _t, top_tets, _ = assemble(mesh, geo)
    top = sorted({fid for (_e, fid) in top_tets})
    return dict(
        Q_mid=Q_mid,
        Q_rad=Q_rad,
        M_mid=M_mid,
        M_rad=M_rad,
        n_dofs=nfr,
        n_tets=nt,
        n_vol_ok=n_ok,
        h_max_upper=h_max_u,
        worst_I1_rel=worst_I1_rel,
        worst_I3_rel=worst_I3_rel,
        top=top,
        taylor_p=taylor_p,
    )


def rump_shift_scan(
    Q: np.ndarray,
    M: np.ndarray,
    t_values: Sequence[float],
    label: str = "full",
) -> List[Dict[str, Any]]:
    """Rump PSD certificate attempts on A = Q - t M for each t."""
    rows = []
    for t in t_values:
        A = Q - float(t) * M
        A = 0.5 * (A + A.T)
        ok, c = rump_psd_certificate(A, extra=0.0)
        try:
            w = eigh(A, eigvals_only=True, subset_by_index=[0, 0])[0]
        except Exception:
            w = float("nan")
        rows.append(
            dict(
                label=label,
                t=float(t),
                rump_psd=bool(ok),
                rump_c=float(c),
                float_lam_min=float(w),
            )
        )
    return rows


def run(
    N1: int = 4,
    N2: int = 2,
    N3: int = 2,
    Y: float = DEFAULT_Y,
    taylor_p: int = 4,
) -> Dict[str, Any]:
    t0 = time.time()
    print(
        f"=== Route A Taylor+Rump  mesh {N1}x{N2}x{N3} Y={Y} p={taylor_p} ==="
    )
    pack = assemble_taylor_mid_rad(N1, N2, N3, Y=Y, taylor_p=taylor_p)
    Q, M = pack["Q_mid"], pack["M_mid"]
    top = pack["top"]

    # float spectrum on mid
    wN = eigh(Q, M, eigvals_only=True, subset_by_index=[0, min(5, Q.shape[0] - 1)])
    free = np.ones(Q.shape[0], dtype=bool)
    free[top] = False
    idx = np.where(free)[0]
    Qf, Mf = Q[np.ix_(idx, idx)], M[np.ix_(idx, idx)]
    wD = eigh(Qf, Mf, eigvals_only=True, subset_by_index=[0, min(4, idx.size - 1)])

    first_pos = next((float(x) for x in wN if x > 1e-6), float("nan"))
    first_D = float(wD[0]) if len(wD) else float("nan")

    # Neumann Q has a near-null constant mode ⇒ strict PSD at t=0 fails.
    # (a) tiny mass shift Q + ε M  (certifies λ ≥ -ε in the mass metric)
    # (b) Dirichlet restriction (no constant null) — primary Rump target
    eps = 1e-10
    rump_neumann = rump_shift_scan(Q, M, [0.0], label="neumann")
    rump_shifted = rump_shift_scan(Q + eps * M, M, [0.0], label="neumann_plus_epsM")
    ts_D = [0.0, 1.0]
    if math.isfinite(first_D):
        ts_D += [0.25 * first_D, 0.5 * first_D, 0.9 * first_D, max(first_D - 1e-3, 0.0)]
    rump_dirichlet = rump_shift_scan(Qf, Mf, ts_D, label="dirichlet")

    glb = (
        diagnostic_glb_lower(first_pos, pack["h_max_upper"])
        if first_pos == first_pos
        else float("nan")
    )
    rump_rows = rump_neumann + rump_shifted + rump_dirichlet
    dirichlet_t0 = next(r for r in rump_dirichlet if r["t"] == 0.0)

    out = dict(
        language=(
            "Taylor-moment Arb integrals for y^{-1}, y^{-3} on Picard tets; "
            "float CR mid matrices; Rump PSD on Dirichlet Q-tM and on "
            "Neumann Q+εM. Neumann t=0 fails by design (constant null). "
            "Not certified N(λ)."
        ),
        mesh=dict(N1=N1, N2=N2, N3=N3, Y=Y),
        taylor_p=taylor_p,
        n_dofs=pack["n_dofs"],
        n_tets=pack["n_tets"],
        n_vol_ok=pack["n_vol_ok"],
        h_max_upper=pack["h_max_upper"],
        worst_I1_rel=pack["worst_I1_rel"],
        worst_I3_rel=pack["worst_I3_rel"],
        lam_N=[float(x) for x in wN],
        lam_D=[float(x) for x in wD],
        first_positive_N=first_pos,
        first_D=first_D,
        glb_sketch=glb,
        KAPPA1=float(KAPPA1),
        rump_scan=rump_rows,
        rump_dirichlet_t0=bool(dirichlet_t0["rump_psd"]),
        rump_any_dirichlet_positive_t=any(
            r["rump_psd"] and r["t"] > 1e-12 for r in rump_dirichlet
        ),
        Q_rad_fro=float(np.linalg.norm(pack["Q_rad"], ord="fro")),
        M_rad_fro=float(np.linalg.norm(pack["M_rad"], ord="fro")),
        seconds=time.time() - t0,
        status="YELLOW",
        hard_map_note="counting_certified remains false",
    )

    print(
        f"tets={out['n_tets']} dofs={out['n_dofs']}  "
        f"vol ok={out['n_vol_ok']}/{out['n_tets']}"
    )
    print(
        f"Taylor p={taylor_p}  worst rel I1={out['worst_I1_rel']:.3e}  "
        f"I3={out['worst_I3_rel']:.3e}"
    )
    print(f"‖Q_rad‖_F={out['Q_rad_fro']:.3e}  ‖M_rad‖_F={out['M_rad_fro']:.3e}")
    print(f"λ^N={out['lam_N']}")
    print(f"first pos N={first_pos:.6f}  first D={first_D:.6f}  GLB sketch={glb:.6f}")
    print("Rump PSD:")
    for r in rump_rows:
        print(
            f"  [{r['label']:18s}] t={r['t']:.6f}  rump={r['rump_psd']}  "
            f"float_λmin={r['float_lam_min']:.6e}"
        )
    print(f"({out['seconds']:.2f}s)  status={out['status']}")
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
    path = Path(args.json_out) if args.json_out else _HERE / "route_A_rump_result.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {path}")
    ok = out["n_vol_ok"] == out["n_tets"] and bool(out.get("rump_dirichlet_t0"))
    print(
        "PASS"
        if ok
        else "FAIL (need volumes + Dirichlet t=0 Rump PSD)"
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
