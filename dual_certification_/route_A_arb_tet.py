#!/usr/bin/env python3
"""
Route A: Arb per-tet weight enclosures for Picard core Q/M (partial Arb assembly).

Uses m3_certify.tet_arb_data + flint arb for:
  - tet volume > 0 (certified)
  - h_T enclosure
  - quadrature weight mid/rad for wQ ~ 1/y, wM ~ 1/y^3 (vertex min/max height)

Assembles float mid + explicit radius matrices for sparse CR assembly on the
level-1 Picard mesh (route A counting core — single copy, free Neumann top).

This is a step up from relative float mid/rad: geometric quantities are
Arb-enclosed. Full Rump PSD / certified eigenspectrum remains open.

Language: addresses Route A Arb path; not certified counting.

Usage:
  python route_A_arb_tet.py
  python route_A_arb_tet.py --N1 4 --N2 2 --N3 2
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
from scipy.sparse import coo_matrix

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
_IE = _ROOT / "independent_exclusion"
sys.path.insert(0, str(_IE))
sys.path.insert(0, str(_HERE))

import flint  # noqa: E402
from flint import arb  # noqa: E402

from cr_prototype import build_mesh, geometry, KAPPA1  # noqa: E402
from m3_certify import tet_arb_data, mid_rad, upper, amax, amin  # noqa: E402
from route_A_counting import (  # noqa: E402
    float_neumann_eigs,
    float_dirichlet_eigs,
    candidate_brackets,
    diagnostic_glb_lower,
    DEFAULT_Y,
)


flint.ctx.prec = 128


def arb_weights_on_tet(P: np.ndarray) -> Dict[str, Any]:
    """
    Enclose hyperbolic CR weights on one tet with vertex matrix P (4x3).

    wQ ∝ ∫_T y^{-1}, wM ∝ ∫_T y^{-3} — use degree-1 safe bounds:
      vol * min(y)^{-k} ≤ ∫ y^{-k} ≤ vol * max(y)^{-k} wait no:
      for y^{-1}: vol/ymax ≤ ∫ y^{-1} ≤ vol/ymin
      for y^{-3}: vol/ymax^3 ≤ ∫ ≤ vol/ymin^3
    Then scale into element Se, Me as in ref_elements / CR assembly.
    """
    vol_a, grads, hT = tet_arb_data(P)
    ys = [arb(float(P[i, 2])) for i in range(4)]
    ymin, ymax = ys[0], ys[0]
    for y in ys[1:]:
        ymin = amin(ymin, y)
        ymax = amax(ymax, y)
    # safe integral bounds
    I1_lo = vol_a / ymax
    I1_hi = vol_a / ymin
    I3_lo = vol_a / (ymax ** 3)
    I3_hi = vol_a / (ymin ** 3)
    # mid/rad for use as wQ, wM scalars (mean of bounds)
    def ball_from_lohi(lo, hi):
        mid = (lo + hi) / 2
        rad = (hi - lo) / 2
        return mid, rad

    wQ_m, wQ_r = ball_from_lohi(I1_lo, I1_hi)
    wM_m, wM_r = ball_from_lohi(I3_lo, I3_hi)
    # grads float mid for assembly (geometry exact doubles)
    g_mid = np.array(
        [[float(grads[a][k].mid()) for k in range(3)] for a in range(4)],
        dtype=float,
    )
    vol_m, vol_r = mid_rad(vol_a)
    h_m, h_r = mid_rad(hT)
    return dict(
        vol_mid=vol_m,
        vol_rad=vol_r,
        hT_mid=h_m,
        hT_rad=h_r,
        wQ_mid=float(wQ_m.mid()) if hasattr(wQ_m, "mid") else float(wQ_m),
        wQ_rad=float(wQ_r.mid()) + float(wQ_r.rad()) if hasattr(wQ_r, "rad") else float(wQ_r),
        wM_mid=float(wM_m.mid()) if hasattr(wM_m, "mid") else float(wM_m),
        wM_rad=float(wM_r.mid()) + float(wM_r.rad()) if hasattr(wM_r, "rad") else float(wM_r),
        grads_mid=g_mid,
        ymin=float(ymin.mid()),
        ymax=float(ymax.mid()),
        vol_positive=bool(vol_a > 0),
    )


def _to_float_mid_rad(a) -> Tuple[float, float]:
    if hasattr(a, "mid"):
        return mid_rad(a)
    return float(a), 0.0


def assemble_arb_mid_rad(
    N1: int,
    N2: int,
    N3: int,
    Y: float = DEFAULT_Y,
) -> Dict[str, Any]:
    """
    Build CR Q,M mid + radius (dense for small meshes) from Arb tet weights.

    Face table matches cr_prototype.assemble structure.
    """
    mesh = build_mesh(N1, N2, N3, Y, curved=True)
    geo = geometry(mesh)
    X, tets = mesh["X"], mesh["tets"]
    nt = len(tets)
    nfr = geo["nfaces"]
    tf = geo["tet_faces"]

    # local CR mass pattern
    Mloc = np.full((4, 4), -1 / 20.0) + np.eye(4) * (9 / 20.0)

    Q_mid = np.zeros((nfr, nfr))
    Q_rad = np.zeros((nfr, nfr))
    M_mid = np.zeros((nfr, nfr))
    M_rad = np.zeros((nfr, nfr))

    h_max_mid = 0.0
    h_max_rad = 0.0
    n_pos = 0
    worst_wQ_rel = 0.0
    worst_wM_rel = 0.0

    for e in range(nt):
        P = X[list(tets[e])]
        d = arb_weights_on_tet(P)
        if d["vol_positive"]:
            n_pos += 1
        h_max_mid = max(h_max_mid, d["hT_mid"])
        h_max_rad = max(h_max_rad, d["hT_mid"] + d["hT_rad"])

        wQ, rQ = d["wQ_mid"], d["wQ_rad"]
        wM, rM = d["wM_mid"], d["wM_rad"]
        if wQ > 0:
            worst_wQ_rel = max(worst_wQ_rel, rQ / wQ)
        if wM > 0:
            worst_wM_rel = max(worst_wM_rel, rM / wM)

        gphi = -3.0 * d["grads_mid"]  # CR gradient factor as float
        # Se = wQ * gphi @ gphi.T  (note: full assembly uses vol in weights;
        # our I1 already includes vol, so Se = I1 * (gphi@gphi.T) / something?
        # Match cr_prototype: Se = wQ * vol * (gphi@gphi.T) with wQ=mean 1/y.
        # Here I1 ≈ ∫ y^{-1} = wQ_mean * vol, so Se = I1 * (gphi@gphi.T).
        GGT = gphi @ gphi.T
        Se_mid = wQ * GGT
        Se_rad = rQ * np.abs(GGT)
        Me_mid = wM * Mloc
        Me_rad = rM * np.abs(Mloc)

        faces = tf[e]
        for a in range(4):
            for b in range(4):
                ia, ib = int(faces[a]), int(faces[b])
                Q_mid[ia, ib] += Se_mid[a, b]
                Q_rad[ia, ib] += Se_rad[a, b]
                M_mid[ia, ib] += Me_mid[a, b]
                M_rad[ia, ib] += Me_rad[a, b]

    # top face dofs
    top = sorted({fid for (_e, fid) in geo.get("top_tets", [])}) if "top_tets" in geo else []
    # re-get top from assemble if needed
    from cr_prototype import assemble as cr_assemble

    _Q, _M, _a, _t, top_tets, _floor = cr_assemble(mesh, geo)
    top = sorted({fid for (_e, fid) in top_tets})

    return dict(
        Q_mid=Q_mid,
        Q_rad=Q_rad,
        M_mid=M_mid,
        M_rad=M_rad,
        n_dofs=nfr,
        n_tets=nt,
        n_vol_positive=n_pos,
        h_max_mid=h_max_mid,
        h_max_encl_upper=h_max_rad,
        worst_wQ_rel=worst_wQ_rel,
        worst_wM_rel=worst_wM_rel,
        top_dofs=top,
        mesh=dict(N1=N1, N2=N2, N3=N3, Y=Y),
        KAPPA1=float(KAPPA1),
    )


def run(
    N1: int = 4,
    N2: int = 2,
    N3: int = 2,
    Y: float = DEFAULT_Y,
    neigs: int = 5,
) -> Dict[str, Any]:
    t0 = time.time()
    print(f"=== Route A Arb tet assembly  mesh {N1}x{N2}x{N3} Y={Y} ===")
    pack = assemble_arb_mid_rad(N1, N2, N3, Y=Y)
    Q, Rq = pack["Q_mid"], pack["Q_rad"]
    M, Rm = pack["M_mid"], pack["M_rad"]
    top = pack["top_dofs"]

    # float spectrum on mid
    lam_N = float_neumann_eigs(Q, M, neigs=neigs)
    lam_D = float_dirichlet_eigs(Q, M, top, neigs=neigs)
    brackets = candidate_brackets(lam_N, lam_D)
    h_u = pack["h_max_encl_upper"]
    glb_N = [diagnostic_glb_lower(float(x), h_u) for x in lam_N]

    # interval residual for first positive mode
    from scipy.linalg import eigh

    w, vecs = eigh(Q, M, subset_by_index=[0, min(neigs, Q.shape[0] - 1)])
    idx = 0
    for i, wi in enumerate(w):
        if wi > 1e-8:
            idx = i
            break
    v = vecs[:, idx]
    av = np.abs(v)
    r_mid = Q @ v - float(w[idx]) * (M @ v)
    r_rad = Rq @ av + abs(float(w[idx])) * (Rm @ av)

    out = dict(
        language=(
            "Partial Arb assembly: tet vol/hT certified positive; y^{-1},y^{-3} "
            "integrals enclosed by min/max height bounds; CR mid/rad matrices "
            "assembled. Float D/N on mid. Not Rump / not certified N(λ)."
        ),
        mesh=pack["mesh"],
        n_dofs=pack["n_dofs"],
        n_tets=pack["n_tets"],
        n_vol_positive=pack["n_vol_positive"],
        h_max_mid=pack["h_max_mid"],
        h_max_encl_upper=pack["h_max_encl_upper"],
        worst_wQ_rel=pack["worst_wQ_rel"],
        worst_wM_rel=pack["worst_wM_rel"],
        KAPPA1=pack["KAPPA1"],
        lam_N=[float(x) for x in lam_N],
        lam_D=[float(x) for x in lam_D],
        glb_sketch_N=glb_N,
        brackets=[(k, lo, hi) for k, lo, hi in brackets],
        interval_residual=dict(
            mid_norm=float(np.linalg.norm(r_mid)),
            rad_norm=float(np.linalg.norm(r_rad)),
            encl_norm=float(np.linalg.norm(np.abs(r_mid) + r_rad)),
        ),
        Q_rad_fro=float(np.linalg.norm(Rq, ord="fro")),
        M_rad_fro=float(np.linalg.norm(Rm, ord="fro")),
        seconds=time.time() - t0,
        status="YELLOW",
        hard_map_note="counting_certified remains false",
    )

    print(
        f"tets={out['n_tets']} dofs={out['n_dofs']}  "
        f"vol>0 certified: {out['n_vol_positive']}/{out['n_tets']}"
    )
    print(
        f"h_max mid={out['h_max_mid']:.5f}  encl_upper={out['h_max_encl_upper']:.5f}"
    )
    print(
        f"worst rel rad wQ={out['worst_wQ_rel']:.3e}  wM={out['worst_wM_rel']:.3e}"
    )
    print(f"‖Q_rad‖_F={out['Q_rad_fro']:.3e}  ‖M_rad‖_F={out['M_rad_fro']:.3e}")
    print(f"{'k':>4}  {'λ^N':>12}  {'GLB':>12}  {'λ^D':>12}")
    for (k, lo, hi), glb in zip(out["brackets"], out["glb_sketch_N"]):
        print(f"{k:4d}  {lo:12.6f}  {glb:12.6f}  {hi:12.6f}")
    print(
        f"interval residual encl_norm={out['interval_residual']['encl_norm']:.3e}"
    )
    print(f"({out['seconds']:.2f}s)  status={out['status']}")
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--N1", type=int, default=4)
    p.add_argument("--N2", type=int, default=2)
    p.add_argument("--N3", type=int, default=2)
    p.add_argument("--Y", type=float, default=DEFAULT_Y)
    p.add_argument("--neigs", type=int, default=5)
    p.add_argument("--json-out", type=str, default="")
    args = p.parse_args(argv)

    out = run(args.N1, args.N2, args.N3, Y=args.Y, neigs=args.neigs)
    path = Path(args.json_out) if args.json_out else _HERE / "route_A_arb_tet_result.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {path}")
    ok = (
        out["n_vol_positive"] == out["n_tets"]
        and all(lo <= hi + 1e-6 for _, lo, hi in out["brackets"])
    )
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
