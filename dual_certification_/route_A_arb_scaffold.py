#!/usr/bin/env python3
"""
Route A Arb / D–N scaffolding (still not certified counting).

Closes the next engineering step after truncation_constants_scaffold:

  1. Float Picard D/N spectrum (reuse route_A_counting)
  2. Interval mid/rad enclosures of Q,M via relative float radii
  3. Interval residual diagnostics for matvec
  4. Float CR-GLB sketch (KAPPA1) on Neumann/Dirichlet
  5. Candidate integer-interval *shape* for N(λ) from float brackets
     (YELLOW — not Arb/Rump certified)

Language: addresses counting gap; does NOT claim N(λ)=[0,0] or dual GREEN.

Usage:
  python route_A_arb_scaffold.py
  python route_A_arb_scaffold.py --N1 6 --N2 3 --N3 3 --neigs 6
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

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from route_A_counting import (  # noqa: E402
    KAPPA1,
    DEFAULT_Y,
    SMOKE_MESH,
    build_picard_core,
    float_dirichlet_eigs,
    float_neumann_eigs,
    candidate_brackets,
    diagnostic_glb_lower,
    truncation_constants_scaffold,
    default_status,
)


def float_to_mid_rad(
    A: np.ndarray,
    rel_rad: float = 1e-12,
    abs_rad: float = 1e-15,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Naive interval enclosure of a float matrix:

      A ∈ [mid − rad, mid + rad],  rad = rel*|A| + abs.

    Placeholder for true Arb assembly of hyperbolic weights.
    """
    mid = np.array(A, dtype=float, copy=True)
    rad = rel_rad * np.abs(mid) + abs_rad
    return mid, rad


def interval_matvec_residual(
    Qmid: np.ndarray,
    Qrad: np.ndarray,
    Mmid: np.ndarray,
    Mrad: np.ndarray,
    v: np.ndarray,
    lam: float,
) -> Dict[str, float]:
    """
    Enclose r = (Q − λ M) v componentwise with mid/rad arithmetic:

      r_mid = Qmid v − λ Mmid v
      r_rad = Qrad |v| + |λ| Mrad |v|
    """
    v = np.asarray(v, dtype=float).ravel()
    av = np.abs(v)
    r_mid = Qmid @ v - lam * (Mmid @ v)
    r_rad = Qrad @ av + abs(lam) * (Mrad @ av)
    return dict(
        residual_mid_norm=float(np.linalg.norm(r_mid)),
        residual_rad_norm=float(np.linalg.norm(r_rad)),
        residual_encl_norm=float(np.linalg.norm(np.abs(r_mid) + r_rad)),
        max_rad=float(np.max(r_rad)),
    )


def candidate_N_interval(
    lam_N: Sequence[float],
    lam_D: Sequence[float],
    lam: float,
) -> Dict[str, Any]:
    """
    Classical D–N counting *shape* on a fixed compact domain:

      # {k : λ_k^D ≤ λ}  ≤  N(λ)  ≤  # {k : λ_k^N ≤ λ}

    On the truncated core with free side BC this is only a candidate.
    Returns integer bounds as a YELLOW enclosure (not certified).
    """
    n_lo = sum(1 for x in lam_D if x <= lam + 1e-14)
    n_hi = sum(1 for x in lam_N if x <= lam + 1e-14)
    # ensure lo ≤ hi even with float noise
    if n_lo > n_hi:
        n_lo, n_hi = n_hi, n_lo
    return dict(
        lambda_target=float(lam),
        N_lower_candidate=int(n_lo),
        N_upper_candidate=int(n_hi),
        enclosure=[int(n_lo), int(n_hi)],
        status="YELLOW",
        note=(
            "Candidate integer interval from float D/N on K_Y only. "
            "Not Arb/Rump; truncation to Γ\\H³ not enclosed. "
            "Do not claim N(λ)=[0,0] for dual cert."
        ),
    )


def run_scaffold(
    N1: int,
    N2: int,
    N3: int,
    Y: float = DEFAULT_Y,
    neigs: int = 6,
    rel_rad: float = 1e-12,
    lambda_probe: float = 1.0,
) -> Dict[str, Any]:
    t0 = time.time()
    data = build_picard_core(N1, N2, N3, Y=Y, curved=True)
    Q = np.asarray(data["Q"], dtype=float)
    M = np.asarray(data["M"], dtype=float)
    top = data["top_face_dofs"]
    h_max = float(data["geo"]["hT"].max())

    lam_N = float_neumann_eigs(Q, M, neigs=neigs)
    lam_D = float_dirichlet_eigs(Q, M, top, neigs=neigs)
    brackets = candidate_brackets(lam_N, lam_D)

    Qmid, Qrad = float_to_mid_rad(Q, rel_rad=rel_rad)
    Mmid, Mrad = float_to_mid_rad(M, rel_rad=rel_rad)

    # residual of first positive Neumann mode if available
    from scipy.linalg import eigh

    w, vecs = eigh(Q, M, subset_by_index=[0, min(neigs, Q.shape[0] - 1)])
    # pick first mode with λ > 1e-8
    mode_idx = 0
    for i, wi in enumerate(w):
        if wi > 1e-8:
            mode_idx = i
            break
    v = vecs[:, mode_idx]
    resid = interval_matvec_residual(Qmid, Qrad, Mmid, Mrad, v, float(w[mode_idx]))

    glb_N = [diagnostic_glb_lower(float(x), h_max) for x in lam_N]
    glb_D = [diagnostic_glb_lower(float(x), h_max) for x in lam_D]

    # N(λ) candidates at λ=1 (FEM exclusion edge) and at first D mode
    N_at_1 = candidate_N_interval(lam_N, lam_D, lambda_probe)
    N_at_D1 = candidate_N_interval(lam_N, lam_D, float(lam_D[0]))

    trunc = truncation_constants_scaffold(Y=Y)
    status = default_status()
    # upgrade artificial_boundary note to reference this scaffold
    status.set(
        "artificial_boundary",
        False,
        "Scaffold mid/rad + trunc constants present (route_A_arb_scaffold / "
        "truncation_constants_scaffold). Still no proved Δ_trunc enclosure.",
    )
    status.set(
        "dn_guaranteed",
        False,
        "Float D/N + relative mid/rad on Q,M. Not Arb assembly; not Rump.",
    )
    status.set(
        "counting_enclosure",
        False,
        f"Candidate N({lambda_probe})={N_at_1['enclosure']} from float brackets "
        "(YELLOW). Not certified integer interval.",
    )

    out = dict(
        language=(
            "Route A Arb scaffold: float D/N + relative mid/rad + GLB sketch + "
            "candidate N(λ) shape + truncation formulae. "
            "Addresses counting gap; not certified."
        ),
        mesh=dict(N1=N1, N2=N2, N3=N3, Y=Y, n_dofs=int(Q.shape[0]), h_max=h_max),
        KAPPA1=float(KAPPA1),
        rel_rad=rel_rad,
        lam_N=[float(x) for x in lam_N],
        lam_D=[float(x) for x in lam_D],
        glb_sketch_N=glb_N,
        glb_sketch_D=glb_D,
        brackets=[(k, lo, hi) for k, lo, hi in brackets],
        interval_residual_mode=resid,
        N_lambda_candidates=dict(at_1=N_at_1, at_first_D=N_at_D1),
        truncation=trunc,
        checklist={k: status.items[k] for k in status.items},
        seconds=time.time() - t0,
        hard_map_note="counting_certified remains false",
    )
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--N1", type=int, default=SMOKE_MESH[0])
    p.add_argument("--N2", type=int, default=SMOKE_MESH[1])
    p.add_argument("--N3", type=int, default=SMOKE_MESH[2])
    p.add_argument("--Y", type=float, default=DEFAULT_Y)
    p.add_argument("--neigs", type=int, default=6)
    p.add_argument("--rel-rad", type=float, default=1e-12)
    p.add_argument("--lambda-probe", type=float, default=1.0)
    p.add_argument("--json-out", type=str, default="")
    args = p.parse_args(argv)

    print("=== Route A Arb / D–N scaffold (YELLOW counting) ===")
    out = run_scaffold(
        args.N1,
        args.N2,
        args.N3,
        Y=args.Y,
        neigs=args.neigs,
        rel_rad=args.rel_rad,
        lambda_probe=args.lambda_probe,
    )
    m = out["mesh"]
    print(
        f"mesh {m['N1']}x{m['N2']}x{m['N3']} Y={m['Y']}  "
        f"dofs={m['n_dofs']}  h_max={m['h_max']:.4e}"
    )
    print(f"{'k':>4}  {'λ^N':>12}  {'λ^D':>12}  {'GLB(N)':>12}")
    for (k, lo, hi), glb in zip(out["brackets"], out["glb_sketch_N"]):
        print(f"{k:4d}  {lo:12.6f}  {hi:12.6f}  {glb:12.6f}")
    print("\nN(λ) candidates:")
    for name, rec in out["N_lambda_candidates"].items():
        print(
            f"  {name}: enclosure={rec['enclosure']}  "
            f"status={rec['status']}"
        )
    print(
        f"\ninterval residual (first pos mode): "
        f"encl_norm={out['interval_residual_mode']['residual_encl_norm']:.3e}"
    )
    print(
        f"trunc Δ_model≈{out['truncation']['delta_trunc_model']:.3e}  "
        f"collar_Y≈{out['truncation']['collar_floor_Y']:.3e}"
    )
    print(f"({out['seconds']:.2f}s)  {out['language']}")

    path = Path(args.json_out) if args.json_out else _HERE / "route_A_arb_result.json"
    path.write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    print(f"wrote {path}")

    # infrastructure: Neumann ≤ Dirichlet, finite
    ok = all(lo <= hi + 1e-8 for _, lo, hi in out["brackets"])
    print("PASS infrastructure" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
