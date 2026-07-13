#!/usr/bin/env python3
"""
CR guaranteed-lower-bound (GLB) sketch on V_h^{P1,per} first positive mode.

Uses the classical CR GLB shape (Carstensen–Gedicke / Liu style):

  λ ≥ λ_h / (1 + κ₁² h_max² λ_h)

with κ₁ = KAPPA1 from independent_exclusion.cr_prototype (CZZ/CP22),
and h_max = max tet diameter from the *reference* cell mesh (same scale
as multi-copy assembly).

Language (honest)
-----------------
- Applied to multi-copy pairing-conforming CR eigenvalues from
  v_h_p1_per_spectrum (float eigsh).
- GLB formula is the *same shape* as level-1 exclusion; applying it on
  Γ₀ multi-copy + weighted hyperbolic form is engineering until the
  full CR theory hypotheses are re-verified for this space (as m3p does
  for exclusion).
- Does NOT claim certified continuum λ₁ or dual GREEN.
- Hard map unchanged.

Usage:
  python cr_glb_p1_per.py
  python cr_glb_p1_per.py --N1 4 --N2 2 --N3 2 --neigs 6
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
_IE = _ROOT / "independent_exclusion"
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_IE))

from cr_prototype import KAPPA1  # noqa: E402
from congruence_prototype import build_reference, ref_geometry, set_level  # noqa: E402
from v_h_p1_per_spectrum import assemble_periodic, float_neumann_eigs, float_dirichlet_eigs  # noqa: E402


def reference_h_max(N1: int, N2: int, N3: int) -> float:
    """Max Euclidean edge length over reference-cell tets (diameter proxy)."""
    X, tets, _ = build_reference(N1, N2, N3)
    h_max = 0.0
    for tet in tets:
        P = X[list(tet)]
        for i in range(4):
            for j in range(i + 1, 4):
                d = float(np.linalg.norm(P[i] - P[j]))
                h_max = max(h_max, d)
    return h_max


def cr_glb(lam_h: float, h_max: float, kappa: float = float(KAPPA1)) -> float:
    """λ ≥ λ_h / (1 + κ² h² λ_h). Returns 0 if λ_h ≤ 0."""
    if lam_h <= 0:
        return 0.0
    den = 1.0 + (kappa ** 2) * (h_max ** 2) * lam_h
    return float(lam_h / den)


def run(
    N1: int = 4,
    N2: int = 2,
    N3: int = 2,
    neigs: int = 6,
) -> Dict[str, Any]:
    t0 = time.time()
    set_level("(2+i)")
    h_max = reference_h_max(N1, N2, N3)
    kappa = float(KAPPA1)

    data = assemble_periodic(N1, N2, N3, verbose=True)
    Q, M = data["Q"], data["M"]
    top = data["top_dofs"]
    lam_N = float_neumann_eigs(Q, M, neigs=neigs)
    lam_D = float_dirichlet_eigs(Q, M, top, neigs=neigs)

    rows = []
    for k in range(min(len(lam_N), len(lam_D))):
        ln, ld = float(lam_N[k]), float(lam_D[k])
        rows.append(
            dict(
                k=k + 1,
                lam_N=ln,
                lam_D=ld,
                glb_N=cr_glb(ln, h_max, kappa),
                glb_D=cr_glb(ld, h_max, kappa),
            )
        )

    # first positive Neumann
    first_pos = None
    first_pos_glb = None
    for row in rows:
        if row["lam_N"] > 1e-6:
            first_pos = row["lam_N"]
            first_pos_glb = row["glb_N"]
            break

    out = dict(
        language=(
            "CR GLB sketch λ ≥ λ_h/(1+κ₁² h_max² λ_h) on V_h^{P1,per} float "
            "spectrum. Engineering — full multi-copy CR theory not re-certified "
            "here. Hard map unchanged."
        ),
        KAPPA1=kappa,
        h_max=h_max,
        mesh=dict(N1=N1, N2=N2, N3=N3),
        n_dofs=data["n_dofs"],
        rows=rows,
        first_positive_Neumann=first_pos,
        first_positive_GLB=first_pos_glb,
        glb_factor=float(first_pos_glb / first_pos) if first_pos and first_pos_glb else None,
        seconds=time.time() - t0,
        warnings=[
            "GLB requires CR interpolation hypotheses on this mesh/space.",
            "Not a dual upper bound; first positive N is upper-flavor discrete.",
            "Do not set counting_certified or rung4_certified from this alone.",
        ],
    )
    print(f"\n=== CR GLB on V_h^{{P1,per}} ===")
    print(f"κ₁={kappa:.6f}  h_max={h_max:.6f}")
    print(f"{'k':>4}  {'λ^N':>12}  {'GLB(N)':>12}  {'λ^D':>12}  {'GLB(D)':>12}")
    for row in rows:
        print(
            f"{row['k']:4d}  {row['lam_N']:12.6f}  {row['glb_N']:12.6f}  "
            f"{row['lam_D']:12.6f}  {row['glb_D']:12.6f}"
        )
    print(
        f"\nfirst positive Neumann ≈ {first_pos:.6f}  "
        f"GLB ≈ {first_pos_glb:.6f}  "
        f"(factor {out['glb_factor']:.4f})"
    )
    print(f"({out['seconds']:.2f}s)")
    for w in out["warnings"]:
        print(f"  WARN: {w}")
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--N1", type=int, default=4)
    p.add_argument("--N2", type=int, default=2)
    p.add_argument("--N3", type=int, default=2)
    p.add_argument("--neigs", type=int, default=6)
    p.add_argument("--json-out", type=str, default="")
    args = p.parse_args(argv)

    out = run(args.N1, args.N2, args.N3, neigs=args.neigs)
    path = Path(args.json_out) if args.json_out else _HERE / "cr_glb_p1_per_result.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {path}")
    ok = (
        out["first_positive_Neumann"] is not None
        and out["first_positive_GLB"] is not None
        and out["first_positive_GLB"] <= out["first_positive_Neumann"] + 1e-12
        and out["first_positive_GLB"] > 0
    )
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
