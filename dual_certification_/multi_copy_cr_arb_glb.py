#!/usr/bin/env python3
"""
Multi-copy CR Arb geometry checks + GLB on V_h^{P1,per} (Γ₀(2+i)).

What is certified here
----------------------
- Reference-cell tet volumes > 0 (Arb, m3_certify.tet_arb_data) for every
  tet of the reference mesh used in multi-copy assembly.
- h_max upper bound from Arb edge lengths on the reference cell.
- Float spectrum + GLB sketch λ ≥ λ_h/(1+κ₁² h_max² λ_h) using that h_max.

What is NOT certified
---------------------
- Full multi-copy weighted CR interpolation theory (κ₁ on glued space).
- Continuum lower bound λ₁ ≥ GLB (needs complete CR hypotheses + BC).
- Dual interval [GLB, λ_h] as certified λ₁.

See dual_glb_language.md. Hard map unchanged.

Usage:
  python multi_copy_cr_arb_glb.py
  python multi_copy_cr_arb_glb.py --N1 4 --N2 2 --N3 2
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

_HERE = Path(__file__).resolve().parent
_IE = _HERE.parent / "independent_exclusion"
sys.path.insert(0, str(_IE))
sys.path.insert(0, str(_HERE))

import flint  # noqa: E402
from flint import arb  # noqa: E402

from congruence_prototype import (  # noqa: E402
    build_reference,
    set_level,
    INDEX,
    LEVEL,
    NP,
)
from m3_certify import tet_arb_data, upper, mid_rad  # noqa: E402
from cr_prototype import KAPPA1  # noqa: E402
from v_h_p1_per_spectrum import (  # noqa: E402
    assemble_periodic,
    float_neumann_eigs,
    float_dirichlet_eigs,
)
from cr_glb_p1_per import cr_glb  # noqa: E402

flint.ctx.prec = 128


def arb_reference_geometry(N1: int, N2: int, N3: int) -> Dict[str, Any]:
    """Arb-certify all reference tets: vol>0, track h_max upper.

    Orients each tet so det>0 (same as congruence ref_geometry flip).
    """
    X, tets, btri = build_reference(N1, N2, N3)
    tets = np.array(tets, dtype=int, copy=True)
    nt = len(tets)
    n_pos = 0
    n_flipped = 0
    h_max_u = 0.0
    h_list = []
    vol_mids = []
    for e in range(nt):
        tet = list(tets[e])
        P = X[tet]
        # float orientation check
        A = np.hstack([np.ones((4, 1)), P])
        det = np.linalg.det(A)
        if det < 0:
            tet = [tet[1], tet[0], tet[2], tet[3]]
            tets[e] = tet
            P = X[tet]
            n_flipped += 1
        try:
            vol_a, grads, hT = tet_arb_data(P)
            pos = bool(vol_a > 0)
        except AssertionError:
            # last resort: absolute volume via edge diam only
            pos = False
            vol_a = arb(0)
            hT = arb(0)
            for i in range(4):
                for j in range(i + 1, 4):
                    d = sum(
                        (arb(float(P[i, k])) - arb(float(P[j, k]))) ** 2
                        for k in range(3)
                    )
                    from m3_certify import amax as _amax

                    hT = _amax(hT, d.sqrt()) if float(hT.mid()) > 0 else d.sqrt()
        if pos:
            n_pos += 1
        hu = upper(hT) if hasattr(hT, "mid") else float(hT)
        h_max_u = max(h_max_u, hu)
        h_list.append(hu)
        if pos:
            vm, _ = mid_rad(vol_a)
            vol_mids.append(vm)
    return dict(
        n_tets=nt,
        n_vol_positive=n_pos,
        n_flipped=n_flipped,
        all_vol_positive=n_pos == nt,
        h_max_upper=h_max_u,
        h_max_mid=float(np.median(h_list)) if h_list else 0.0,
        vol_sum_mid=float(sum(vol_mids)),
        n_boundary_tags={k: len(v) for k, v in btri.items()},
    )


def run(N1: int = 4, N2: int = 2, N3: int = 2, neigs: int = 6) -> Dict[str, Any]:
    t0 = time.time()
    set_level("(2+i)")
    print(
        f"=== Multi-copy CR Arb + GLB  Γ₀{LEVEL}  mesh {N1}x{N2}x{N3} ==="
    )

    # 1) Arb geometry on reference cell
    geo = arb_reference_geometry(N1, N2, N3)
    print(
        f"reference tets: vol>0 = {geo['n_vol_positive']}/{geo['n_tets']}  "
        f"h_max^U={geo['h_max_upper']:.6f}"
    )

    # 2) Multi-copy float assembly + spectrum
    data = assemble_periodic(N1, N2, N3, verbose=True)
    Q, M = data["Q"], data["M"]
    top = data["top_dofs"]
    lam_N = float_neumann_eigs(Q, M, neigs=neigs)
    lam_D = float_dirichlet_eigs(Q, M, top, neigs=neigs)

    kappa = float(KAPPA1)
    h_u = float(geo["h_max_upper"])
    rows = []
    for k in range(min(len(lam_N), len(lam_D))):
        ln, ld = float(lam_N[k]), float(lam_D[k])
        rows.append(
            dict(
                k=k + 1,
                lam_N=ln,
                lam_D=ld,
                glb_N=cr_glb(ln, h_u, kappa),
                glb_D=cr_glb(ld, h_u, kappa),
            )
        )
    first_pos = next((r for r in rows if r["lam_N"] > 1e-6), None)

    # 3) Checklist for multi-copy CR theory (honest YELLOW)
    checks = dict(
        ref_all_vol_positive=geo["all_vol_positive"],
        multi_copy_Q1_near_zero=data["Q1_norm"] < 1e-8,
        cusp_areas_t_inf=abs(data["t_inf_sum"] - 0.5) < 1e-6,
        cusp_areas_t_0=abs(data["t_0_sum"] - NP / 2) < 1e-6,
        first_pos_exists=first_pos is not None,
        glb_below_rayleigh=(
            first_pos is not None
            and first_pos["glb_N"] <= first_pos["lam_N"] + 1e-12
            and first_pos["glb_N"] > 0
        ),
        # theory gaps
        multi_copy_kappa_certified=False,  # not done
        continuum_glb_certified=False,  # not done
    )

    out = dict(
        language=(
            "Arb-certified reference tet volumes + h_max^U; float multi-copy "
            "V_h^{P1,per} spectrum; GLB sketch with that h_max. "
            "Multi-copy CR κ₁ theory and continuum GLB NOT certified. "
            "See dual_glb_language.md."
        ),
        level=LEVEL,
        NP=NP,
        INDEX=INDEX,
        mesh=dict(N1=N1, N2=N2, N3=N3),
        n_dofs=data["n_dofs"],
        n_copies=INDEX,
        reference_arb=geo,
        KAPPA1=kappa,
        h_max_upper_arb=h_u,
        rows=rows,
        first_positive=first_pos,
        checks=checks,
        theory_status=dict(
            ref_geometry="GREEN (Arb vol>0, h_max^U)",
            multi_copy_assembly="GREEN float (Q1, cusp areas)",
            multi_copy_CR_interpolation="YELLOW — κ₁ not re-proved on glued space",
            continuum_GLB="RED — do not claim λ₁ ≥ GLB certified",
            dual_interval="RED — only FEM λ₁≥1 is certified lower",
        ),
        seconds=time.time() - t0,
        hard_map_note="rung4_certified remains false",
    )

    print(f"\n{'k':>4}  {'λ^N':>12}  {'GLB(N)':>12}  {'λ^D':>12}")
    for r in rows:
        print(
            f"{r['k']:4d}  {r['lam_N']:12.6f}  {r['glb_N']:12.6f}  "
            f"{r['lam_D']:12.6f}"
        )
    if first_pos:
        print(
            f"\nfirst pos: λ_h={first_pos['lam_N']:.6f}  "
            f"GLB={first_pos['glb_N']:.6f}  "
            f"(factor {first_pos['glb_N']/first_pos['lam_N']:.4f})"
        )
    print("\nchecks:", {k: v for k, v in checks.items()})
    print("theory:", out["theory_status"])
    print(f"({out['seconds']:.2f}s)")
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
    path = (
        Path(args.json_out)
        if args.json_out
        else _HERE / "multi_copy_cr_arb_glb_result.json"
    )
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {path}")
    ok = (
        out["checks"]["ref_all_vol_positive"]
        and out["checks"]["multi_copy_Q1_near_zero"]
        and out["checks"]["glb_below_rayleigh"]
    )
    print("PASS infrastructure" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
