#!/usr/bin/env python3
"""
Rung 4 — dual-certification pipeline for Γ₀(2+i), N𝔭=5.

AgentReady §6 Rung 4 checklist orchestration:

  [x] FEM lower: Theorem 2 — no eigenvalue in (0,1) (re-verified by import/status)
  [x] Two-cusp Hejhal uses Rung 3 coupling (two_cusp_hejhal_N5)
  [ ] Overlap + counting with ε<0.1  — overlap computed; counting status explicit
  [x] F5 residue unit tests (residue_F5_tests)
  [x] Gluing of 6 copies exact (same)

Certification language (AgentReady §12–13):
  - FEM exclusion of (0,1) is **machine-certified** (independent_exclusion Thm 2).
  - Hejhal candidate + D(K) defect bound give an **existence interval under
    Assumptions H,A,S** when residual η is small enough — not a full §13 cert
    until counting + reproducibility suite pass.
  - This runner reports GREEN/YELLOW/RED per item; only emits "Rung4 certified"
    if every hard item is green (counting currently remains YELLOW unless Route A/B cert).

Usage:
  python rung4_N5_dual.py
  python rung4_N5_dual.py --M 48 --r-grid 6.0:12.0:25
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_ROOT / "independent_exclusion"))

from two_cusp_hejhal_N5 import (  # noqa: E402
    build_block_system, coupled_mid_rad, precondition_right, R_THEN, Y0_DEFAULT,
)
from defect_bound_arb import defect_to_lambda_error, ball_mid  # noqa: E402


# ---------------------------------------------------------------------------
# FEM side (Theorem 2)
# ---------------------------------------------------------------------------

def fem_N5_status() -> Dict[str, Any]:
    """
    Re-use certified Theorem 2: no Laplace eigenvalue in (0,1) on Γ₀(2+i).

    Does not re-run the full Rump cert (expensive); records provenance and
    optional light import checks.
    """
    out = dict(
        level="(2+i)",
        NP=5,
        index=6,
        claim="no eigenvalue in (0,1)",
        source="independent_exclusion/PROOF.md Theorem 2; m3p_certify.py",
        frozen=dict(
            mesh="8x4x3 24-split",
            n_dofs=28400,
            Y=1.25,
            params="(θ,θ2,α,θ4,ρ̃)=(0.6,0.9,0.2,0.85,9.0)",
            nu_star=1.02,
            float_pencil_min_eig=1.44,
            certified_date="2026-07-10",
        ),
        # Interval for dual loop: rigorous lower endpoint from exclusion
        fem_lower_open=1.0,  # spectrum in (0,1) empty ⇒ λ1 ≥ 1
        fem_upper_engineering=None,  # no certified FEM upper near 45 yet
        language=(
            "FEM certifies exclusion of (0,1) only. Any upper bound near the "
            "first eigenvalue remains an engineering target until a conforming "
            "Rayleigh certificate exists."
        ),
        status="GREEN",
    )
    # Light re-verify residue + gluing used by the FEM multi-copy path
    try:
        import congruence_prototype as cp
        pts, glue, cc = cp.build_gluing("(2+i)")
        assert len(pts) == 6 and cc.count(0) == 1 and cc.count(1) == 5
        out["gluing_reverified"] = True
        out["glue_generators"] = {k: v for k, v in glue.items()}
    except Exception as e:
        out["gluing_reverified"] = False
        out["gluing_error"] = str(e)
        out["status"] = "RED"
    return out


# ---------------------------------------------------------------------------
# Hejhal residual scan (two-cusp, Rung 3 operator)
# ---------------------------------------------------------------------------

def smallest_singular(M: np.ndarray) -> float:
    s = np.linalg.svd(M, compute_uv=False)
    return float(s[-1]) if len(s) else float("inf")


def _sigma_at(M: int, r: float, Y0: float, theta: float) -> Tuple[float, float, float]:
    """Return (sigma_min, delta, tau) residual proxies at r."""
    sysb = build_block_system(M, float(r), Y0, theta=theta)
    mid, rad = coupled_mid_rad(sysb)
    mid_p, _, _ = precondition_right(mid, rad, sysb.w_inf, sysb.w_0)
    col = np.maximum(np.linalg.norm(mid_p, axis=0), 1e-300)
    G = mid_p / col[None, :]
    # relative residual = σ_min / σ_max (scale-free)
    svals = np.linalg.svd(G, compute_uv=False)
    sig = float(svals[-1])
    sig_max = float(svals[0]) if len(svals) else 1.0
    rel = sig / max(sig_max, 1e-300)
    # Use relative residual as automorphy defect; τ scaled by elliptic factor
    delta = rel
    tau = rel * math.sqrt(1.0 + float(r) ** 2)
    return sig, delta, tau


def hejhal_residual_scan(
    M: int = 48,
    Y0: float = 1.25,
    r_min: float = 6.0,
    r_max: float = 14.0,
    n_r: int = 33,
    theta: float = 0.5,
    refine: bool = True,
) -> Dict[str, Any]:
    """
    Scan r ↦ relative σ_min for the two-cusp collocation operator, then
    optionally golden-section refine. Candidate spectral parameter under H,A,S.
    """
    rs = np.linspace(r_min, r_max, n_r)
    rows = []
    best = None
    t0 = time.time()
    for r in rs:
        sig, delta, tau = _sigma_at(M, float(r), Y0, theta)
        rec = dict(
            r=float(r),
            sigma_min=sig,
            delta=delta,
            tau=tau,
            eta=delta + tau,
        )
        rows.append(rec)
        if best is None or rec["eta"] < best["eta"]:
            best = rec
    # Local refine around best
    if refine and best is not None:
        a = max(r_min, best["r"] - (r_max - r_min) / max(n_r, 2))
        b = min(r_max, best["r"] + (r_max - r_min) / max(n_r, 2))
        phi = (math.sqrt(5) - 1) / 2
        for _ in range(12):
            r1 = b - phi * (b - a)
            r2 = a + phi * (b - a)
            _, d1, t1 = _sigma_at(M, r1, Y0, theta)
            _, d2, t2 = _sigma_at(M, r2, Y0, theta)
            e1, e2 = d1 + t1, d2 + t2
            if e1 < e2:
                b, best = r2, dict(r=r1, sigma_min=0.0, delta=d1, tau=t1, eta=e1)
                best["sigma_min"] = _sigma_at(M, r1, Y0, theta)[0]
            else:
                a, best = r1, dict(r=r2, sigma_min=0.0, delta=d2, tau=t2, eta=e2)
                best["sigma_min"] = _sigma_at(M, r2, Y0, theta)[0]
    elapsed = time.time() - t0
    return dict(
        M=M,
        Y0=Y0,
        r_grid=[float(x) for x in rs],
        rows=rows,
        best=best,
        seconds=elapsed,
        language=(
            "relative sigma_min is a discrete collocation residual proxy, "
            "not a certified Maass residual; D(K) applied under H,A,S."
        ),
    )


def apply_defect_bound(best: Dict[str, Any], Y: float = 1.25, Y0: float = 1.25) -> Dict[str, Any]:
    """Apply Theorem D(K) plug-in to (δ,τ) at best r."""
    r = best["r"]
    delta = best["delta"]
    tau = best["tau"]
    out = defect_to_lambda_error(
        delta, tau, Y=Y, Y0=Y0, r=r, field="i", M=400, theta=0.5, C_H=1.0,
        sharp_geom=True, U_norm=1.0, classical_CK=False,
    )
    lam = 1.0 + r * r
    dlam = float(ball_mid(out["lambda_error_bound"]))
    dr = dlam / max(2.0 * r, 1e-12)
    # Practical width target for Rung 4: need dlam < eps_tol (default 0.1)
    # Collar floor C2 e^{-2π Y0}:
    C2 = float(ball_mid(out["C2"]))
    collar = C2 * math.exp(-2 * math.pi * Y0)
    return dict(
        r_candidate=r,
        lambda_candidate=lam,
        delta=delta,
        tau=tau,
        eta=delta + tau,
        lambda_error_bound=dlam,
        r_error_bound=dr,
        collar_floor=collar,
        L2_error_bound=float(ball_mid(out["L2_error_bound"])) if out.get("L2_error_bound") is not None else None,
        eta0=float(ball_mid(out["eta0"])),
        eta_le_eta0=bool(out.get("eta_le_eta0")),
        eta_over_eta0=float(out.get("eta_over_eta0", float("nan"))),
        eta_orders_gap=float(out.get("eta_orders_of_magnitude_gap", float("nan"))),
        C1=float(ball_mid(out["C1"])),
        C2=C2,
        sharp_geom=bool(out.get("sharp_geom", True)),

        assumptions=out.get("assumptions"),
        hejhal_interval_lambda=(lam - dlam, lam + dlam),
        hejhal_interval_r=(r - dr, r + dr),
        Y0_used=Y0,
    )


# ---------------------------------------------------------------------------
# Dual overlap + counting status
# ---------------------------------------------------------------------------

def dual_overlap(
    fem: Dict[str, Any],
    hej: Dict[str, Any],
    eps_tol: float = 0.1,
) -> Dict[str, Any]:
    """
    Dual loop:
      FEM: spectrum ∩ (0,1)=∅  ⇒  λ1 ≥ 1  (certified)
      Hejhal+D(K): λ ∈ [λ̃ - dλ, λ̃ + dλ] under H,A,S

    Rung 4 wants diam < eps_tol, disjoint from (0,1), + counting.
    """
    lo_h, hi_h = hej["hejhal_interval_lambda"]
    # Clamp nonsense negative lows from huge defect bounds for reporting
    lo_h_rep, hi_h_rep = lo_h, hi_h
    fem_lo = fem["fem_lower_open"]
    ov_lo = max(fem_lo, lo_h)
    ov_hi = hi_h
    nonempty = ov_lo <= ov_hi
    hej_width = max(0.0, hi_h - lo_h)
    width = max(0.0, ov_hi - ov_lo) if nonempty else float("inf")
    disjoint_from_unit = lo_h >= 1.0 - 1e-15
    # width_ok: the *Hejhal* existence interval is narrower than tol
    width_ok = hej_width < eps_tol and hej_width > 0
    counting = dict(
        status="YELLOW",
        note=(
            "Route A prototype (route_A_counting.py) addresses the counting "
            "gap but is not yet a certified N(λ)=[0,0] on (1,λ1)."
        ),
    )
    # Practical dual green: narrow Hejhal interval above 1, overlaps [1,∞)
    practical_green = (
        nonempty and disjoint_from_unit and width_ok and hej.get("eta_le_eta0", False)
    )
    return dict(
        fem_lower=fem_lo,
        hejhal_interval=(lo_h_rep, hi_h_rep),
        hejhal_width=hej_width,
        overlap=(ov_lo, ov_hi) if nonempty else None,
        nonempty=nonempty,
        width=width,
        width_tol=eps_tol,
        width_ok=width_ok,
        disjoint_from_unit_interval=disjoint_from_unit,
        eta_le_eta0=hej.get("eta_le_eta0"),
        collar_floor=hej.get("collar_floor"),
        counting=counting,
        dual_interval_status=(
            "GREEN" if practical_green else
            "YELLOW" if nonempty else
            "RED"
        ),
        language=(
            "FEM gives certified λ≥1. Hejhal+D(K) gives a conditional existence "
            "interval under H,A,S. Full §13 first-eigenvalue cert needs "
            "counting + η≤η0 + width<tol + reproducibility."
        ),
    )


def counting_float_probe() -> Dict[str, Any]:
    """Optional float Route A smoke — not certified."""
    try:
        from route_A_counting import run_smoke  # type: ignore
        return dict(status="YELLOW", smoke="available", note="see route_A_status.md")
    except Exception:
        return dict(
            status="YELLOW",
            note="route_A_counting present; run `python route_A_counting.py --smoke`",
        )


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def demo_tight_eta(
    r: float = 8.0,
    Y0_defect: float = 1.5,
    eta: float = 1e-10,
    eps_tol: float = 0.1,
) -> Dict[str, Any]:
    """
    Pipeline self-test: if a Hejhal run *had* residual η this small,
    D(K) would produce width < eps_tol. Used to validate dual arithmetic.
    """
    fem = fem_N5_status()
    delta = 0.5 * eta
    tau = 0.5 * eta
    best = dict(r=r, sigma_min=eta, delta=delta, tau=tau, eta=eta)
    hej = apply_defect_bound(best, Y=1.25, Y0=Y0_defect)
    dual = dual_overlap(fem, hej, eps_tol=eps_tol)
    return dict(demo=True, hej=hej, dual=dual, target_eta=eta, Y0_defect=Y0_defect)


def run_rung4(
    M: int = 48,
    Y0: float = 1.25,
    r_min: float = 6.0,
    r_max: float = 14.0,
    n_r: int = 25,
    eps_tol: float = 0.1,
    Y0_defect: Optional[float] = None,
    run_demo: bool = True,
) -> Dict[str, Any]:
    """
    Y0: collocation height for residual scan.
    Y0_defect: height fed to D(K) collar term (default max(Y0, 1.25)).
    """
    if Y0_defect is None:
        Y0_defect = max(Y0, 1.25)
    print("=" * 68)
    print("Rung 4 — Dual certification pipeline for Γ₀(2+i), N𝔭=5")
    print("=" * 68)

    # 1) F5 + gluing
    print("\n[1] Residue F5 + exact 6-copy gluing")
    import subprocess
    rc = subprocess.call([sys.executable, "-u", str(_HERE / "residue_F5_tests.py")])
    f5_ok = rc == 0
    print(f"  residue_F5_tests exit={rc}  {'PASS' if f5_ok else 'FAIL'}")

    # 2) FEM
    print("\n[2] FEM Theorem 2 reuse")
    fem = fem_N5_status()
    print(f"  status={fem['status']}  claim={fem['claim']}")
    print(f"  gluing_reverified={fem.get('gluing_reverified')}")

    # 3) Hejhal residual scan + D(K)
    print("\n[3] Two-cusp Hejhal residual scan (Rung 3 operator)")
    print(f"  M={M}  r∈[{r_min},{r_max}]  n_r={n_r}  Y0_colloc={Y0}")
    scan = hejhal_residual_scan(M=M, Y0=Y0, r_min=r_min, r_max=r_max, n_r=n_r)
    best = scan["best"]
    print(f"  best r≈{best['r']:.6f}  sigma_min={best['sigma_min']:.6e}  "
          f"eta={best['eta']:.6e}  ({scan['seconds']:.1f}s)")
    print("\n[4] Theorem D(K) defect bound")
    print(f"  Y0_defect={Y0_defect} (collar floor C2 e^{{-2π Y0}})")
    hej = apply_defect_bound(best, Y=1.25, Y0=Y0_defect)
    print(f"  λ̃={hej['lambda_candidate']:.6f}  dλ≤{hej['lambda_error_bound']:.6e}")
    print(f"  collar_floor≈{hej['collar_floor']:.6e}  C1≈{hej['C1']:.6e}")
    print(f"  Hejhal λ-interval ≈ [{hej['hejhal_interval_lambda'][0]:.6f}, "
          f"{hej['hejhal_interval_lambda'][1]:.6f}]")
    print(f"  eta≤eta0? {hej['eta_le_eta0']}  (eta0≈{hej['eta0']:.3e})")
    # Quantitative gap (honest blocker)
    gap = hej["eta"] / max(hej["eta0"], 1e-300)
    orders = math.log10(gap) if gap > 0 else 0.0
    print(f"  eta/eta0 ≈ {gap:.3e}  (~{orders:.1f} orders of magnitude gap)")
    print(f"  hard map must stay: eta_le_eta0=false until residual closes this gap")

    # 5) Dual overlap
    print("\n[5] Dual overlap")
    dual = dual_overlap(fem, hej, eps_tol=eps_tol)
    print(f"  dual_interval_status={dual['dual_interval_status']}")
    print(f"  nonempty={dual['nonempty']} width={dual['width']:.6e} "
          f"width_ok={dual['width_ok']} (tol={eps_tol})")
    print(f"  counting={dual['counting']['status']}: {dual['counting']['note'][:70]}...")

    # 6) Overall
    hard = {
        "F5_gluing": f5_ok,
        "FEM_exclusion_01": fem["status"] == "GREEN",
        "two_cusp_hejhal_scan": best is not None,
        "defect_bound_applied": True,
        "eta_le_eta0": bool(hej["eta_le_eta0"]),
        "overlap_nonempty_disjoint_unit": bool(
            dual["nonempty"] and dual["disjoint_from_unit_interval"]
        ),
        "width_lt_tol": bool(dual["width_ok"]),
        "counting_certified": dual["counting"]["status"] == "GREEN",
    }
    # AgentReady full stop needs counting GREEN and width + eta
    rung4_certified = all(
        [
            hard["F5_gluing"],
            hard["FEM_exclusion_01"],
            hard["eta_le_eta0"],
            hard["overlap_nonempty_disjoint_unit"],
            hard["width_lt_tol"],
            hard["counting_certified"],
        ]
    )
    # Partial success: pipeline green except counting / width
    pipeline_ok = hard["F5_gluing"] and hard["FEM_exclusion_01"] and hard["two_cusp_hejhal_scan"]

    print("\n" + "=" * 68)
    print("Rung 4 hard items:")
    for k, v in hard.items():
        print(f"  [{'YES' if v else ' NO'}] {k}")
    print(f"\n  pipeline_ok (infra)     = {pipeline_ok}")
    print(f"  Rung4 §13 CERTIFIED     = {rung4_certified}")
    if not rung4_certified:
        print("  Reason: counting not certified and/or defect η>η0 / width≥tol.")
        print("  Dual loop still useful under H,A,S with engineering Hejhal residual.")

    demo = None
    if run_demo:
        print("\n[6] Dual-loop arithmetic demo (synthetic tight η)")
        # η needed for dλ < 0.1 at Y0=1.5: roughly (0.1 - C2 e^{-2π*1.5})/C1
        demo = demo_tight_eta(r=max(best["r"], 7.0), Y0_defect=1.5, eta=1e-14, eps_tol=eps_tol)
        print(f"  synthetic eta=1e-14  dλ≤{demo['hej']['lambda_error_bound']:.6e}  "
              f"width_ok={demo['dual']['width_ok']}  "
              f"dual={demo['dual']['dual_interval_status']}  "
              f"eta≤η0={demo['hej']['eta_le_eta0']}")
        print(f"  (shows D(K)+overlap arithmetic; not a physical residual)")

    result = dict(
        hard=hard,
        pipeline_ok=pipeline_ok,
        rung4_certified=rung4_certified,
        fem=fem,
        hejhal_scan_best=best,
        hejhal_defect=hej,
        dual=dual,
        dual_demo_tight_eta=demo,
        assumptions=["H", "A", "S"],
        level="(2+i)",
        NP=5,
        blockers=[
            k for k, v in hard.items() if not v
        ],
    )

    out_json = _HERE / "rung4_result.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, default=str)
    print(f"\n  wrote {out_json}")
    return result


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Rung 4 N=5 dual pipeline")
    p.add_argument("--M", type=int, default=48)
    p.add_argument("--Y0", type=float, default=1.25, help="collocation height")
    p.add_argument("--Y0-defect", type=float, default=None,
                   help="D(K) collar height (default max(Y0,1.25))")
    p.add_argument("--r-grid", type=str, default="6.0:14.0:25",
                   help="r_min:r_max:n_points")
    p.add_argument("--eps", type=float, default=0.1, help="width tolerance")
    args = p.parse_args(argv)
    rmin, rmax, nr = args.r_grid.split(":")
    res = run_rung4(
        M=args.M,
        Y0=args.Y0,
        r_min=float(rmin),
        r_max=float(rmax),
        n_r=int(nr),
        eps_tol=args.eps,
        Y0_defect=args.Y0_defect,
    )
    # Exit 0 if infrastructure ok; exit 2 if full cert; exit 1 on hard fail
    if res["rung4_certified"]:
        return 0
    if res["pipeline_ok"]:
        return 0  # partial success still exit 0 for CI of infrastructure
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
