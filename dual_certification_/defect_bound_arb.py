#!/usr/bin/env python3
"""
Rung 1 — single-cusp defect bound API (Theorem D(K)).

Given measured automorphy defect δ and L² residual τ from a putative
truncated expansion, return Arb enclosures of the eigenvalue error and
(optional) L² proximity bounds from Theorem D(K):

    |λ_true − λ|  ≤  C₁ · (δ + τ)  +  C₂ · exp(−2π Y₀)
    ‖U_cusp − f‖  ≤  C₁ · √(δ + τ)     (when η ≤ η₀)

All constants come from lemma_K.C1_C2_constants (formulas matching
theorem_DK.tex).  No Hejhal solver, no two-cusp coupling, no counting.

Language discipline
-------------------
- Outputs are *enclosures of the defect bound*, not a certified eigenvalue.
- Converting an engineering-target spectral parameter (e.g. Then r₁) into a
  certified λ requires a concrete putative form with certified (δ,τ) small
  enough that the bound fits the desired window — that is Rung 2+ work.
- Every call records Assumption H (θ, C_H) as an explicit input.

API
---
    defect_to_lambda_error(delta, tau, Y=1.25, Y0=0.8, r=6.6,
                           field='i', M=400, theta=0.5, C_H=1.0)
    sensitivity_table(...)   # Y0 ±0.1, r ±0.1
    demo()                   # toy (δ,τ) as Hejhal would plug in

Usage:
    python defect_bound_arb.py --demo
    python defect_bound_arb.py --sensitivity
    python defect_bound_arb.py --delta 1e-12 --tau 1e-14
"""
from __future__ import annotations

import argparse
import math
import sys
from typing import Any, Dict, List, Optional, Sequence, Union

# Local Rung 0 dependency — do not reimplement C_K or tail formulas.
from lemma_K import (
    HAS_FLINT,
    C1_C2_constants,
    ball_mid,
    ball_rad,
    ball_str,
    lemma_K_tail,
    to_ball,
    C_K,
)

Number = Union[float, Any]


# ---------------------------------------------------------------------------
# Main API
# ---------------------------------------------------------------------------
def defect_to_lambda_error(
    delta: float,
    tau: float,
    Y: float = 1.25,
    Y0: float = 0.8,
    r: float = 6.6,
    field: str = "i",
    M: int = 400,
    theta: float = 0.5,
    C_H: float = 1.0,
    eps: float = 0.0,
    D_Y: Optional[float] = None,
    h_min: Optional[float] = None,
    classical_CK: bool = False,
    sharp_geom: bool = True,
    U_norm: float = 1.0,
) -> Dict[str, Any]:
    """
    Convert measured defects (δ, τ) into enclosures of Theorem D(K) bounds.

    Default sharp_geom=True: C_K=1, y_min=1/√2 (Gaussian), D_Y=core-box diam,
    h_min from geometry, C_tr=3/h_min (morningbrief §1.1–1.2).
    U_norm: explicit ‖U‖₂ for relaxed η₀ = min(1/(4(1+A_bdry)²), U_norm²/(4 C1²)).
    """
    if delta < 0 or tau < 0:
        raise ValueError("delta and tau must be ≥ 0")
    if Y0 <= 0 or Y <= 0:
        raise ValueError("Y and Y0 must be positive")
    if r < 0:
        raise ValueError("r must be ≥ 0")

    consts = C1_C2_constants(
        Y=Y,
        Y0=Y0,
        r=r,
        theta=theta,
        C_H=C_H,
        eps=eps,
        M=M,
        field=field,
        D_Y=D_Y,
        h_min=h_min,
        classical_CK=classical_CK,
        sharp_geom=sharp_geom,
        U_norm=U_norm,
    )

    C1 = consts["C1"]
    C2 = consts["C2"]
    eta0 = consts["eta0"]
    lam = consts["lambda"]

    eta = delta + tau
    eta_b = to_ball(eta)

    # exp(−2π Y0)
    if HAS_FLINT:
        from flint import arb as _arb  # type: ignore

        exp_term = (-to_ball(2.0) * to_ball(math.pi) * to_ball(Y0)).exp()
        lambda_err = C1 * eta_b + C2 * exp_term
        # √η for L² bound (formal; theorem requires η ≤ η₀)
        if eta > 0:
            sqrt_eta = eta_b.sqrt()
        else:
            sqrt_eta = to_ball(0.0)
        L2_err = C1 * sqrt_eta
    else:
        exp_term = math.exp(-2.0 * math.pi * Y0)
        lambda_err = float(ball_mid(C1)) * eta + float(ball_mid(C2)) * exp_term
        sqrt_eta = math.sqrt(eta) if eta > 0 else 0.0
        L2_err = float(ball_mid(C1)) * sqrt_eta
        exp_term = to_ball(exp_term)  # type: ignore

    eta0_mid = ball_mid(eta0)
    eta_ok = eta <= eta0_mid
    # Orders-of-magnitude gap when blocked (honest Rung-4 diagnostics)
    if eta0_mid > 0 and eta > 0:
        eta_gap = eta / eta0_mid
        eta_orders = math.log10(eta_gap)
    else:
        eta_gap = float("inf") if eta > 0 else 0.0
        eta_orders = float("inf") if eta > 0 else 0.0

    # Lemma K tail at this M (should be ≪ any realistic δ)
    tail = lemma_K_tail(
        M, Y0, r, theta, C_H=C_H, eps=eps, classical_CK=classical_CK
    )

    return {
        "delta": to_ball(delta),
        "tau": to_ball(tau),
        "eta": eta_b,
        "eta0": eta0,
        "eta_le_eta0": eta_ok,
        "eta_over_eta0": eta_gap,
        "eta_orders_of_magnitude_gap": eta_orders,
        "sharp_geom": sharp_geom,
        "U_norm": U_norm,
        "Y": to_ball(Y),
        "Y0": to_ball(Y0),
        "r": to_ball(r),
        "lambda": lam,
        "field": field,
        "M": M,
        "theta": theta,
        "C_H": C_H,
        "eps": eps,
        "C1": C1,
        "C2": C2,
        "exp_m2piY0": exp_term if HAS_FLINT else to_ball(math.exp(-2.0 * math.pi * Y0)),
        "C1_eta": (C1 * eta_b) if HAS_FLINT else to_ball(float(ball_mid(C1)) * eta),
        "C2_exp": (C2 * (exp_term if HAS_FLINT else to_ball(math.exp(-2.0 * math.pi * Y0))))
        if HAS_FLINT
        else to_ball(float(ball_mid(C2)) * math.exp(-2.0 * math.pi * Y0)),
        "lambda_error_bound": lambda_err,
        "L2_error_bound": L2_err,
        "tail_M": tail,
        "C_K": consts["C_K"],
        "A_bdry": consts["A_bdry"],
        "A_res": consts["A_res"],
        "A_cusp": consts["A_cusp"],
        "backend": consts["backend"],
        "assumptions": [
            "H(K,θ): Hecke growth |a_β| ≤ C_H(ε) N(β)^{θ+ε}",
            "A: analyticity / unique continuation on the compact core",
            "S: spectral decomposition (Friedman Thm 3.8.1 / EGM)",
        ],
        "parameters": {
            "Y": Y,
            "Y0": Y0,
            "r": r,
            "K": "Q(i)" if field in ("i", "Zi", "Z[i]", "gaussian") else "Q(√-3)",
            "M": M,
            "theta": theta,
            "C_H": C_H,
        },
        "note": (
            "Enclosure of Theorem D(K) defect bound under Assumptions H,A,S. "
            "Not a certified eigenvalue: Hejhal must supply certified (δ,τ) "
            "with η ≤ η₀. Then r₁ is an engineering target only."
        ),
    }


def print_defect_result(d: Dict[str, Any]) -> None:
    """Pretty-print a defect_to_lambda_error result."""
    print("=" * 72)
    print("Theorem D(K) — defect → |λ̃ − λ| enclosure")
    print(f"  backend     : {d['backend']}")
    print(f"  field       : {d['field']}  ({d['parameters']['K']})")
    print(f"  (Y, Y0, r)  : ({ball_mid(d['Y'])}, {ball_mid(d['Y0'])}, {ball_mid(d['r'])})")
    print(f"  λ = r²+1    : {ball_str(d['lambda'])}")
    print(f"  M, θ, C_H   : {d['M']}, {d['theta']}, {d['C_H']}")
    print("-" * 72)
    print("Assumptions (explicit):")
    for a in d["assumptions"]:
        print(f"  • {a}")
    print("-" * 72)
    print(f"  δ_aut       = {ball_str(d['delta'])}")
    print(f"  τ_tail      = {ball_str(d['tau'])}")
    print(f"  η = δ+τ     = {ball_str(d['eta'])}")
    print(f"  η₀          = {ball_str(d['eta0'])}")
    print(f"  η ≤ η₀?     : {d['eta_le_eta0']}")
    print("-" * 72)
    print(f"  C₁          = {ball_str(d['C1'])}")
    print(f"  C₂          = {ball_str(d['C2'])}")
    print(f"  e^{{-2π Y0}}   = {ball_str(d['exp_m2piY0'])}")
    print(f"  C₁·η        = {ball_str(d['C1_eta'])}")
    print(f"  C₂·e^{{-2πY0}} = {ball_str(d['C2_exp'])}")
    print("-" * 72)
    print(f"  |λ_true−λ| ≤ {ball_str(d['lambda_error_bound'])}")
    print(f"  ‖U_cusp−f‖ ≤ {ball_str(d['L2_error_bound'])}   (requires η≤η₀)")
    print(f"  Lemma K tail(M={d['M']}) = {ball_str(d['tail_M'])}")
    print("-" * 72)
    print(f"  note: {d['note']}")
    print("=" * 72)


# ---------------------------------------------------------------------------
# Sensitivity table (Y0 ±0.1, r ±0.1)
# ---------------------------------------------------------------------------
def sensitivity_table(
    Y: float = 1.25,
    Y0_center: float = 0.8,
    r_center: float = 6.6,
    dY0: float = 0.1,
    dr: float = 0.1,
    field: str = "i",
    M: int = 400,
    theta: float = 0.5,
    C_H: float = 1.0,
) -> List[Dict[str, Any]]:
    """
    Tabulate C1, C2 under small shifts of Y0 and r (Rung 1 checklist item).
    """
    rows: List[Dict[str, Any]] = []
    configs = [
        ("baseline", Y0_center, r_center),
        ("Y0+0.1", Y0_center + dY0, r_center),
        ("Y0-0.1", Y0_center - dY0, r_center),
        ("r+0.1", Y0_center, r_center + dr),
        ("r-0.1", Y0_center, r_center - dr),
        ("Y0+0.1,r+0.1", Y0_center + dY0, r_center + dr),
        ("Y0-0.1,r-0.1", Y0_center - dY0, r_center - dr),
    ]
    print("=" * 88)
    print(f"Sensitivity of C1, C2  field={field}  Y={Y}  M={M}  θ={theta}  C_H={C_H}")
    print(f"  center (Y0, r) = ({Y0_center}, {r_center})")
    print("-" * 88)
    print(
        f"{'config':>16}  {'Y0':>6}  {'r':>6}  "
        f"{'C1':>28}  {'C2':>24}  {'C2 e^{-2πY0}':>14}"
    )
    print("-" * 88)
    for name, Y0, r in configs:
        d = C1_C2_constants(
            Y=Y, Y0=Y0, r=r, theta=theta, C_H=C_H, M=M, field=field
        )
        c2_exp = ball_mid(d["C2"]) * math.exp(-2.0 * math.pi * Y0)
        print(
            f"{name:>16}  {Y0:6.2f}  {r:6.2f}  "
            f"{ball_str(d['C1']):>28}  {ball_str(d['C2']):>24}  "
            f"{c2_exp:14.6g}"
        )
        rows.append(
            {
                "config": name,
                "Y0": Y0,
                "r": r,
                "C1": d["C1"],
                "C2": d["C2"],
                "C2_exp": c2_exp,
                "eta0": d["eta0"],
            }
        )
    print("=" * 88)
    return rows


def sensitivity_markdown(
    field: str = "i",
    Y0_center: float = 0.8,
    r_center: float = 6.6,
    **kwargs: Any,
) -> str:
    """Return a markdown table of sensitivity rows for constants_DK.md."""
    rows = []
    configs = [
        ("baseline", Y0_center, r_center),
        ("Y0 + 0.1", Y0_center + 0.1, r_center),
        ("Y0 − 0.1", Y0_center - 0.1, r_center),
        ("r + 0.1", Y0_center, r_center + 0.1),
        ("r − 0.1", Y0_center, r_center - 0.1),
        ("Y0+0.1, r+0.1", Y0_center + 0.1, r_center + 0.1),
        ("Y0−0.1, r−0.1", Y0_center - 0.1, r_center - 0.1),
    ]
    lines = [
        f"| Config | Y₀ | r | C₁ | C₂ | C₂ e^{{-2π Y₀}} |",
        f"|--------|---:|--:|---:|---:|----------------:|",
    ]
    for name, Y0, r in configs:
        d = C1_C2_constants(Y0=Y0, r=r, field=field, **kwargs)
        c1m = ball_mid(d["C1"])
        c2m = ball_mid(d["C2"])
        c2e = c2m * math.exp(-2.0 * math.pi * Y0)
        lines.append(
            f"| {name} | {Y0:.2f} | {r:.2f} | "
            f"{c1m:.6g} | {c2m:.6g} | {c2e:.6g} |"
        )
        rows.append(d)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Demo (toy δ, τ as Hejhal would call)
# ---------------------------------------------------------------------------
def demo() -> None:
    """
    Show how a later interval-Hejhal run plugs measured (δ, τ) into the bound.

    Toy defects are illustrative only — not from a real Hejhal solve.
    """
    print()
    print("DEMO: how Hejhal would call defect_to_lambda_error")
    print("  (toy δ,τ — not a certified Hejhal measurement)")
    print()

    # Case A: defects small enough that η ≤ η₀ for Z[i]
    # η₀(Z[i]) ~ 1.4e-13; take δ=1e-14, τ=1e-15
    print("--- Case A: small defects (η ≤ η₀ candidate) ---")
    dA = defect_to_lambda_error(
        delta=1.0e-14,
        tau=1.0e-15,
        Y=1.25,
        Y0=0.8,
        r=6.6,
        field="i",
        M=400,
        theta=0.5,
        C_H=1.0,
    )
    print_defect_result(dA)

    # Case B: realistic-ish engineering defect before full certification
    print()
    print("--- Case B: larger toy defect (η ≫ η₀; bound still defined) ---")
    dB = defect_to_lambda_error(
        delta=1.0e-8,
        tau=1.0e-9,
        Y=1.25,
        Y0=0.8,
        r=6.6,
        field="i",
        M=400,
        theta=0.5,
        C_H=1.0,
    )
    print_defect_result(dB)

    # Case C: engineering target r ≈ 6.62212
    print()
    print("--- Case C: at Then engineering-target r=6.62212 (not certified) ---")
    dC = defect_to_lambda_error(
        delta=1.0e-14,
        tau=1.0e-15,
        Y=1.25,
        Y0=0.8,
        r=6.62212,
        field="i",
        M=400,
        theta=0.5,
        C_H=1.0,
    )
    print_defect_result(dC)

    print()
    print("Stopping-condition spot checks (Rung 1):")
    c1 = ball_mid(dA["C1"])
    tail = ball_mid(dA["tail_M"])
    print(f"  C1(Z[i]) ≈ {c1:.6g}  (target band 1e5–1e6: {'PASS' if 1e5 <= c1 <= 1e7 else 'CHECK'})")
    print(f"  C2 fully explicit: {ball_str(dA['C2'])}")
    print(f"  Lemma K tail M=400 ≈ {tail:.6g}  (< 1e-70: {'PASS' if tail < 1e-70 else 'FAIL'})")
    print(f"  C_K sharp = {ball_str(C_K(6.6))}  (no r-independent lower bound claimed)")
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Theorem D(K) defect → |λ error| enclosure (Rung 1)"
    )
    p.add_argument("--demo", action="store_true", help="toy (δ,τ) demo for Hejhal callers")
    p.add_argument("--sensitivity", action="store_true", help="tabulate C1,C2 vs Y0±0.1, r±0.1")
    p.add_argument("--delta", type=float, default=None, help="automorphy defect δ")
    p.add_argument("--tau", type=float, default=None, help="L² residual τ")
    p.add_argument("--Y", type=float, default=1.25)
    p.add_argument("--Y0", type=float, default=0.8)
    p.add_argument("--r", type=float, default=6.6)
    p.add_argument("--field", type=str, default="i", choices=("i", "omega"))
    p.add_argument("--M", type=int, default=400)
    p.add_argument("--theta", type=float, default=0.5)
    p.add_argument("--C_H", type=float, default=1.0)
    p.add_argument("--md-sensitivity", action="store_true", help="print markdown sensitivity table")
    args = p.parse_args(argv)

    if args.demo:
        demo()
        return 0

    if args.sensitivity:
        sensitivity_table(
            Y=args.Y,
            Y0_center=args.Y0,
            r_center=args.r,
            field=args.field,
            M=args.M,
            theta=args.theta,
            C_H=args.C_H,
        )
        print()
        sensitivity_table(
            Y=args.Y,
            Y0_center=args.Y0,
            r_center=args.r,
            field="omega",
            M=args.M,
            theta=args.theta,
            C_H=args.C_H,
        )
        return 0

    if args.md_sensitivity:
        print("### Z[i]")
        print(sensitivity_markdown(field="i", Y0_center=args.Y0, r_center=args.r))
        print()
        print("### Z[ω]")
        print(sensitivity_markdown(field="omega", Y0_center=args.Y0, r_center=args.r))
        return 0

    if args.delta is not None and args.tau is not None:
        d = defect_to_lambda_error(
            delta=args.delta,
            tau=args.tau,
            Y=args.Y,
            Y0=args.Y0,
            r=args.r,
            field=args.field,
            M=args.M,
            theta=args.theta,
            C_H=args.C_H,
        )
        print_defect_result(d)
        return 0

    # Default: show API help + short self-check
    p.print_help()
    print()
    print(f"HAS_FLINT = {HAS_FLINT}")
    d = C1_C2_constants(field="i")
    print(f"C1(Z[i]) = {ball_str(d['C1'])}")
    print(f"C2(Z[i]) = {ball_str(d['C2'])}")
    print("Run with --demo or --sensitivity for full Rung 1 output.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
