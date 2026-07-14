#!/usr/bin/env python3
r"""Arb lower bound for the Track-B rotation-odd oldspace witness.

This certifies the projected-mass input of Theorem D(K) for one fixed finite
Whittaker trial.  It constructs the cusp-local witness

    W_B = P_{Q,-} P_{-z,+} P_old F,

where P_old averages the six induced components, P_{-z,+} enforces the exact
Picard relation z ~ -z, and P_{Q,-}=(I-Q)/2 with Qz=iz.  The final global
quasimode U will satisfy

    ||P_{Q,-} P_old U|| >= ||W_B||,

because the explicit C2 cusp cutoff is identically one on the witness slab.
The full theorem remains uncertified until the smooth-gluing residual is
bounded.

All coefficient algebra, Bessel values, Gram eigenvalues, and y integration
are evaluated with python-flint Arb/Acb.  Two norm lower bounds are computed:

1. exact full-torus Parseval divided by two, justified by W_B(-z)=W_B(z);
2. the direct half-torus Fourier Gram with Rump-certified eigenvalues.

A separate pointwise check compares the collapsed Fourier formula with a
direct six-component evaluation before the norm integration.
"""
from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from flint import acb, acb_mat, arb, ctx

from continuum_defect_arb import TrialEvaluator, lower, parse_trial, upper


def projected_coefficients(ev: TrialEvaluator) -> dict[str, Any]:
    """Collapse six copies, enforce -z invariance, then project Q-odd."""
    inf_index = {(a, b): k for k, (a, b, _nn) in enumerate(ev.modes_inf)}
    zero_index = {
        (a, b): ev.ni + k for k, (a, b, _nn) in enumerate(ev.modes_0)
    }
    if len(inf_index) != len(ev.modes_inf):
        raise AssertionError("duplicate infinity mode")
    if len(zero_index) != len(ev.modes_0):
        raise AssertionError("duplicate zero-cusp mode")

    old: list[acb] = []
    lift_records = []
    for k, (u, v, nn) in enumerate(ev.modes_inf):
        # Invert ((2a-b)/5,(a+2b)/5)=(u,v).
        a0, b0 = 2 * u + v, -u + 2 * v
        if a0 * a0 + b0 * b0 != 5 * nn:
            raise AssertionError("incorrect zero-cusp lift norm")
        kz = zero_index.get((a0, b0))
        if kz is None:
            raise AssertionError(f"missing lifted zero-cusp mode {(a0, b0)}")
        # Sum over translations t=0,...,4.  The lifted frequency u is an
        # integer, so every phase is one and the sum is exactly five.  All
        # other zero-cusp modes cancel by the fifth-root geometric sum.
        old.append((ev.coeff[k] + 5 * ev.coeff[kz]) / 6)
        lift_records.append(
            {
                "infinity_mode": [u, v],
                "zero_mode": [a0, b0],
                "norm_relation": f"{a0*a0+b0*b0}=5*{nn}",
            }
        )

    even: list[acb] = []
    for k, (u, v, _nn) in enumerate(ev.modes_inf):
        even.append((old[k] + old[inf_index[(-u, -v)]]) / 2)

    odd: list[acb] = []
    for k, (u, v, _nn) in enumerate(ev.modes_inf):
        # The coefficient of exp(2*pi*i*(u*x1+v*x2)) in f(Qz) is
        # the coefficient of mode (-v,u) in f.
        odd.append((even[k] - even[inf_index[(-v, u)]]) / 2)

    symmetry_defects = {"minus_z_even": arb(0), "quarter_turn_odd": arb(0)}
    for k, (u, v, _nn) in enumerate(ev.modes_inf):
        symmetry_defects["minus_z_even"] = max(
            symmetry_defects["minus_z_even"],
            abs(odd[k] - odd[inf_index[(-u, -v)]]).upper(),
        )
        symmetry_defects["quarter_turn_odd"] = max(
            symmetry_defects["quarter_turn_odd"],
            abs(odd[k] + odd[inf_index[(-v, u)]]).upper(),
        )

    # Count the modes selected by the five-translation average independently.
    selected = []
    for a, b, _nn in ev.modes_0:
        if (2 * a - b) % 5 == 0 and (a + 2 * b) % 5 == 0:
            selected.append(((2 * a - b) // 5, (a + 2 * b) // 5))
    if set(selected) != set(inf_index):
        raise AssertionError("translation-selected zero-cusp lattice mismatch")

    return {
        "old": old,
        "even": even,
        "odd": odd,
        "lift_records": lift_records,
        "translation_selected_modes": len(selected),
        "symmetry_defect_balls": {k: str(v) for k, v in symmetry_defects.items()},
        "symmetry_defect_uppers": {k: upper(v) for k, v in symmetry_defects.items()},
    }


def half_gram_lowers(ev: TrialEvaluator) -> tuple[dict[int, arb], dict[str, Any]]:
    groups: dict[int, list[int]] = {}
    for k, (a, _b, _nn) in enumerate(ev.modes_inf):
        groups.setdefault(a, []).append(k)
    lowers: dict[int, arb] = {}
    records: dict[str, Any] = {}
    for a, ids in groups.items():
        gram = acb_mat(len(ids), len(ids))
        for ii, k in enumerate(ids):
            b = ev.modes_inf[k][1]
            for jj, ell in enumerate(ids):
                d = b - ev.modes_inf[ell][1]
                if d == 0:
                    gram[ii, jj] = acb(arb(1) / 2)
                elif d % 2 == 0:
                    gram[ii, jj] = acb(0)
                else:
                    gram[ii, jj] = acb(0, arb(1) / (ev.pi * d))
        eigs = gram.eig(algorithm="rump")
        if any(not eig.imag.contains(0) for eig in eigs):
            raise ArithmeticError(f"nonreal half-Gram eigenvalue in block {a}")
        lam = eigs[0].real.lower()
        for eig in eigs[1:]:
            candidate = eig.real.lower()
            if bool(candidate < lam):
                lam = candidate
        if not bool(lam > 0):
            raise ArithmeticError(f"half-Gram block {a} not certified positive")
        lowers[a] = lam
        records[str(a)] = {"size": len(ids), "lambda_min_lower": lower(lam)}
    return lowers, {"groups": groups, "records": records}


def integrate_projected_norm(
    ev: TrialEvaluator, odd: list[acb], ny: int, y_min: str, y_max: str
) -> dict[str, Any]:
    """Return two certified vector-valued norm lower bounds."""
    ya, yb = arb(y_min), arb(y_max)
    if not bool(ya > 1) or not bool(yb > ya):
        raise ValueError("require 1 < y_min < y_max in the embedded cusp")
    edges = [ya + (yb - ya) * j / ny for j in range(ny + 1)]
    gram_lowers, gram_data = half_gram_lowers(ev)
    groups: dict[int, list[int]] = gram_data.pop("groups")

    scalar_parseval_half = arb(0)
    scalar_direct_half = arb(0)
    positive_parseval = 0
    positive_gram = 0
    for j in range(ny):
        ybox = edges[j].union(edges[j + 1])
        energies: list[arb] = []
        for k, mode in enumerate(ev.modes_inf):
            _u, _v, mag = ev._frequency(0, mode)
            radial = odd[k] * ev._kir(0, k, mag, ybox)
            mag_lo = abs(radial).lower()
            energies.append(mag_lo * mag_lo if bool(mag_lo > 0) else arb(0))

        # W_B is exactly even under z -> -z.  Its integral over the chosen
        # half torus is therefore one half of full-torus Parseval.
        parseval_planar = sum(energies, arb(0)) / 2
        parseval_lo = (parseval_planar / ybox).lower()
        if bool(parseval_lo > 0):
            scalar_parseval_half += (edges[j + 1] - edges[j]) * parseval_lo
            positive_parseval += 1

        gram_planar = arb(0)
        for a, ids in groups.items():
            gram_planar += gram_lowers[a] * sum(
                (energies[k] for k in ids), arb(0)
            )
        gram_lo = (gram_planar / ybox).lower()
        if bool(gram_lo > 0):
            scalar_direct_half += (edges[j + 1] - edges[j]) * gram_lo
            positive_gram += 1

    # Every component of P_old F equals the scalar old average, so the
    # six-component bundle norm is six times the scalar norm squared.
    vector_parseval = 6 * scalar_parseval_half
    vector_gram = 6 * scalar_direct_half
    if not bool(vector_parseval > 0) or not bool(vector_gram > 0):
        raise ArithmeticError("failed to certify positive Track-B witness mass")
    mu_parseval = vector_parseval.sqrt().lower()
    mu_gram = vector_gram.sqrt().lower()
    return {
        "domain": (
            "Picard half torus [-1/2,1/2]x[0,1/2], "
            f"y in [{y_min},{y_max}]"
        ),
        "hyperbolic_measure_reduction": "|y K|^2 dx1 dx2 dy/y^3 = |K|^2 dx1 dx2 dy/y",
        "y_segments": ny,
        "positive_segments": {
            "parseval": positive_parseval,
            "direct_half_gram": positive_gram,
        },
        "parseval_via_evenness": {
            "method": "full-torus Parseval divided by 2 using W_B(-z)=W_B(z)",
            "vector_norm2_lower": lower(vector_parseval),
            "vector_norm_lower": lower(mu_parseval),
            "vector_norm2_ball": str(vector_parseval),
        },
        "direct_half_gram": {
            "method": "half-period Fourier Gram; eigenvalues certified by Rump",
            "blocks": gram_data["records"],
            "vector_norm2_lower": lower(vector_gram),
            "vector_norm_lower": lower(mu_gram),
            "vector_norm2_ball": str(vector_gram),
        },
        "certified_witness_norm_lower": max(lower(mu_parseval), lower(mu_gram)),
    }


def fourier_value(
    ev: TrialEvaluator, coeff: list[acb], x1: arb, x2: arb, y: arb
) -> acb:
    total = acb(0)
    for k, mode in enumerate(ev.modes_inf):
        u, v, mag = ev._frequency(0, mode)
        phase = acb(0, 2 * ev.pi * (u * x1 + v * x2)).exp()
        total += coeff[k] * y * ev._kir(0, k, mag, y) * phase
    return total


def direct_old(ev: TrialEvaluator, x1: arb, x2: arb, y: arb) -> acb:
    return sum((ev.component(c, x1, x2, y) for c in range(6)), acb(0)) / 6


def pointwise_crosscheck(ev: TrialEvaluator, odd: list[acb]) -> dict[str, Any]:
    points = [
        (arb("0.125"), arb("0.0625"), arb("1.125")),
        (arb("-0.2"), arb("0.3"), arb("1.2")),
    ]
    records = []
    max_difference = arb(0)
    all_overlap = True
    for x1, x2, y in points:
        def even_old(a: arb, b: arb) -> acb:
            return (direct_old(ev, a, b, y) + direct_old(ev, -a, -b, y)) / 2

        direct = (even_old(x1, x2) - even_old(-x2, x1)) / 2
        collapsed = fourier_value(ev, odd, x1, x2, y)
        diff = direct - collapsed
        overlaps = diff.real.contains(0) and diff.imag.contains(0)
        all_overlap = all_overlap and overlaps
        max_difference = max(max_difference, abs(diff).upper())
        records.append(
            {
                "point": [str(x1), str(x2), str(y)],
                "direct_ball": str(direct),
                "collapsed_ball": str(collapsed),
                "difference_ball": str(diff),
                "difference_contains_zero": overlaps,
            }
        )
    if not all_overlap:
        raise ArithmeticError("direct/collapsed Track-B pointwise check disagrees")
    return {
        "method_a": "direct average of six TrialEvaluator components at z,-z,Qz,-Qz",
        "method_b": "collapsed integral-lattice Fourier coefficients",
        "all_difference_balls_contain_zero": all_overlap,
        "max_difference_abs_upper": upper(max_difference),
        "records": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--trial", type=Path, default=Path("six_copy_hejhal_balanced_coeffs.json")
    )
    parser.add_argument("--bits", type=int, default=160)
    parser.add_argument("--y-segments", type=int, default=256)
    parser.add_argument("--y-min", default="1.20")
    parser.add_argument("--y-max", default="1.45")
    parser.add_argument("--cutoff-support-y", default="1.01")
    parser.add_argument("--plateau-y", default="1.20")
    parser.add_argument(
        "--json-out", type=Path, default=Path("track_b_projected_mass_arb_result.json")
    )
    parser.add_argument("--quiet", action="store_true")
    ns = parser.parse_args()
    ctx.prec = ns.bits
    cutoff_support = Decimal(ns.cutoff_support_y)
    plateau_y = Decimal(ns.plateau_y)
    witness_y_min = Decimal(ns.y_min)
    plateau_compatible = bool(
        cutoff_support > Decimal(1)
        and plateau_y > cutoff_support
        and witness_y_min >= plateau_y
    )
    if not plateau_compatible:
        raise ValueError(
            "require 1 < cutoff_support_y < plateau_y <= witness y_min"
        )

    data, row = parse_trial(ns.trial)
    ev = TrialEvaluator(data["parameters"]["M"], str(row["r"]), row["coefficients"])
    algebra = projected_coefficients(ev)
    norm = integrate_projected_norm(
        ev, algebra["odd"], ns.y_segments, ns.y_min, ns.y_max
    )
    crosscheck = pointwise_crosscheck(ev, algebra["odd"])

    digest = hashlib.sha256(ns.trial.read_bytes()).hexdigest()
    result = {
        "status": "Arb Track-B projected-mass witness for fixed finite trial",
        "witness_certified": True,
        "theorem_DK_projected_mass_admissible": True,
        "theorem_DK_remaining_blocker": (
            "certify the Track-B smooth-gluing residual: transition value/gradient "
            "defects and partition constants; projected mass is closed"
        ),
        "trial": str(ns.trial.resolve()),
        "trial_sha256": digest,
        "parameters": {
            "M_infinity": ev.M,
            "M_zero": 5 * ev.M,
            "r": str(row["r"]),
            "precision_bits": ns.bits,
            "y_segments": ns.y_segments,
            "y_min": ns.y_min,
            "y_max": ns.y_max,
        },
        "plateau_construction": {
            "embedded_horoball": "y > 1 modulo the Picard parabolic subgroup",
            "cutoff_support_y": ns.cutoff_support_y,
            "cutoff_plateau_y": ns.plateau_y,
            "transition_variable": (
                f"t=(y-{ns.cutoff_support_y})/"
                f"({ns.plateau_y}-{ns.cutoff_support_y})"
            ),
            "C2_smoothstep": "10*t^3-15*t^4+6*t^5",
            "cutoff": "0 below support, smoothstep on transition, 1 on plateau",
            "other_partition_weights": (
                "(1-chi_B)*phi_j for a normalized subordinate family sum phi_j=1"
            ),
            "witness_slab_inside_plateau": plateau_compatible,
            "projected_automorphization_error_on_witness_slab": 0,
            "certified_mu_B_lower": norm["certified_witness_norm_lower"],
        },
        "projector": {
            "formula": "P_(Q,-) P_(-z,+) P_old",
            "old_copy_average": "(F_0+...+F_5)/6",
            "minus_z_average": "(f(z)+f(-z))/2",
            "rotation_odd": "(f(z)-f(i z))/2",
            "bundle_norm_factor": 6,
            "simplicity_assumption": False,
        },
        "finite_algebra": {
            "translation_selected_modes": algebra["translation_selected_modes"],
            "infinity_modes": len(ev.modes_inf),
            "lift_bijection_proved": algebra["translation_selected_modes"]
            == len(ev.modes_inf),
            "symmetry_defect_balls": algebra["symmetry_defect_balls"],
            "symmetry_defect_uppers": algebra["symmetry_defect_uppers"],
            "lift_records": algebra["lift_records"],
        },
        "pointwise_independent_crosscheck": crosscheck,
        "norm_certificate": norm,
        "bessel_enclosures": {
            "direct_arb": ev.direct_k_count,
            "mean_value_fallback": ev.mean_value_k_count,
        },
        "hard_map_changed": False,
    }
    ns.json_out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    if not ns.quiet:
        print(json.dumps(result, indent=2))
    else:
        print(
            json.dumps(
                {
                    "json_out": str(ns.json_out.resolve()),
                    "bits": ns.bits,
                    "y_segments": ns.y_segments,
                    "witness_norm_lower": norm["certified_witness_norm_lower"],
                    "crosscheck": crosscheck["all_difference_balls_contain_zero"],
                }
            )
        )
    print(ns.json_out.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
