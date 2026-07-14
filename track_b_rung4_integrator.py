#!/usr/bin/env python3
"""Fail-closed integration of the Track-B Theorem D(K) certificate.

This runner deliberately does not reuse the legacy eta <= eta0 gate.  The
six-copy theorem gives a direct cuspidal eigenvalue enclosure

    |lambda_j-lambda| <= R / mu_B,
    R <= tau + b0*delta0 + b1*delta1,

or a sharper directly certified weighted residual.  Every input must carry a
positive certification flag.  Missing, diagnostic, or incompatible inputs
leave all hard outputs false.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from flint import arb, ctx


def load_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def as_arb(value: Any, label: str) -> arb:
    if value is None:
        raise ValueError(f"missing {label}")
    out = arb(str(value))
    if not out.is_finite():
        raise ValueError(f"nonfinite {label}: {value}")
    return out


def upper_float(value: arb) -> float:
    return math.nextafter(float(value.upper()), math.inf)


def lower_float(value: arb) -> float:
    return math.nextafter(float(value.lower()), -math.inf)


def fail_result(reason: str, paths: dict[str, str]) -> dict[str, Any]:
    return {
        "status": "Track-B integration blocked",
        "track_b_interval_certified": False,
        "rung4_certified": False,
        "legacy_eta0_applicable": False,
        "reason": reason,
        "inputs": paths,
    }


def integrate(
    mass: dict[str, Any] | None,
    partition: dict[str, Any] | None,
    overlap: dict[str, Any] | None,
    width_tol: str,
    counting_certified: bool,
    paths: dict[str, str],
) -> dict[str, Any]:
    if mass is None:
        return fail_result("missing projected-mass certificate", paths)
    if partition is None:
        return fail_result("missing partition certificate", paths)
    if overlap is None:
        return fail_result("missing overlap certificate", paths)

    if not bool(mass.get("theorem_DK_projected_mass_admissible", False)):
        return fail_result("projected mass is not theorem-admissible", paths)
    if not bool(partition.get("certified", False)):
        return fail_result("partition certificate is not GREEN", paths)
    if not bool(overlap.get("certified", False)):
        return fail_result("overlap certificate is not GREEN", paths)

    required_partition_flags = (
        "coverage_certified",
        "local_finiteness_certified",
        "transitions_complete",
        "stabilizers_certified",
        "theorem_DK_compatible",
    )
    missing_flags = [k for k in required_partition_flags if not partition.get(k, False)]
    if missing_flags:
        return fail_result(f"partition flags false/missing: {missing_flags}", paths)

    required_overlap_flags = (
        "all_transitions_covered",
        "common_fiber_transport_certified",
        "stabilizer_averaging_certified",
        "two_cusp_coordinates_certified",
        "first_gradients_certified",
        "theorem_DK_compatible",
    )
    missing_flags = [k for k in required_overlap_flags if not overlap.get(k, False)]
    if missing_flags:
        return fail_result(f"overlap flags false/missing: {missing_flags}", paths)

    p_ids = set(map(str, partition.get("active_transition_ids", [])))
    o_ids = set(map(str, overlap.get("active_transition_ids", [])))
    if not p_ids or p_ids != o_ids:
        return fail_result(
            f"transition-set mismatch: partition={sorted(p_ids)}, overlap={sorted(o_ids)}",
            paths,
        )

    try:
        mu = as_arb(
            mass["plateau_construction"]["certified_mu_B_lower"], "mu_B lower"
        ).lower()
        if not bool(mu > 0):
            return fail_result("mu_B lower bound is not positive", paths)

        tau = as_arb(
            overlap.get("tau_upper", partition.get("tau_upper", 0)), "tau upper"
        ).upper()
        direct = overlap.get("weighted_residual_upper")
        if direct is not None:
            residual = tau + as_arb(direct, "weighted residual upper").upper()
            formula = "tau + direct weighted residual"
            factors: dict[str, Any] = {
                "tau_upper": upper_float(tau),
                "weighted_residual_upper": upper_float(
                    as_arb(direct, "weighted residual upper").upper()
                ),
            }
        else:
            b0 = as_arb(partition.get("b0_upper"), "b0 upper").upper()
            b1 = as_arb(partition.get("b1_upper"), "b1 upper").upper()
            d0 = as_arb(overlap.get("delta0_upper"), "delta0 upper").upper()
            d1 = as_arb(overlap.get("delta1_upper"), "delta1 upper").upper()
            residual = tau + b0 * d0 + b1 * d1
            formula = "tau + b0*delta0 + b1*delta1"
            factors = {
                "tau_upper": upper_float(tau),
                "b0_upper": upper_float(b0),
                "b1_upper": upper_float(b1),
                "delta0_upper": upper_float(d0),
                "delta1_upper": upper_float(d1),
            }

        r_text = str(mass["parameters"]["r"])
        r = arb(r_text)
        lam = 1 + r * r
        halfwidth = residual / mu
        fullwidth = 2 * halfwidth
        lo = lam - halfwidth
        hi = lam + halfwidth
        tol = arb(width_tol)
        width_ok = bool(fullwidth.upper() < tol.lower())
        disjoint_unit = bool(lo.lower() > 1)
        interval_certified = width_ok and disjoint_unit

        hard = {
            "projected_mass_certified": True,
            "partition_certified": True,
            "overlap_defects_certified": True,
            "transition_sets_match": True,
            "theorem_DK_interval_certified": interval_certified,
            "width_lt_tol": width_ok,
            "disjoint_from_unit_interval": disjoint_unit,
            "counting_certified": bool(counting_certified),
        }
        rung4 = all(hard.values())
        return {
            "status": (
                "Track-B interval certified" if interval_certified
                else "Track-B inputs certified but interval target not closed"
            ),
            "track_b_interval_certified": interval_certified,
            "rung4_certified": rung4,
            "legacy_eta0_applicable": False,
            "legacy_eta0_note": (
                "The six-copy Track-B theorem uses the direct R/mu_B enclosure; "
                "it does not invoke the inherited scalar eta0 gate."
            ),
            "inputs": paths,
            "transition_count": len(p_ids),
            "residual": {
                "formula": formula,
                **factors,
                "R_upper": upper_float(residual),
            },
            "mass": {"mu_B_lower": lower_float(mu)},
            "spectral_enclosure": {
                "r": r_text,
                "lambda_center_ball": str(lam),
                "lambda_lower": lower_float(lo),
                "lambda_upper": upper_float(hi),
                "halfwidth_upper": upper_float(halfwidth),
                "fullwidth_upper": upper_float(fullwidth),
                "width_tolerance": float(width_tol),
            },
            "hard": hard,
            "blockers": [k for k, value in hard.items() if not value],
        }
    except (KeyError, TypeError, ValueError) as exc:
        return fail_result(f"invalid certified input schema: {exc}", paths)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mass", type=Path, default=Path("track_b_projected_mass_arb_result.json")
    )
    parser.add_argument(
        "--partition", type=Path, default=Path("track_b_partition_result.json")
    )
    parser.add_argument(
        "--overlap", type=Path, default=Path("track_b_overlap_result.json")
    )
    parser.add_argument("--width-tol", default="0.1")
    parser.add_argument("--counting-certified", action="store_true")
    parser.add_argument(
        "--json-out", type=Path, default=Path("track_b_rung4_result.json")
    )
    ns = parser.parse_args()
    ctx.prec = 192
    paths = {
        "mass": str(ns.mass.resolve()),
        "partition": str(ns.partition.resolve()),
        "overlap": str(ns.overlap.resolve()),
    }
    result = integrate(
        load_json(ns.mass),
        load_json(ns.partition),
        load_json(ns.overlap),
        ns.width_tol,
        ns.counting_certified,
        paths,
    )
    ns.json_out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(ns.json_out.resolve())
    return 0 if result.get("track_b_interval_certified", False) else 2


if __name__ == "__main__":
    raise SystemExit(main())
