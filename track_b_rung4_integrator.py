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
import hashlib
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


def canonical_hash(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def validate_partition_artifact(
    partition: dict[str, Any], floor: dict[str, Any] | None
) -> str | None:
    """Return a fail-closed incompatibility reason, or ``None``.

    This is deliberately stronger than checking the top-level GREEN label:
    the serialized definition and configuration are rehashed, and every
    audit leaf must carry the same definition/geometry/formula/configuration
    identifiers as the result that consumes it.
    """
    required = (
        "global_partition_certified", "global_weight_bounds_certified",
        "coverage_certified", "local_finiteness_certified",
        "denominator_positive_certified", "partition_sum_certified",
        "stabilizer_averages_certified", "weight_gradients_certified",
        "weight_laplacians_certified", "floor_weight_consistency_certified",
        "geometry_incidence_certified", "stability_check_passed",
        "partition_constants_certified",
    )
    missing = [name for name in required if not partition.get(name, False)]
    if missing:
        return f"partition flags false/missing: {missing}"
    if partition.get("provisional", True):
        return "partition result is provisional or lacks a stable closing run"
    if int(partition.get("unresolved_fallback_count", -1)) != 0:
        return "partition derivative ledger used an unresolved fallback"
    if not partition.get("per_active_weight_enclosures_recorded", False):
        return "partition lacks per-active-weight interval enclosures"

    definition = partition.get("partition_definition")
    definition_hash = partition.get("partition_definition_hash")
    if not isinstance(definition, dict) or canonical_hash(definition) != definition_hash:
        return "partition definition hash mismatch"
    dims = partition.get("subdivision")
    config = {
        "subdivision": dims,
        "Taylor_degree": partition.get("partition_degree"),
        "arb_bits": partition.get("arb_bits"),
        "s_cap": partition.get("coverage_proof", {}).get("s_cap"),
    }
    if canonical_hash(config) != partition.get("configuration_hash"):
        return "partition configuration hash mismatch"
    identifiers = {
        "partition_definition_hash": definition_hash,
        "geometry_certificate_dependency": partition.get("geometry_incidence_hash"),
        "weight_formulas_hash": partition.get("weight_formulas_hash"),
        "configuration_hash": partition.get("configuration_hash"),
    }
    if any(not isinstance(value, str) or len(value) != 64
           for value in identifiers.values()):
        return "partition contains a missing or malformed deterministic identifier"

    audit_text = partition.get("audit_jsonl")
    if not audit_text:
        return "partition result has no audit ledger path"
    audit_path = Path(audit_text)
    if not audit_path.is_file():
        return "partition audit ledger is missing"
    digest = hashlib.sha256(audit_path.read_bytes()).hexdigest()
    if digest != partition.get("audit_sha256"):
        return "partition audit ledger hash mismatch"
    try:
        count = 0
        ids: set[str] = set()
        with audit_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                row = json.loads(line)
                count += 1
                ids.add(str(row.get("cell_id")))
                if any(row.get(key) != value for key, value in identifiers.items()):
                    return "partition audit/result identifier mismatch"
                if row.get("derivative_bound_method") != (
                    "exact factorwise product-rule majorant"
                ):
                    return "partition contains midpoint-only or unknown derivative bounds"
                if not row.get("per_active_weight_enclosures_recorded", False):
                    return "partition audit omits per-active-weight enclosures"
                for weight in row.get("per_active_weight_enclosures", []):
                    for field in (
                        "grad_chi_hyperbolic_certified_upper",
                        "abs_Delta_chi_certified_upper",
                    ):
                        if not as_arb(weight.get(field), field).is_finite():
                            return "partition audit contains an unbounded derivative"
        expected = int(dims[0]) * int(dims[1]) * int(dims[2])
        if count != expected or len(ids) != expected:
            return "partition audit is incomplete or contains duplicate cells"
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
        return f"invalid partition audit ledger: {exc}"

    if floor is not None:
        try:
            pwidth = as_arb(partition.get("floor_width"), "partition floor width")
            fwidth = as_arb(floor.get("floor_width"), "floor certificate width")
        except ValueError as exc:
            return str(exc)
        if not (pwidth - fwidth).contains(0):
            return "floor width mismatch between partition and local certificate"
        if partition.get("floor_weight_formula_id") != floor.get(
            "floor_weight_formula_id"
        ):
            return "floor weight formula mismatch"
        if partition.get("floor_geometry_incidence_hash") != floor.get(
            "floor_geometry_incidence_hash"
        ):
            return "floor/global geometry incidence hash mismatch"
    return None


def integrate(
    mass: dict[str, Any] | None,
    partition: dict[str, Any] | None,
    overlap: dict[str, Any] | None,
    width_tol: str,
    counting_certified: bool,
    paths: dict[str, str],
    floor: dict[str, Any] | None = None,
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
    partition_problem = validate_partition_artifact(partition, floor)
    if partition_problem is not None:
        return fail_result(partition_problem, paths)
    if not bool(overlap.get("certified", False)):
        return fail_result("overlap certificate is not GREEN", paths)

    required_partition_flags = (
        "coverage_certified",
        "local_finiteness_certified",
        "transitions_complete",
        "stabilizers_certified",
        "theorem_DK_compatible",
        "partition_constants_certified",
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
        floor_gate = None
        if floor is not None:
            required_floor_flags = (
                "floor_residual_certified",
                "continuum_remainder_certified",
                "stability_check_passed",
                "geometry_incidence_certified",
                "projected_symmetries_certified",
            )
            missing_floor = [k for k in required_floor_flags if not floor.get(k, False)]
            if missing_floor:
                return fail_result(
                    f"floor certificate flags false/missing: {missing_floor}", paths
                )
            if int(floor.get("bessel_fallback_count", -1)) != 0:
                return fail_result("floor certificate used a Bessel fallback", paths)
            floor_upper = as_arb(floor.get("floor_l2_upper"), "floor L2 upper").upper()
            budget_lower = as_arb(
                floor.get("allowed_budget_lower"), "floor allowed budget lower"
            ).lower()
            if not bool(floor_upper < budget_lower):
                return fail_result(
                    "floor endpoint comparison does not prove upper(residual) < lower(budget)",
                    paths,
                )
            direct = floor_upper
            floor_gate = {
                "floor_l2_upper_ball": str(floor_upper),
                "allowed_budget_lower_ball": str(budget_lower),
                "endpoint_comparison_certified": True,
            }
        else:
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
        if floor is not None:
            hard["floor_residual_certified"] = True
            hard["floor_stability_certified"] = True
            hard["floor_bessel_fallback_zero"] = True
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
                "arb_upper_ball": str(residual.upper()),
            },
            "floor_comparison": floor_gate,
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
    parser.add_argument("--floor", type=Path, default=None)
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
        "floor": None if ns.floor is None else str(ns.floor.resolve()),
    }
    result = integrate(
        load_json(ns.mass),
        load_json(ns.partition),
        load_json(ns.overlap),
        ns.width_tol,
        ns.counting_certified,
        paths,
        load_json(ns.floor) if ns.floor is not None else None,
    )
    ns.json_out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(ns.json_out.resolve())
    return 0 if result.get("track_b_interval_certified", False) else 2


if __name__ == "__main__":
    raise SystemExit(main())
