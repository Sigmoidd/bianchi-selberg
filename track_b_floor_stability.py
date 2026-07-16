#!/usr/bin/env python3
"""Fail-closed verifier for a closing floor run and its adaptive stability run."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from flint import arb, ctx


FLOOR_WEIGHT_FORMULA_ID = "track-b-floor-quintic-logrho-width-0.30/v1"


def _canonical_hash(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--closing-result", type=Path, required=True)
    parser.add_argument("--adaptive-result", type=Path, required=True)
    parser.add_argument("--geometry", type=Path,
                        default=Path("track_b_partition_result.json"))
    parser.add_argument("--global-partition", type=Path, default=None)
    parser.add_argument("--json-out", type=Path,
                        default=Path("track_b_floor_stability_result.json"))
    parser.add_argument("--bits", type=int, default=192)
    ns = parser.parse_args()
    ctx.prec = max(128, ns.bits)
    closing = json.loads(ns.closing_result.read_text(encoding="utf-8"))
    adaptive = json.loads(ns.adaptive_result.read_text(encoding="utf-8"))
    geometry = json.loads(ns.geometry.read_text(encoding="utf-8"))
    global_partition = (
        json.loads(ns.global_partition.read_text(encoding="utf-8"))
        if ns.global_partition is not None else None
    )
    floor_geometry = geometry.get("floor_collar_incidence")
    floor_geometry_hash = _canonical_hash(floor_geometry)
    leaves = _jsonl(Path(adaptive["leaf_audit_jsonl"]))
    events = _jsonl(Path(adaptive["parent_child_audit_jsonl"]))
    base_records = _jsonl(Path(adaptive["base_audit"]))

    leaf_ids = [q["adaptive_id"] for q in leaves]
    parent_ids = [q["parent_id"] for q in events]
    child_ids = [child for event in events for child in event["children"]]
    tree_exact = bool(
        adaptive.get("partition_exact_once_asserted", False)
        and len(leaf_ids) == len(set(leaf_ids))
        and len(parent_ids) == len(set(parent_ids))
        and all(parent not in leaf_ids for parent in parent_ids)
        and all(child in leaf_ids for child in child_ids)
        and len(leaves) == len(base_records) - len(events) + len(child_ids)
    )
    fallback_count = sum(int(q.get("bessel_fallback_count", 0)) for q in leaves)
    closing_upper = arb(closing["floor_l2_upper"]).upper()
    adaptive_upper = arb(adaptive["adaptive"]["floor_l2_upper"]).upper()
    budget = min(
        arb(closing["allowed_budget_lower"]).lower(),
        arb(adaptive["adaptive"]["allowed_budget_lower"]).lower(),
    )
    maximum = max(closing_upper, adaptive_upper)
    smaller_margin = (budget - maximum).lower()
    difference = abs(closing_upper - adaptive_upper).upper()
    allowed_difference = (arb("0.10") * smaller_margin).lower()
    refinement_requirement = bool(
        float(adaptive.get("squared_contribution_fraction", 0)) >= 0.90
        and int(adaptive.get("maximum_depth", 0)) >= 1
        and len(events) > 0
    )
    local_inputs = {
        "closing_strict_inequality": bool(closing_upper < budget),
        "adaptive_strict_inequality": bool(adaptive_upper < budget),
        "continuum_remainder_certified": bool(
            closing.get("continuum_remainder_certified", False)
        ),
        "geometry_incidence_certified": bool(
            closing.get("geometry_incidence_certified", False)
        ),
        "projected_symmetries_certified": bool(
            closing.get("projected_symmetries_certified", False)
        ),
        "witness_mass_certified": bool(closing.get("witness_mass_certified", False)),
        "closing_bessel_fallback_zero": int(
            closing.get("bessel_fallback_count", -1)
        ) == 0,
        "adaptive_bessel_fallback_zero": fallback_count == 0,
        "adaptive_tree_exact_once": tree_exact,
        "refines_at_least_90_percent_squared_remainder": refinement_requirement,
        "stability_difference_within_ten_percent_smaller_margin": bool(
            smaller_margin > 0 and difference <= allowed_difference
        ),
        "floor_width_exactly_0.30": bool(
            (arb(closing.get("floor_width", "nan")) - arb("0.30")).contains(0)
        ),
        "floor_geometry_dependency_certified": bool(
            floor_geometry and floor_geometry.get("certified", False)
        ),
        "floor_geometry_hash_matches_when_recorded": closing.get(
            "floor_geometry_incidence_hash", floor_geometry_hash
        ) == floor_geometry_hash,
    }
    stable = all(local_inputs.values())
    global_compatible = bool(
        global_partition is not None
        and global_partition.get("global_partition_certified", False)
        and global_partition.get("global_weight_bounds_certified", False)
        and global_partition.get("stability_check_passed", False)
        and not global_partition.get("provisional", True)
        and global_partition.get("floor_width") == "0.30"
        and global_partition.get("floor_weight_formula_id")
            == FLOOR_WEIGHT_FORMULA_ID
        and global_partition.get("floor_geometry_incidence_hash")
            == floor_geometry_hash
    )
    result = {
        "schema": "track-b-floor-stability/v1",
        "closing_result": str(ns.closing_result.resolve()),
        "adaptive_result": str(ns.adaptive_result.resolve()),
        "trial": closing.get("trial"),
        "trial_sha256": closing.get("trial_sha256"),
        "spectral_parameter": closing.get("spectral_parameter"),
        "closing_floor_l2_upper": str(closing_upper),
        "adaptive_floor_l2_upper": str(adaptive_upper),
        # Consumer-facing conservative local certificate fields.  The larger
        # closing bound is retained; the refined bound is only the stability
        # witness.
        "floor_l2_upper": str(maximum),
        "allowed_budget_lower": str(budget),
        "budget_lower": str(budget),
        "smaller_positive_margin_lower": str(smaller_margin),
        "difference_upper": str(difference),
        "allowed_difference_upper": str(allowed_difference),
        "local_conditions": local_inputs,
        "stability_check_passed": stable,
        "floor_residual_certified": stable,
        "continuum_remainder_certified": local_inputs[
            "continuum_remainder_certified"
        ],
        "geometry_incidence_certified": local_inputs[
            "geometry_incidence_certified"
        ],
        "projected_symmetries_certified": local_inputs[
            "projected_symmetries_certified"
        ],
        "bessel_fallback_count": fallback_count,
        "floor_width": "0.30",
        "floor_weight_formula_id": FLOOR_WEIGHT_FORMULA_ID,
        "floor_geometry_incidence_hash": floor_geometry_hash,
        "global_partition_certified": global_compatible,
        "global_weight_bounds_certified": global_compatible,
        "rung4_integrator_comparison_certified": False,
        "rung4_certified": False,
    }
    ns.json_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if stable else 2


if __name__ == "__main__":
    raise SystemExit(main())
