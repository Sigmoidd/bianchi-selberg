#!/usr/bin/env python3
"""Independent fail-closed verifier for a global Track-B Hejhal defect."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from flint import arb, ctx

from track_b_two_cusp_data import canonical_hash


ROOT = Path(__file__).resolve().parent
EXPECTED_CHANNELS = {
    "floor",
    "cusp_infinity_continuum",
    "cusp_zero_continuum",
    "nonfloor_face_value",
    "nonfloor_face_gradient",
    "unknown_infinite_fourier_tail",
    "reprojection",
}
ALLOWED_RULES = {"triangle", "disjoint_l2", "exact_orthogonality", "certified_gram"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reconstruct_total(result: dict[str, Any]) -> tuple[arb | None, list[str]]:
    errors: list[str] = []
    channels = result.get("channel_ledger", [])
    ids = [q.get("id") for q in channels]
    if len(ids) != len(set(ids)):
        errors.append("duplicate channel contribution")
    if set(ids) != EXPECTED_CHANNELS:
        errors.append("required channel set is incomplete or contains an unknown channel")
    by_id = {q.get("id"): q for q in channels}
    for channel in EXPECTED_CHANNELS:
        record = by_id.get(channel)
        if record is None or not record.get("certified", False) or record.get("upper") is None:
            errors.append(f"missing or uncertified channel: {channel}")
    if errors:
        return None, errors

    plan = result.get("aggregation_plan", [])
    used: list[str] = []
    group_bounds: list[arb] = []
    for group in plan:
        rule = group.get("rule")
        members = group.get("channels", [])
        if rule not in ALLOWED_RULES:
            errors.append(f"unjustified aggregation rule: {rule}")
            continue
        if rule in {"disjoint_l2", "exact_orthogonality", "certified_gram"} and not group.get(
            "justification_certified", False
        ):
            errors.append(f"quadratic/Gram aggregation lacks proof: {rule}")
            continue
        if any(member not in EXPECTED_CHANNELS for member in members):
            errors.append("aggregation plan names an unknown channel")
            continue
        used.extend(members)
        values = [arb(str(by_id[member]["upper"])).upper() for member in members]
        if rule == "triangle":
            group_bounds.append(sum(values, arb(0)).upper())
        elif rule in {"disjoint_l2", "exact_orthogonality"}:
            group_bounds.append(sum((q * q for q in values), arb(0)).sqrt().upper())
        else:
            gram_upper = group.get("certified_gram_upper")
            if gram_upper is None:
                errors.append("certified_gram rule lacks certified_gram_upper")
            else:
                group_bounds.append(arb(str(gram_upper)).upper())
    if sorted(used) != sorted(EXPECTED_CHANNELS):
        errors.append("aggregation plan does not consume every channel exactly once")
    if errors:
        return None, errors
    return sum(group_bounds, arb(0)).upper(), errors


def structural_checks(result: dict[str, Any]) -> dict[str, bool]:
    total, errors = reconstruct_total(result)
    serialized_total = result.get("global_defect_upper")
    total_match = False
    if total is not None and serialized_total is not None:
        total_match = (arb(str(serialized_total)) - total).contains(0)
    threshold = result.get("allowed_defect_lower")
    strict = bool(
        total is not None and threshold is not None and total < arb(str(threshold)).lower()
    )
    return {
        "schema_exact": result.get("schema") == "track-b-global-hejhal-defect/v1",
        "expected_channels_exact_once": not any("channel" in q and ("duplicate" in q or "incomplete" in q) for q in errors),
        "all_channels_certified": total is not None,
        "gradient_pullback_certified": bool(result.get("gradient_pullback_certified", False)),
        "fourier_tail_certified": bool(result.get("fourier_tails_certified", False)),
        "tail_majorant_retained_mode_count_zero": int(result.get("tail_majorant_retained_mode_count", -1)) == 0,
        "reprojection_hashes_match": bool(result.get("reprojection_hashes_match", False)),
        "floor_certificate_compatible": bool(result.get("floor_certificate_compatible", False)),
        "collocation_not_substituted": result.get("physical_collocation_is_theorem_channel") is False,
        "aggregation_reconstructed": total is not None,
        "serialized_total_matches": total_match,
        "global_threshold_certified": bool(result.get("global_threshold_certified", False)),
        "threshold_not_local_substitution": bool(result.get("threshold_trial_compatibility_certified", False)),
        "cutoff_stability_passed": bool(result.get("cutoff_stability_passed", False)),
        "strict_endpoint_inequality": strict,
        "resolved_retained_bessel_fallback_zero": int(result.get("resolved_bessel_fallback_count", -1)) == 0,
    }


def verify(ns: argparse.Namespace) -> dict[str, Any]:
    ctx.prec = max(128, ns.bits)
    result = json.loads(ns.result.read_text(encoding="utf-8"))
    theorem = json.loads(ns.theorem_defect_definition.read_text(encoding="utf-8"))
    partition = json.loads(ns.global_partition.read_text(encoding="utf-8"))
    trial = json.loads(ns.two_cusp_result.read_text(encoding="utf-8"))
    floor = json.loads(ns.floor_certificate.read_text(encoding="utf-8"))
    checks = structural_checks(result)
    checks.update({
        "theorem_hash_matches": result.get("theorem_defect_definition_hash")
            == theorem.get("theorem_defect_definition_hash"),
        "partition_hash_matches": result.get("partition_hash")
            == partition.get("partition_definition_hash"),
        "weight_hash_matches": result.get("weight_hash")
            == partition.get("weight_formulas_hash"),
        "floor_artifact_hash_matches": result.get("floor_certificate_hash") == sha256(ns.floor_certificate),
        "two_cusp_assembly_hash_matches": result.get("two_cusp_assembly_hash") == trial.get("assembly_hash"),
        "coefficient_vector_hash_matches": result.get("coefficient_vector_hash") == trial.get("coefficient_vector_hash"),
        "coefficient_to_field_map_hash_matches": result.get("coefficient_to_field_map_hash")
            == trial.get("coefficient_to_field_map_hash"),
        "channel_ledger_hash_matches": result.get("aggregation_ledger_hash")
            == canonical_hash(result.get("channel_ledger", [])),
        "serialized_claim_is_fail_closed": not result.get("global_hejhal_defect_certified", False),
        "rung4_claim_is_fail_closed": not result.get("rung4_certified", False),
        "hejhal_existence_false": result.get("hejhal_existence_certified") is False,
        "dual_certification_false": result.get("dual_certification") is False,
    })
    certified = all(checks.values())
    return {
        "schema": "track-b-global-hejhal-independent-verification/v1",
        "result_path": str(ns.result.resolve()),
        "result_sha256": sha256(ns.result),
        "checks": checks,
        "failed_checks": [key for key, value in checks.items() if not value],
        "independent_global_verification": certified,
        "global_hejhal_defect_certified": certified,
        "rung4_certified": False,
        "hejhal_existence_certified": False,
        "dual_certification": False,
        "status": "GREEN" if certified else "RED",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", type=int, default=192)
    parser.add_argument("--result", type=Path, default=ROOT / "track_b_global_hejhal_defect_result.json")
    parser.add_argument("--theorem-defect-definition", type=Path, default=ROOT / "track_b_theorem_defect_definition.json")
    parser.add_argument("--global-partition", type=Path, default=ROOT / "track_b_global_partition_result.json")
    parser.add_argument("--floor-certificate", type=Path, default=ROOT / "track_b_floor_stability_d10_result.json")
    parser.add_argument("--two-cusp-result", type=Path, default=ROOT / "track_b_two_cusp_result.json")
    parser.add_argument("--json-out", type=Path, default=ROOT / "track_b_global_hejhal_verification.json")
    ns = parser.parse_args()
    verification = verify(ns)
    ns.json_out.write_text(json.dumps(verification, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "independent_global_verification": verification["independent_global_verification"],
        "failed_checks": verification["failed_checks"],
    }, indent=2))
    return 0 if verification["independent_global_verification"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
