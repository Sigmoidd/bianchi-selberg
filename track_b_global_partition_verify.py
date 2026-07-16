#!/usr/bin/env python3
"""Independent fail-closed verifier for a global Track-B partition ledger."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from flint import arb, ctx


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()]


def _arb(value: Any) -> arb:
    try:
        out = arb(str(value))
    except (ValueError, TypeError):
        return arb("nan")
    return out


def verify(result: dict[str, Any], records: list[dict[str, Any]]) -> dict[str, Any]:
    dims = tuple(int(q) for q in result.get("subdivision", []))
    expected = dims[0] * dims[1] * dims[2] if len(dims) == 3 else -1
    ids = [str(q.get("cell_id")) for q in records]
    expected_ids = ({f"{i},{j},{k}" for i in range(dims[0])
                     for j in range(dims[1]) for k in range(dims[2])}
                    if len(dims) == 3 else set())
    exact_once = (
        len(records) == expected and len(ids) == len(set(ids))
        and set(ids) == expected_ids
    )
    hashes_match = bool(records) and all(
        q.get("partition_definition_hash") == result.get("partition_definition_hash")
        and q.get("geometry_certificate_dependency") == result.get("geometry_incidence_hash")
        and q.get("weight_formulas_hash") == result.get("weight_formulas_hash")
        and q.get("configuration_hash") == result.get("configuration_hash")
        for q in records
    )
    minimum = min((_arb(q.get("Phi_lower")).lower() for q in records),
                  default=arb("-inf"))
    maximum_gradient = max(
        (_arb(q.get("maximum_grad_chi_upper")).upper() for q in records),
        default=arb("+inf"),
    )
    maximum_laplacian = max(
        (_arb(q.get("maximum_Delta_chi_upper")).upper() for q in records),
        default=arb("+inf"),
    )
    fallbacks = sum(int(q.get("fallback_count", -1)) for q in records)
    cell_flags = bool(records) and all(
        all(q.get("coverage_flags", {}).values())
        and q.get("floor_weight_consistency", False)
        and _arb(q.get("sum_chi_interval", "nan")).contains(1)
        and _arb(q.get("partition_deviation_upper", "nan")).contains(0)
        and bool(q.get("geometry_region_type"))
        and _arb(q.get("maximum_grad_chi_upper")).is_finite()
        and _arb(q.get("maximum_Delta_chi_upper")).is_finite()
        for q in records
    )
    details_complete = bool(records) and all(
        q.get("per_active_weight_enclosures_recorded", False)
        and len(q.get("per_active_weight_enclosures", []))
            == len(q.get("active_normalized_weights", []))
        and all(
            _arb(w.get(field)).is_finite()
            for w in q.get("per_active_weight_enclosures", [])
            for field in (
                "phi_interval", "grad_phi_hyperbolic_certified_upper",
                "Delta_phi_certified_interval", "chi_interval",
                "grad_chi_hyperbolic_certified_upper",
                "Delta_chi_certified_interval", "abs_Delta_chi_certified_upper",
            )
        )
        for q in records
    )
    endpoint_agreement = bool(
        _arb(result.get("minimum_denominator_lower")).contains(minimum)
        and _arb(result.get("maximum_gradient_upper")).contains(maximum_gradient)
        and _arb(result.get("maximum_laplacian_upper")).contains(maximum_laplacian)
    )
    stabilizer = result.get("stabilizer_certificate", {})
    stabilizers_exact = bool(
        stabilizer.get("certified", False)
        and stabilizer.get("global_stratum_group_reindexing_exact", False)
        and stabilizer.get("vertex_group_reindexing_exact")
        and all(stabilizer.get("vertex_group_reindexing_exact", {}).values())
        and not stabilizer.get("sampled_equality_used", True)
    )
    reconstructed_payload = "".join(
        json.dumps(q, sort_keys=True) + "\n" for q in records
    ).encode("utf-8")
    audit_hash_match = bool(
        result.get("audit_records_sha256")
        and hashlib.sha256(reconstructed_payload).hexdigest()
            == result.get("audit_records_sha256")
    )
    required_result_flags = (
        "coverage_certified", "denominator_positive_certified",
        "partition_sum_certified", "stabilizer_averages_certified",
        "weight_gradients_certified", "weight_laplacians_certified",
        "floor_weight_consistency_certified", "geometry_incidence_certified",
        "global_partition_certified", "global_weight_bounds_certified",
        "stability_check_passed",
    )
    flags_green = all(result.get(name, False) for name in required_result_flags)
    verified = bool(
        exact_once and hashes_match and minimum > 0 and fallbacks == 0
        and cell_flags and details_complete and endpoint_agreement
        and stabilizers_exact and audit_hash_match and flags_green
        and result.get("rung4_certified") is False
    )
    return {
        "schema": "track-b-global-partition-verification/v1",
        "verified": verified,
        "cell_count": len(records),
        "expected_cell_count": expected,
        "exact_once": exact_once,
        "hashes_match": hashes_match,
        "audit_hash_match": audit_hash_match,
        "reconstructed_minimum_denominator_lower": str(minimum),
        "reconstructed_maximum_gradient_upper": str(maximum_gradient),
        "reconstructed_maximum_laplacian_upper": str(maximum_laplacian),
        "fallback_count": fallbacks,
        "cell_flags_green": cell_flags,
        "per_weight_enclosures_complete": details_complete,
        "stabilizers_exact": stabilizers_exact,
        "endpoint_agreement": endpoint_agreement,
        "result_flags_green": flags_green,
        "global_partition_certified": verified,
        "global_weight_bounds_certified": verified,
        "rung4_certified": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--audit-jsonl", type=Path, default=None)
    parser.add_argument("--json-out", type=Path,
                        default=Path("track_b_global_partition_verification.json"))
    parser.add_argument("--bits", type=int, default=192)
    ns = parser.parse_args()
    ctx.prec = max(128, ns.bits)
    result = json.loads(ns.result.read_text(encoding="utf-8"))
    audit = ns.audit_jsonl or Path(result.get("audit_jsonl", ""))
    if not audit.is_file():
        output = {
            "verified": False,
            "reason": "missing per-cell audit ledger",
            "global_partition_certified": False,
            "global_weight_bounds_certified": False,
            "rung4_certified": False,
        }
    else:
        output = verify(result, load_jsonl(audit))
    ns.json_out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))
    return 0 if output.get("verified", False) else 2


if __name__ == "__main__":
    raise SystemExit(main())
