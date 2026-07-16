#!/usr/bin/env python3
"""Independent verifier for an observable oversampled Track-B cutoff result."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from flint import acb, arb, ctx

from track_b_two_cusp_data import GENERATORS, canonical_hash, cusp_modes
from track_b_two_cusp_hejhal import acb_from_json, matrix_vector
from track_b_two_cusp_verify import slow_row


ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]


def summary(residual: list[acb], indices: list[int]) -> dict[str, arb]:
    bounds = [abs(residual[i]).upper() for i in indices]
    return {
        "component": max(bounds) if bounds else arb(0),
        "l2": sum((q * q for q in bounds), arb(0)).sqrt().upper(),
    }


def endpoint_contains(recomputed: arb, serialized: str) -> bool:
    return bool(recomputed <= arb(serialized).upper())


def verify(result_path: Path, row_path: Path, enrichment_path: Path) -> dict[str, Any]:
    result = json.loads(result_path.read_text(encoding="utf-8"))
    ledger = load_jsonl(row_path)
    enrichment = load_jsonl(enrichment_path)
    ctx.prec = max(128, int(result["arb_bits"]))
    cutoff = int(result["fourier_cutoff"])
    low, high = [arb(q) for q in result["spectral_parameter_input"].split(",")]
    r = low.union(high)
    modes_inf, modes_zero = cusp_modes(cutoff)
    n = len(modes_inf) + len(modes_zero)
    cache: dict[tuple[str, str], acb] = {}
    reconstructed = [
        slow_row(record, modes_inf, modes_zero, r, cache)[0] for record in ledger
    ]
    coefficients = [acb_from_json(q) for q in result["physical_coefficient_interval_vector"]]
    residual = matrix_vector(reconstructed, coefficients)
    selected = [int(q) for q in result["normalization"]["selected_physical_rows"]]
    selected_set = set(selected)
    omitted = [i for i in range(len(ledger)) if i not in selected_set]
    full_summary = summary(residual, list(range(len(ledger))))
    selected_summary = summary(residual, selected)
    omitted_summary = summary(residual, omitted)
    norm_index = int(result["normalization"]["coefficient_index"])
    normalization_ok = bool(
        (coefficients[norm_index] - 1).real.contains(0)
        and (coefficients[norm_index] - 1).imag.contains(0)
    )
    row_ids = [q["row_id"] for q in ledger]
    expected_rows = result["candidate_point_count"] * len(GENERATORS) * 6
    all_blocks = {q["matrix_block"] for q in ledger} == {
        "infinity->infinity", "infinity->zero", "zero->infinity", "zero->zero"
    }
    enrichment_indices = [q["row_index"] for q in enrichment]
    checks = {
        "schema": result.get("schema") == "track-b-observable-two-cusp-hejhal/v1",
        "row_ledger_hash": sha256(row_path) == result.get("full_row_ledger_hash"),
        "enrichment_ledger_hash": sha256(enrichment_path) == result.get("enrichment_ledger_hash"),
        "row_count": len(ledger) == expected_rows == result.get("candidate_physical_row_count"),
        "row_ids_exact_once": len(row_ids) == len(set(row_ids)),
        "row_order": [q["row_index"] for q in ledger] == list(range(len(ledger))),
        "all_four_blocks": all_blocks,
        "selected_exact_once": len(selected) == n - 1 and len(selected) == len(selected_set),
        "omitted_exact_complement": len(omitted) == len(ledger) - (n - 1),
        "enrichment_matches_selected": enrichment_indices == selected,
        "all_enrichment_rows_physical": all(
            q.get("physical_row") and not q.get("synthetic_row") for q in enrichment
        ),
        "coefficient_count": len(coefficients) == n,
        "coefficient_hash": canonical_hash(result["physical_coefficient_interval_vector"])
            == result.get("coefficient_hash"),
        "normalization_exact": normalization_ok,
        "full_component_endpoint": endpoint_contains(
            full_summary["component"], result["full_physical_residual"]["component_upper"]
        ),
        "full_l2_endpoint": endpoint_contains(
            full_summary["l2"], result["full_physical_residual"]["l2_upper"]
        ),
        "selected_component_endpoint": endpoint_contains(
            selected_summary["component"], result["selected_physical_residual"]["component_upper"]
        ),
        "selected_l2_endpoint": endpoint_contains(
            selected_summary["l2"], result["selected_physical_residual"]["l2_upper"]
        ),
        "omitted_component_endpoint": endpoint_contains(
            omitted_summary["component"], result["omitted_physical_residual"]["component_upper"]
        ),
        "omitted_l2_endpoint": endpoint_contains(
            omitted_summary["l2"], result["omitted_physical_residual"]["l2_upper"]
        ),
        "phase_alias_free": result["phase_alias_audit"]["no_unexplained_exact_phase_aliases"],
        "candidate_rank_full": result["observability_diagnostics"]["candidate_normalized_numerical_rank"] == n,
        "selected_rank_full": result["observability_diagnostics"]["selected_numerical_rank"] == n,
        "verified_interval_solve": bool(result.get("verified_solve_certified", False)),
        "cutoff_rank_certified": bool(result.get("cutoff_rank_certified", False)),
        "bessel_fallback_zero": int(result.get("resolved_bessel_fallback_count", -1)) == 0,
        "trial_not_frozen": result.get("trial_frozen") is False,
        "rung4_false": result.get("rung4_certified") is False,
        "dual_false": result.get("dual_certification") is False,
    }
    verified = all(checks.values())
    return {
        "schema": "track-b-observable-two-cusp-independent-verification/v1",
        "fourier_cutoff": cutoff,
        "result_path": str(result_path.resolve()),
        "result_sha256": sha256(result_path),
        "checks": checks,
        "failed_checks": [q for q, value in checks.items() if not value],
        "independent_bessel_cache_entries": len(cache),
        "independent_oversampled_verification": verified,
        "cutoff_rank_certified": verified,
        "trial_frozen": False,
        "global_hejhal_defect_certified": False,
        "rung4_certified": False,
        "dual_certification": False,
        "status": "GREEN" if verified else "RED",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--row-ledger", type=Path, required=True)
    parser.add_argument("--enrichment-ledger", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    ns = parser.parse_args()
    output = verify(ns.result, ns.row_ledger, ns.enrichment_ledger)
    ns.json_out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "fourier_cutoff": output["fourier_cutoff"],
        "independent_oversampled_verification": output["independent_oversampled_verification"],
        "failed_checks": output["failed_checks"],
        "trial_frozen": False,
    }, indent=2))
    return 0 if output["independent_oversampled_verification"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
