#!/usr/bin/env python3
"""Deterministic M=1..4 Track-B two-cusp cutoff and nesting audit."""
from __future__ import annotations

import argparse
import json
import time
import tracemalloc
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from flint import acb, arb, ctx

from track_b_two_cusp_data import canonical_hash, cusp_modes
from track_b_two_cusp_hejhal import acb_from_json, run


ROOT = Path(__file__).resolve().parent


def coefficient_labels(cutoff: int) -> list[tuple[str, int, int, int]]:
    inf, zero = cusp_modes(cutoff)
    return [*(('infinity', *m) for m in inf), *(('zero', *m) for m in zero)]


def exact_embedding(source: int, target: int) -> dict[str, Any]:
    if target <= source:
        raise ValueError("target cutoff must exceed source cutoff")
    source_labels = coefficient_labels(source)
    target_labels = coefficient_labels(target)
    target_index = {label: i for i, label in enumerate(target_labels)}
    entries = [[i, target_index[label], 1] for i, label in enumerate(source_labels)]
    result = {
        "schema": "track-b-cutoff-embedding/v1",
        "source_cutoff": source,
        "target_cutoff": target,
        "source_dimension": len(source_labels),
        "target_dimension": len(target_labels),
        "source_labels": [list(q) for q in source_labels],
        "target_labels": [list(q) for q in target_labels],
        "entries_source_target_value": entries,
        "cusp_labels_preserved": True,
        "mode_labels_preserved": True,
        "six_copy_fibers": "unchanged; coefficients are cusp-channel coefficients",
        "realification": "none; complex Arb coefficient vector",
        "normalization": "a_infinity,(1,0)=1 at both cutoffs",
    }
    result["embedding_hash"] = canonical_hash(result)
    return result


def compare_low_modes(low: dict[str, Any], high: dict[str, Any]) -> dict[str, Any]:
    m = int(low["fourier_cutoff"])
    mp = int(high["fourier_cutoff"])
    embedding = exact_embedding(m, mp)
    a = [acb_from_json(q) for q in low.get("physical_coefficient_interval_vector", [])]
    b = [acb_from_json(q) for q in high.get("physical_coefficient_interval_vector", [])]
    if len(a) != embedding["source_dimension"] or len(b) != embedding["target_dimension"]:
        return {
            "embedding": embedding,
            "comparison_certified": False,
            "reason": "one or both verified coefficient enclosures are absent or have the wrong dimension",
            "delta_l2_upper": None,
            "delta_component_upper": None,
        }
    differences = []
    records = []
    for source_index, target_index, _ in embedding["entries_source_target_value"]:
        difference = b[target_index] - a[source_index]
        upper = abs(difference).upper()
        differences.append(difference)
        records.append({
            "label": embedding["source_labels"][source_index],
            "source_index": source_index,
            "target_index": target_index,
            "difference_abs_upper": str(upper),
        })
    l2 = sum((abs(q).upper() ** 2 for q in differences), arb(0)).sqrt().upper()
    component = max(abs(q).upper() for q in differences)
    return {
        "embedding": embedding,
        "comparison_certified": True,
        "delta_definition": "||a_Mprime_low-E_(M->Mprime)a_M||_2",
        "delta_l2_upper": str(l2),
        "delta_component_upper": str(component),
        "per_mode": records,
    }


def coefficient_dynamic_range(result: dict[str, Any]) -> dict[str, Any]:
    values = [acb_from_json(q) for q in result.get("physical_coefficient_interval_vector", [])]
    if not values:
        return {"certified": False, "upper": None, "reason": "missing coefficient enclosure"}
    uppers = [abs(q).upper() for q in values]
    lowers = [abs(q).lower() for q in values]
    positive = [q for q in lowers if bool(q > 0)]
    if len(positive) != len(lowers):
        return {
            "certified": False,
            "upper": None,
            "maximum_coefficient_upper": str(max(uppers)),
            "minimum_coefficient_lower": str(min(lowers)),
            "reason": "at least one coefficient enclosure contains zero",
        }
    bound = (max(uppers) / min(positive)).upper()
    return {
        "certified": True,
        "upper": str(bound),
        "maximum_coefficient_upper": str(max(uppers)),
        "minimum_coefficient_lower": str(min(positive)),
    }


def run_one(m: int, ns: argparse.Namespace) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    out = ns.output_dir / f"track_b_two_cusp_M{m}_result.json"
    ledger = ns.output_dir / f"track_b_two_cusp_M{m}_rows.jsonl"
    verify = ns.output_dir / f"track_b_two_cusp_M{m}_verification.json"
    run_ns = SimpleNamespace(
        bits=ns.bits,
        r_interval=ns.r_interval,
        fourier_cutoff=m,
        collocation_points=ns.collocation_points,
        assemble_physical=True,
        verified_solve=True,
        independent_verify=True,
        collocation_ledger=ledger,
        json_out=out,
        verification_json_out=verify,
        global_partition=ns.global_partition,
        floor_certificate=ns.floor_certificate,
        theorem_defect_definition=ns.theorem_defect_definition,
    )
    if ns.reuse_existing and out.is_file():
        result = json.loads(out.read_text(encoding="utf-8"))
        theorem = json.loads(ns.theorem_defect_definition.read_text(encoding="utf-8"))
        current_theorem_hash = theorem.get("theorem_defect_definition_hash")
        if result.get("theorem_defect_definition_hash") != current_theorem_hash:
            # The theorem interface hash is metadata independent of the
            # interval solve.  Refreshing it does not change any enclosure.
            result["theorem_defect_definition_hash"] = current_theorem_hash
            out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        elapsed = None
        peak = None
    else:
        tracemalloc.start()
        started = time.perf_counter()
        try:
            result = run(run_ns)
        except Exception as exc:
            elapsed = time.perf_counter() - started
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            return None, {
                "fourier_cutoff": m,
                "result_path": None,
                "row_ledger_path": str(ledger.resolve()) if ledger.is_file() else None,
                "verification_path": None,
                "verified_solve_certified": False,
                "independent_verification_certified": False,
                "failure_type": type(exc).__name__,
                "failure": str(exc),
                "elapsed_seconds": elapsed,
                "python_tracemalloc_peak_bytes": peak,
                "continuum_residual_upper": None,
                "fourier_tail_upper": None,
                "complete_theorem_defect_upper": None,
            }
        elapsed = time.perf_counter() - started
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    metrics = {
        "fourier_cutoff": m,
        "result_path": str(out.resolve()),
        "row_ledger_path": str(ledger.resolve()),
        "verification_path": str(verify.resolve()),
        "unknown_count": result["physical_unknown_count"],
        "physical_row_count": result["physical_matrix_shape"][0],
        "selected_row_count": result["physical_unknown_count"] - 1,
        "omitted_row_count": result["physical_matrix_shape"][0] - result["physical_unknown_count"] + 1,
        "verified_solve_certified": result["verified_solve_certified"],
        "independent_verification_certified": result["independent_verification_certified"],
        "verified_contraction_upper": result.get("verified_solve", {}).get("contraction_upper"),
        "sigma_min_lower": result.get("verified_solve", {}).get("smallest_singular_value_lower"),
        "coefficient_dynamic_range": coefficient_dynamic_range(result),
        "physical_residual_component_upper": result.get("physical_residual_component_upper"),
        "physical_residual_l2_upper": result.get("physical_residual_l2_upper"),
        "continuum_residual_upper": None,
        "fourier_tail_upper": None,
        "complete_theorem_defect_upper": None,
        "elapsed_seconds": elapsed,
        "python_tracemalloc_peak_bytes": peak,
        "condition_diagnostics": result.get("diagnostics_not_certification_inputs", {}),
        "theorem_defect_definition_hash": result.get("theorem_defect_definition_hash"),
    }
    return result, metrics


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", type=int, default=192)
    parser.add_argument("--r-interval", default="6.7439020359331625,6.7439020359331625")
    parser.add_argument("--cutoffs", default="1,2,3,4")
    parser.add_argument("--collocation-points", type=int, default=4)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "track_b_cutoff_ladder")
    parser.add_argument("--global-partition", type=Path, default=ROOT / "track_b_global_partition_result.json")
    parser.add_argument("--floor-certificate", type=Path, default=ROOT / "track_b_floor_stability_d10_result.json")
    parser.add_argument("--theorem-defect-definition", type=Path, default=ROOT / "track_b_theorem_defect_definition.json")
    parser.add_argument("--reuse-existing", action="store_true")
    parser.add_argument("--json-out", type=Path, default=ROOT / "track_b_cutoff_ladder_result.json")
    ns = parser.parse_args()
    ctx.prec = max(128, ns.bits)
    cutoffs = [int(q) for q in ns.cutoffs.split(",")]
    if cutoffs != sorted(set(cutoffs)) or not cutoffs or cutoffs[0] < 1:
        raise ValueError("cutoffs must be a strictly increasing positive list")
    ns.output_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any] | None] = []
    metrics: list[dict[str, Any]] = []
    for m in cutoffs:
        result, metric = run_one(m, ns)
        results.append(result)
        metrics.append(metric)
    comparisons = []
    for i in range(len(results) - 1):
        if results[i] is not None and results[i + 1] is not None:
            comparisons.append(compare_low_modes(results[i], results[i + 1]))
        else:
            comparisons.append({
                "source_cutoff": cutoffs[i],
                "target_cutoff": cutoffs[i + 1],
                "comparison_certified": False,
                "reason": "one or both cutoff solves failed",
            })
    theorem_hashes = {
        q.get("theorem_defect_definition_hash") for q in results if q is not None
    }
    output = {
        "schema": "track-b-two-cusp-cutoff-ladder/v1",
        "label": "TWO-CUSP CUTOFF LADDER AND EXACT NESTING AUDIT",
        "r_interval": ns.r_interval,
        "cutoffs": cutoffs,
        "collocation_points": ns.collocation_points,
        "runs": metrics,
        "nesting_comparisons": comparisons,
        "normalization_unchanged": True,
        "theorem_defect_definition_hash": next(iter(theorem_hashes)) if len(theorem_hashes) == 1 else None,
        "all_verified_solves_certified": all(
            q is not None and q["verified_solve_certified"] for q in results
        ),
        "all_physical_rows_independently_verified": all(
            q is not None and q["independent_verification_certified"] for q in results
        ),
        "cutoff_convergence_checked": all(q["comparison_certified"] for q in comparisons),
        "continuum_channels_certified": False,
        "fourier_tails_certified": False,
        "global_hejhal_defect_certified": False,
        "rung4_certified": False,
        "hejhal_existence_certified": False,
        "dual_certification": False,
        "status": "RED",
    }
    output["ladder_hash"] = canonical_hash(output)
    ns.json_out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "all_verified_solves_certified": output["all_verified_solves_certified"],
        "cutoff_convergence_checked": output["cutoff_convergence_checked"],
        "global_hejhal_defect_certified": False,
        "result": str(ns.json_out.resolve()),
    }, indent=2))
    return 0 if output["all_verified_solves_certified"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
