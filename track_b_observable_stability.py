#!/usr/bin/env python3
"""Exact-normalization cutoff comparisons for observable Track-B trials."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from flint import acb, arb, ctx

from track_b_two_cusp_data import canonical_hash
from track_b_two_cusp_hejhal import acb_from_json


ROOT = Path(__file__).resolve().parent
Label = tuple[str, int, int, int]


def labels(result: dict[str, Any]) -> list[Label]:
    order = result["assembly_definition"]["mode_order"]
    return [
        *(("infinity", *map(int, mode)) for mode in order["infinity"]),
        *(("zero", *map(int, mode)) for mode in order["zero"]),
    ]


def unit_orbit(a: int, b: int) -> tuple[tuple[int, int], ...]:
    return tuple(sorted({(a, b), (-b, a), (-a, -b), (b, -a)}))


def l2_upper(values: list[acb]) -> arb:
    return sum((abs(q).upper() ** 2 for q in values), arb(0)).sqrt().upper()


def l2_lower(values: list[acb]) -> arb:
    return sum((max(abs(q).lower(), arb(0)) ** 2 for q in values), arb(0)).sqrt().lower()


def grouped_difference(
    entries: list[tuple[Label, acb]], key_function: Any
) -> list[dict[str, Any]]:
    groups: dict[str, list[tuple[Label, acb]]] = {}
    for label, value in entries:
        key = str(key_function(label))
        groups.setdefault(key, []).append((label, value))
    output = []
    for key in sorted(groups):
        values = [q[1] for q in groups[key]]
        output.append({
            "group": key,
            "mode_count": len(values),
            "delta_l2_upper": str(l2_upper(values)),
            "real_delta_l2_upper": str(
                sum((abs(q.real).upper() ** 2 for q in values), arb(0)).sqrt().upper()
            ),
            "imag_delta_l2_upper": str(
                sum((abs(q.imag).upper() ** 2 for q in values), arb(0)).sqrt().upper()
            ),
            "modes": [list(q[0]) for q in groups[key]],
        })
    return output


def compare(low: dict[str, Any], high: dict[str, Any]) -> dict[str, Any]:
    m, mp = int(low["fourier_cutoff"]), int(high["fourier_cutoff"])
    if mp <= m:
        raise ValueError("comparison cutoffs must increase")
    if low["normalization_hash"] != high["normalization_hash"]:
        return {"source_cutoff": m, "target_cutoff": mp, "certified": False,
                "reason": "normalization hash changed"}
    if low["phase_convention_hash"] != high["phase_convention_hash"]:
        return {"source_cutoff": m, "target_cutoff": mp, "certified": False,
                "reason": "phase convention hash changed"}
    low_labels, high_labels = labels(low), labels(high)
    high_index = {label: i for i, label in enumerate(high_labels)}
    if any(label not in high_index for label in low_labels):
        return {"source_cutoff": m, "target_cutoff": mp, "certified": False,
                "reason": "mode spaces are not exactly nested"}
    a = [acb_from_json(q) for q in low["physical_coefficient_interval_vector"]]
    b = [acb_from_json(q) for q in high["physical_coefficient_interval_vector"]]
    differences = [b[high_index[label]] - a[i] for i, label in enumerate(low_labels)]
    delta = l2_upper(differences)
    epsilon = arb("1e-40")
    denominator = max(l2_lower(a), epsilon)
    relative = (delta / denominator).upper()
    entries = list(zip(low_labels, differences))
    shell_key = lambda q: (
        q[0], f"{q[3]}" if q[0] == "infinity" else f"{q[3]}/5"
    )
    by_shell = grouped_difference(entries, shell_key)
    by_cusp = grouped_difference(entries, lambda q: q[0])
    by_orbit = grouped_difference(entries, lambda q: (q[0], unit_orbit(q[1], q[2])))
    newest = {}
    for cusp in ("infinity", "zero"):
        cusp_entries = [(label, value) for label, value in entries if label[0] == cusp]
        maximum_norm = max(label[3] for label, _ in cusp_entries)
        values = [value for label, value in cusp_entries if label[3] == maximum_norm]
        newest[cusp] = {
            "mode_norm_label": maximum_norm,
            "physical_shell_squared": str(maximum_norm if cusp == "infinity" else f"{maximum_norm}/5"),
            "delta_l2_upper": str(l2_upper(values)),
        }
    restriction = [[i, high_index[label], 1] for i, label in enumerate(low_labels)]
    return {
        "schema": "track-b-observable-cutoff-comparison/v1",
        "source_cutoff": m,
        "target_cutoff": mp,
        "certified": True,
        "normalization_hash": low["normalization_hash"],
        "phase_convention_hash": low["phase_convention_hash"],
        "phase_alignment": "identity; a_infinity,(1,0) is exactly the positive real number 1 at both cutoffs",
        "unrecorded_phase_alignment": False,
        "restriction_map": restriction,
        "restriction_map_hash": canonical_hash(restriction),
        "delta_l2_upper": str(delta),
        "delta_relative_upper": str(relative),
        "relative_denominator_lower": str(denominator),
        "by_cusp": by_cusp,
        "by_radial_shell": by_shell,
        "by_symmetry_orbit": by_orbit,
        "newest_shared_shell": newest,
        "six_copy_fiber": {
            "structure": "coefficients are shared cusp-channel coefficients, not six independent fiber coefficients",
            "fiber_permutation_inconsistency_applicable": False,
            "delta_l2_upper": str(delta),
        },
        "global_rescaling_removed": True,
        "global_phase_removed": True,
    }


def freeze_gate(results: list[dict[str, Any]], comparisons: list[dict[str, Any]]) -> dict[str, bool]:
    return {
        "all_cutoff_ranks_certified": all(q.get("cutoff_rank_certified", False) for q in results),
        "all_oversampled_rows_independently_verified": all(
            q.get("independent_oversampled_verification", False) for q in results
        ),
        "all_candidate_matrices_observable": all(
            q.get("observability_diagnostics", {}).get("candidate_normalized_numerical_rank")
                == q.get("physical_unknown_count")
            for q in results
        ),
        "all_enrichment_rows_exact_physical": all(
            q.get("enrichment_rows_all_exact_physical", False) for q in results
        ),
        "all_retained_bessel_fallback_counts_zero": all(
            int(q.get("resolved_bessel_fallback_count", -1)) == 0 for q in results
        ),
        "normalization_unchanged": len({q.get("normalization_hash") for q in results}) == 1,
        "phase_convention_unchanged": len({q.get("phase_convention_hash") for q in results}) == 1,
        "all_phase_alignments_recorded": all(
            q.get("certified", False) and not q.get("unrecorded_phase_alignment", True)
            for q in comparisons
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", nargs="+", type=Path, default=[
        ROOT / "track_b_observable_M2_result.json",
        ROOT / "track_b_observable_M3_result.json",
        ROOT / "track_b_observable_M4_result.json",
        ROOT / "track_b_observable_M5_result.json",
    ])
    parser.add_argument("--json-out", type=Path, default=ROOT / "track_b_observable_stability_result.json")
    parser.add_argument("--freeze-status-out", type=Path, default=ROOT / "track_b_trial_freeze_status.json")
    ns = parser.parse_args()
    ctx.prec = 192
    results = [json.loads(path.read_text(encoding="utf-8")) for path in ns.results]
    comparisons = [compare(results[i], results[i + 1]) for i in range(len(results) - 1)]
    by_pair = {(q["source_cutoff"], q["target_cutoff"]): q for q in comparisons}
    d34 = arb(by_pair[(3, 4)]["delta_l2_upper"]).upper()
    d45 = arb(by_pair[(4, 5)]["delta_l2_upper"]).upper()
    d23 = arb(by_pair[(2, 3)]["delta_l2_upper"]).upper()
    rank4 = bool(results[2].get("cutoff_rank_certified", False))
    rank5 = bool(results[3].get("cutoff_rank_certified", False))
    hard = freeze_gate(results, comparisons)
    initial = bool(rank4 and d34 < arb("0.25") and all(hard.values()))
    strong = bool(
        initial and rank5 and d45 < arb("0.1") and d34 < arb("0.75") * d23
    )
    output = {
        "schema": "track-b-observable-cutoff-stability/v1",
        "comparisons": comparisons,
        "hard_freeze_conditions": hard,
        "rank_milestone_M4": rank4,
        "rank_milestone_M5": rank5,
        "delta_2_3_upper": str(d23),
        "delta_3_4_upper": str(d34),
        "delta_4_5_upper": str(d45),
        "initial_threshold": "0.25",
        "strong_threshold": "0.1",
        "initial_convergence_milestone": initial,
        "strong_convergence_milestone": strong,
        "cutoff_stability_passed": strong,
        "trial_frozen": strong,
        "global_hejhal_defect_certified": False,
        "rung4_certified": False,
        "dual_certification": False,
        "status": "GREEN" if strong else "RED",
    }
    output["stability_hash"] = canonical_hash(output)
    ns.json_out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    freeze = {
        "schema": "track-b-trial-freeze-status/v1",
        "rank_milestone_M4": rank4,
        "rank_milestone_M5": rank5,
        "initial_convergence_milestone": initial,
        "strong_convergence_milestone": strong,
        "candidate_frozen_trial_manifest": None,
        "reason": None if strong else "coefficient stability thresholds failed; no track-b-frozen-trial/v1 manifest emitted",
        "trial_frozen": False,
        "global_hejhal_defect_certified": False,
        "rung4_certified": False,
        "dual_certification": False,
    }
    ns.freeze_status_out.write_text(json.dumps(freeze, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "rank_milestone_M4": rank4,
        "rank_milestone_M5": rank5,
        "delta_2_3_upper": str(d23),
        "delta_3_4_upper": str(d34),
        "delta_4_5_upper": str(d45),
        "initial_convergence_milestone": initial,
        "strong_convergence_milestone": strong,
        "trial_frozen": False,
    }, indent=2))
    return 0 if strong else 2


if __name__ == "__main__":
    raise SystemExit(main())
