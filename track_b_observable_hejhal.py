#!/usr/bin/env python3
"""Adaptive exact-row enrichment for observable two-cusp Hejhal systems."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path
from typing import Any

from flint import acb, arb, ctx
import numpy as np

from track_b_collocation_family import (
    build_family, exact_candidate_points, exact_h3_action, f, phase_audit,
)
from track_b_two_cusp_data import GENERATORS, canonical_hash, cusp_data, mode_frequency
from track_b_two_cusp_hejhal import (
    BLOCKS, ROUNDING_SAFETY, ValidatedWhittaker, acb_from_json, acb_json,
    assemble_physical_rows_from_points, exact_power_two_scaling, matrix_vector,
    midpoint_matrix, parse_r_interval, residual_summary, scale_system,
    verified_contraction_solve,
)
from track_b_two_cusp_verify import slow_row


ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalization_index(modes_inf: list[tuple[int, int, int]]) -> int:
    return next(i for i, mode in enumerate(modes_inf) if mode[:2] == (1, 0))


def numerical_rank(matrix: np.ndarray) -> tuple[int, float, np.ndarray]:
    singular = np.linalg.svd(matrix, compute_uv=False)
    tolerance = max(matrix.shape) * np.finfo(float).eps * max(float(singular[0]), 1.0)
    return int(np.sum(singular > tolerance)), tolerance, singular


def least_observed_vector(basis: np.ndarray) -> tuple[np.ndarray, float]:
    _u, singular, vh = np.linalg.svd(basis, full_matrices=True)
    vector = vh[-1].conjugate()
    value = float(singular[-1]) if basis.shape[0] >= basis.shape[1] else 0.0
    return vector, value


def adaptive_select(
    physical: np.ndarray, ledger: list[dict[str, Any]], norm_index: int
) -> tuple[list[int], list[dict[str, Any]], dict[str, Any]]:
    row_count, n = physical.shape
    normal = np.zeros(n, dtype=np.complex128)
    normal[norm_index] = 1
    basis = normal.reshape(1, -1)
    chosen: list[int] = []
    events: list[dict[str, Any]] = []

    def append_best(candidates: list[int], reason: str) -> None:
        nonlocal basis
        vector, least = least_observed_vector(basis)
        scored = []
        for index in candidates:
            if index in chosen:
                continue
            action = float(abs(physical[index] @ vector))
            scored.append((action, str(ledger[index]["row_id"]), index))
        if not scored:
            raise ArithmeticError(f"no unused candidate row for {reason}")
        action, _row_id, index = max(scored, key=lambda q: (q[0], tuple(-ord(c) for c in q[1])))
        old_rank = numerical_rank(basis)[0]
        candidate = np.vstack([basis, physical[index]])
        new_rank = numerical_rank(candidate)[0]
        if new_rank <= old_rank:
            # The single least vector may be numerically ambiguous inside a
            # multidimensional nullspace.  Fall back deterministically to the
            # row with largest total projection onto the full nullspace.
            _u, s, vh = np.linalg.svd(basis, full_matrices=True)
            rank = numerical_rank(basis)[0]
            null = vh[rank:].conjugate().T
            rescored = []
            for j in candidates:
                if j in chosen:
                    continue
                score = float(np.linalg.norm(physical[j] @ null))
                rescored.append((score, str(ledger[j]["row_id"]), j))
            score, _rid, index = max(rescored, key=lambda q: (q[0], tuple(-ord(c) for c in q[1])))
            action = score
            candidate = np.vstack([basis, physical[index]])
            new_rank = numerical_rank(candidate)[0]
        if new_rank <= old_rank:
            raise ArithmeticError(f"candidate family cannot enrich rank at {old_rank}/{n}")
        chosen.append(index)
        basis = candidate
        events.append({
            "step": len(chosen),
            "row_index": index,
            "row_id": ledger[index]["row_id"],
            "matrix_block": ledger[index]["matrix_block"],
            "group_element": ledger[index]["group_element"],
            "source_copy": ledger[index]["source_copy"],
            "least_singular_value_before_diagnostic": least,
            "action_on_least_observed_vector": action,
            "rank_after": new_rank,
            "reason": reason,
            "physical_row": True,
            "synthetic_row": False,
        })

    for block in BLOCKS:
        append_best(
            [i for i, record in enumerate(ledger) if record["matrix_block"] == block],
            f"seed exact physical block {block}",
        )
    while len(chosen) < n - 1:
        append_best(list(range(row_count)), "least-observed-vector enrichment")

    selected = np.vstack([physical[chosen], normal])
    rank, tolerance, singular = numerical_rank(selected)
    full_with_norm = np.vstack([physical, normal])
    full_rank, full_tolerance, full_singular = numerical_rank(full_with_norm)
    diagnostics = {
        "candidate_normalized_shape": list(full_with_norm.shape),
        "candidate_normalized_numerical_rank": full_rank,
        "candidate_rank_tolerance": full_tolerance,
        "candidate_smallest_singular_value": float(full_singular[-1]),
        "selected_normalized_shape": list(selected.shape),
        "selected_numerical_rank": rank,
        "selected_rank_tolerance": tolerance,
        "selected_smallest_singular_value": float(singular[-1]),
        "selected_condition_estimate": float(singular[0] / singular[-1]),
        "meaningful_singular_gap": bool(singular[-1] > 100 * tolerance),
        "all_four_blocks_selected": {ledger[i]["matrix_block"] for i in chosen} == set(BLOCKS),
    }
    return chosen, events, diagnostics


def subset_residual_summary(
    residual: list[acb], indices: list[int], ledger: list[dict[str, Any]]
) -> dict[str, Any]:
    if not indices:
        return {"row_count": 0, "component_upper": "0", "l2_upper": "0", "worst_row": None}
    bounds = [(ROUNDING_SAFETY * abs(residual[i]).upper()).upper() for i in indices]
    worst_local = max(range(len(indices)), key=lambda k: float(bounds[k]))
    worst = indices[worst_local]
    l2 = sum((q * q for q in bounds), arb(0)).sqrt().upper()
    return {
        "row_count": len(indices),
        "component_upper": str(bounds[worst_local]),
        "l2_upper": str(l2),
        "worst_row_index": worst,
        "worst_row_id": ledger[worst]["row_id"],
        "worst_cusp_pair": ledger[worst]["matrix_block"],
        "worst_pairing_relation": ledger[worst]["group_element"],
    }


def independent_rows(
    ledger: list[dict[str, Any]], modes_inf: list[tuple[int, int, int]],
    modes_zero: list[tuple[int, int, int]], r: arb,
) -> tuple[list[list[acb]], int]:
    cache: dict[tuple[str, str], acb] = {}
    rows = [slow_row(record, modes_inf, modes_zero, r, cache)[0] for record in ledger]
    return rows, len(cache)


def radial_profile_audit(
    points: list[dict[str, Any]], modes_inf: list[tuple[int, int, int]],
    modes_zero: list[tuple[int, int, int]], backend: ValidatedWhittaker,
) -> dict[str, Any]:
    per_cusp = {}
    overall_max = 0.0
    overall_worst = None
    for cusp, modes in (("infinity", modes_inf), ("zero", modes_zero)):
        exact_heights = set()
        for record in points:
            point = (f(record["x1"]), f(record["x2"]), f(record["y"]))
            exact_heights.add(point[2])
            for matrix in GENERATORS.values():
                exact_heights.add(exact_h3_action(matrix, point)[2])
        heights = [arb(q.numerator) / q.denominator for q in sorted(exact_heights)]
        shells: dict[str, arb] = {}
        for mode in modes:
            magnitude = mode_frequency(cusp, mode)[2]
            label = str(mode[2] if cusp == "infinity" else f"{mode[2]}/5")
            shells.setdefault(label, magnitude)
        profiles = {}
        for label, magnitude in shells.items():
            profiles[label] = np.array([
                float(backend.radial(magnitude, height)[0].real.mid()) for height in heights
            ])
        maximum = 0.0
        worst = None
        labels = sorted(profiles)
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                a, b = profiles[labels[i]], profiles[labels[j]]
                denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
                coherence = 0.0 if denominator == 0 else float(abs(np.vdot(a, b)) / denominator)
                if coherence > maximum:
                    maximum = coherence
                    worst = [labels[i], labels[j]]
        per_cusp[cusp] = {
            "deduplicated_true_height_count": len(heights),
            "radial_shell_count": len(shells),
            "maximum_profile_coherence_diagnostic": maximum,
            "worst_shell_pair": worst,
            "heights": [str(q) for q in heights],
        }
        if maximum > overall_max:
            overall_max = maximum
            overall_worst = {"cusp": cusp, "shells": worst}
    return {
        "schema": "track-b-radial-profile-observability/v1",
        "per_cusp": per_cusp,
        "maximum_profile_coherence_diagnostic": overall_max,
        "worst_pair": overall_worst,
        "source_and_true_transformed_heights_recorded": True,
        "profiles_use_deduplicated_true_heights": True,
        "diagnostic_only": True,
    }


def run_cutoff(ns: argparse.Namespace) -> dict[str, Any]:
    ctx.prec = max(128, ns.bits)
    started = time.perf_counter()
    r = parse_r_interval(ns.r_interval)
    points = exact_candidate_points(ns.candidate_points)
    family = build_family(ns.candidate_points, max(ns.fourier_cutoff, 5))
    phase = phase_audit(ns.fourier_cutoff, points)
    if not phase["no_unexplained_exact_phase_aliases"]:
        raise ArithmeticError("candidate family has an unexplained exact phase alias")
    backend = ValidatedWhittaker(r)
    rows, ledger, modes_inf, modes_zero = assemble_physical_rows_from_points(
        ns.fourier_cutoff, points, backend, ns.row_ledger
    )
    radial_audit = radial_profile_audit(points, modes_inf, modes_zero, backend)
    physical_mid = midpoint_matrix(rows)
    norm_index = normalization_index(modes_inf)
    selected, events, observability = adaptive_select(physical_mid, ledger, norm_index)
    ns.enrichment_ledger.write_text(
        "".join(json.dumps(q, sort_keys=True) + "\n" for q in events), encoding="utf-8"
    )
    n = len(rows[0])
    normalized_matrix = [rows[i] for i in selected]
    normalization_row = [acb(0) for _ in range(n)]
    normalization_row[norm_index] = acb(1)
    normalized_matrix.append(normalization_row)
    rhs = [acb(0) for _ in selected] + [acb(1)]
    dr, dc = exact_power_two_scaling(normalized_matrix)
    scaled_matrix, scaled_rhs = scale_system(normalized_matrix, rhs, dr, dc)
    solve = verified_contraction_solve(scaled_matrix, scaled_rhs)
    coefficients: list[acb] = []
    if solve.get("certified", False):
        scaled_coeff = [acb_from_json(q) for q in solve["coefficient_enclosure"]]
        coefficients = [dc[j] * scaled_coeff[j] for j in range(n)]
    coefficient_json = [acb_json(q) for q in coefficients]
    coefficients = [acb_from_json(q) for q in coefficient_json]

    residual = matrix_vector(rows, coefficients) if coefficients else []
    selected_set = set(selected)
    omitted = [i for i in range(len(rows)) if i not in selected_set]
    full_summary = residual_summary(residual, ledger) if residual else None
    selected_summary = subset_residual_summary(residual, selected, ledger) if residual else None
    omitted_summary = subset_residual_summary(residual, omitted, ledger) if residual else None

    independent_matrix, independent_bessel_cache = independent_rows(
        ledger, modes_inf, modes_zero, r
    ) if coefficients else ([], 0)
    independent_residual = matrix_vector(independent_matrix, coefficients) if coefficients else []
    independent_summary = residual_summary(independent_residual, ledger) if coefficients else None
    independent_agrees = bool(
        independent_summary and full_summary
        and arb(independent_summary["component_upper"]).upper()
            <= ROUNDING_SAFETY * arb(full_summary["component_upper"]).upper()
        and arb(independent_summary["l2_upper"]).upper()
            <= ROUNDING_SAFETY * arb(full_summary["l2_upper"]).upper()
    )

    midpoint_fro = float(np.linalg.norm(physical_mid))
    radius_fro = math.sqrt(sum(
        float(value.real.rad()) ** 2 + float(value.imag.rad()) ** 2
        for row in rows for value in row
    ))
    bessel = backend.summary()
    theorem = json.loads(ns.theorem_defect_definition.read_text(encoding="utf-8"))
    assembly_definition = {
        "schema": "track-b-observable-two-cusp-assembly/v1",
        "r_interval": ns.r_interval,
        "fourier_cutoff": ns.fourier_cutoff,
        "collocation_family_hash": family["collocation_family_hash"],
        "candidate_point_count": ns.candidate_points,
        "candidate_physical_row_count": len(rows),
        "selected_physical_rows": selected,
        "normalization": {"cusp": "infinity", "mode": [1, 0], "value": "1"},
        "mode_order": {"infinity": modes_inf, "zero": modes_zero},
        "cusp_data_hash": cusp_data()["cusp_data_hash"],
    }
    cutoff_rank = bool(
        solve.get("certified", False)
        and observability["candidate_normalized_numerical_rank"] == n
        and observability["selected_numerical_rank"] == n
        and observability["meaningful_singular_gap"]
    )
    result = {
        "schema": "track-b-observable-two-cusp-hejhal/v1",
        "label": "OBSERVABLE OVERSAMPLED TWO-CUSP HEJHAL SYSTEM",
        "arb_bits": int(ctx.prec),
        "spectral_parameter_input": ns.r_interval,
        "spectral_parameter_interval": str(r),
        "fourier_cutoff": ns.fourier_cutoff,
        "physical_unknown_count": n,
        "candidate_point_count": ns.candidate_points,
        "candidate_physical_row_count": len(rows),
        "selected_square_row_count": n,
        "selected_physical_row_count": len(selected),
        "omitted_verification_row_count": len(omitted),
        "all_four_cusp_blocks_present": {q["matrix_block"] for q in ledger} == set(BLOCKS),
        "phase_alias_audit": phase,
        "transformed_height_audit": family["transformed_height_audit"],
        "radial_profile_observability": radial_audit,
        "observability_diagnostics": observability,
        "normalization": {
            "row": "a_infinity,(1,0)=1",
            "coefficient_index": norm_index,
            "target": "1",
            "phase_convention": "normalizing coefficient is the exact positive real number 1",
            "selected_physical_rows": selected,
        },
        "normalization_hash": canonical_hash({"cusp": "infinity", "mode": [1, 0], "value": "1"}),
        "phase_convention_hash": canonical_hash({"a_infinity_(1,0)": "exact positive real 1"}),
        "scaling": {
            "D_r": [str(q) for q in dr], "D_c": [str(q) for q in dc],
            "formula": "B_tilde=D_r B D_c; a=D_c a_tilde",
            "regularization": None,
        },
        "verified_solve": solve,
        "verified_solve_hash": canonical_hash(solve),
        "verified_solve_certified": bool(solve.get("certified", False)),
        "cutoff_rank_certified": cutoff_rank,
        "certified_sigma_min_lower": solve.get("smallest_singular_value_lower"),
        "contraction_upper": solve.get("contraction_upper"),
        "interval_radius_midpoint_ratio": radius_fro / max(midpoint_fro, 1e-300),
        "condition_estimate": observability["selected_condition_estimate"],
        "physical_coefficient_interval_vector": coefficient_json,
        "coefficient_hash": canonical_hash(coefficient_json),
        "full_physical_residual": full_summary,
        "selected_physical_residual": selected_summary,
        "omitted_physical_residual": omitted_summary,
        "independent_oversampled_residual": independent_summary,
        "independent_oversampled_verification": independent_agrees,
        "independent_bessel_cache_entries": independent_bessel_cache,
        "resolved_bessel_fallback_count": bessel["failed"],
        "bessel_backend": bessel,
        "assembly_definition": assembly_definition,
        "assembly_hash": canonical_hash(assembly_definition),
        "collocation_family_hash": family["collocation_family_hash"],
        "selected_row_hash": canonical_hash(selected),
        "full_row_ledger": str(ns.row_ledger.resolve()),
        "full_row_ledger_hash": sha256(ns.row_ledger),
        "enrichment_ledger": str(ns.enrichment_ledger.resolve()),
        "enrichment_ledger_hash": sha256(ns.enrichment_ledger),
        "enrichment_rows_all_exact_physical": all(
            q["physical_row"] and not q["synthetic_row"] for q in events
        ),
        "cusp_data_hash": cusp_data()["cusp_data_hash"],
        "theorem_defect_definition_hash": theorem["theorem_defect_definition_hash"],
        "elapsed_seconds": time.perf_counter() - started,
        "cutoff_stability_passed": False,
        "trial_frozen": False,
        "global_hejhal_defect_certified": False,
        "rung4_certified": False,
        "dual_certification": False,
        "status": "RANK_GREEN_STABILITY_PENDING" if cutoff_rank and independent_agrees else "RED",
    }
    ns.json_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "fourier_cutoff": ns.fourier_cutoff,
        "candidate_rank": observability["candidate_normalized_numerical_rank"],
        "selected_rank": observability["selected_numerical_rank"],
        "cutoff_rank_certified": cutoff_rank,
        "certified_sigma_min_lower": result["certified_sigma_min_lower"],
        "contraction_upper": result["contraction_upper"],
        "full_physical_residual_l2": None if full_summary is None else full_summary["l2_upper"],
        "independent_oversampled_verification": independent_agrees,
        "trial_frozen": False,
    }, indent=2))
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", type=int, default=192)
    parser.add_argument("--r-interval", default="6.7439020359331625,6.7439020359331625")
    parser.add_argument("--fourier-cutoff", type=int, required=True)
    parser.add_argument("--candidate-points", type=int, default=20)
    parser.add_argument("--theorem-defect-definition", type=Path, default=ROOT / "track_b_theorem_defect_definition.json")
    parser.add_argument("--row-ledger", type=Path, required=True)
    parser.add_argument("--enrichment-ledger", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    ns = parser.parse_args()
    result = run_cutoff(ns)
    return 0 if result["cutoff_rank_certified"] and result["independent_oversampled_verification"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
