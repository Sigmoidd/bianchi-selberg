#!/usr/bin/env python3
"""Independent reconstruction verifier for the Track-B two-cusp system."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from flint import acb, arb, ctx

from track_b_two_cusp_data import (
    GENERATORS, GLUE, arb_fraction, canonical_hash, cusp_data, cusp_modes,
    h3_action, mode_frequency,
)


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()]


def parse_interval(text: str) -> arb:
    lo, hi = [arb(q.strip()) for q in text.split(",")]
    return lo.union(hi)


def decode_acb(value: dict[str, str]) -> acb:
    return acb(arb(value["real"]), arb(value["imag"]))


def slow_component(
    copy: int, point: tuple[arb, arb, arb], modes_inf: list[tuple[int, int, int]],
    modes_zero: list[tuple[int, int, int]], r: arb,
    radial_cache: dict[tuple[str, str], acb] | None = None,
) -> list[acb]:
    x1, x2, y = point
    ni, nz = len(modes_inf), len(modes_zero)
    row = [acb(0) for _ in range(ni + nz)]
    cusp = "infinity" if copy == 0 else "zero"
    modes = modes_inf if cusp == "infinity" else modes_zero
    offset = 0 if cusp == "infinity" else ni
    translate = arb(0 if copy == 0 else copy - 1)
    pi = arb.pi()
    order = acb(0, r)
    for index, mode in enumerate(modes):
        u, v, magnitude = mode_frequency(cusp, mode)
        argument = 2 * pi * magnitude * y
        key = (str(magnitude), str(y))
        bessel = None if radial_cache is None else radial_cache.get(key)
        if bessel is None:
            bessel = acb(argument).bessel_k(order)
            if not bessel.is_finite() or not bessel.imag.contains(0):
                raise ArithmeticError("independent direct Arb K_ir evaluation failed")
            if radial_cache is not None:
                radial_cache[key] = bessel
        phase = acb(0, 2 * pi * (u * (x1 + translate) + v * x2)).exp()
        row[offset + index] = y * bessel * phase
    return row


def slow_row(
    record: dict[str, Any], modes_inf: list[tuple[int, int, int]],
    modes_zero: list[tuple[int, int, int]], r: arb,
    radial_cache: dict[tuple[str, str], acb] | None = None,
) -> tuple[list[acb], tuple[arb, arb, arb]]:
    point_record = record["collocation_point"]
    point = tuple(arb_fraction(point_record[name]) for name in ("x1", "x2", "y"))
    transformed = h3_action(GENERATORS[record["group_element"]], *point)
    left = slow_component(
        record["source_copy"], point, modes_inf, modes_zero, r, radial_cache
    )
    right = slow_component(
        record["target_copy"], transformed, modes_inf, modes_zero, r, radial_cache
    )
    return [a - b for a, b in zip(left, right)], transformed


def matvec(matrix: list[list[acb]], vector: list[acb]) -> list[acb]:
    return [sum((row[j] * vector[j] for j in range(len(vector))), acb(0))
            for row in matrix]


def verify(result: dict[str, Any], ledger: list[dict[str, Any]], ledger_hash: str) -> dict[str, Any]:
    ctx.prec = max(128, int(result.get("arb_bits", 0)))
    cutoff = int(result["fourier_cutoff"])
    r = parse_interval(result["spectral_parameter_input"])
    modes_inf, modes_zero = cusp_modes(cutoff)
    n = len(modes_inf) + len(modes_zero)
    expected = int(result["collocation_point_count"]) * len(GENERATORS) * 6
    ids = [str(row.get("row_id")) for row in ledger]
    exact_once = len(ledger) == expected and len(ids) == len(set(ids))
    index_order = [int(row.get("row_index", -1)) for row in ledger]
    exact_order = index_order == list(range(expected))
    all_blocks = {row.get("matrix_block") for row in ledger} == {
        "infinity->infinity", "infinity->zero", "zero->infinity", "zero->zero"
    }
    glue_exact = all(
        int(row["target_copy"])
        == int(GLUE[row["group_element"]][int(row["source_copy"])])
        for row in ledger
    )

    rows: list[list[acb]] = []
    transform_agreement = True
    radial_cache: dict[tuple[str, str], acb] = {}
    for record in ledger:
        row, transformed = slow_row(
            record, modes_inf, modes_zero, r, radial_cache
        )
        rows.append(row)
        stored = record["transformed_point_interval"]
        transform_agreement = transform_agreement and all(
            arb(stored[name]).contains(value)
            for name, value in zip(("x1", "x2", "y"), transformed)
        )

    coefficients = [decode_acb(value)
                    for value in result.get("physical_coefficient_interval_vector", [])]
    coefficient_count = len(coefficients) == n
    residual = matvec(rows, coefficients) if coefficient_count else []
    residual_bounds = [abs(value).upper() for value in residual]
    residual_finite = bool(residual) and all(value.is_finite() for value in residual)
    worst = max(range(len(residual)), key=lambda i: float(residual_bounds[i])) if residual else None
    residual_l2 = (
        sum((value * value for value in residual_bounds), arb(0)).sqrt().upper()
        if residual else arb("+inf")
    )
    residual_component = residual_bounds[worst] if worst is not None else arb("+inf")

    normalization = result["normalization"]
    norm_index = int(normalization["coefficient_index"])
    norm_residual = coefficients[norm_index] - 1 if coefficient_count else acb("nan")
    normalization_finite = norm_residual.is_finite()

    selected = [int(q) for q in normalization["selected_physical_rows"]]
    normalization_row = [acb(0) for _ in range(n)]
    normalization_row[norm_index] = acb(1)
    selected_matrix = [rows[index] for index in selected] + [normalization_row]
    selected_rhs = [acb(0) for _ in selected] + [acb(1)]
    selected_residual = [a - b for a, b in zip(
        matvec(selected_matrix, coefficients), selected_rhs
    )] if coefficient_count else []
    selected_contains_zero = bool(selected_residual) and all(
        value.real.contains(0) and value.imag.contains(0) for value in selected_residual
    )

    scaling = result["scaling"]
    dc = [arb(text) for text in scaling["D_c"]]
    scaled_coeff = [decode_acb(value) for value in
                    result.get("verified_solve", {}).get("coefficient_enclosure", [])]
    scaling_roundtrip = bool(
        len(scaled_coeff) == n and coefficient_count
        and all(
            (coefficients[j] - dc[j] * scaled_coeff[j]).real.contains(0)
            and (coefficients[j] - dc[j] * scaled_coeff[j]).imag.contains(0)
            for j in range(n)
        )
    )

    assembly_hash = canonical_hash(result["assembly_definition"])
    data = cusp_data()
    deterministic_hashes = bool(
        assembly_hash == result.get("assembly_hash")
        and data["cusp_data_hash"] == result.get("cusp_data_hash")
        and canonical_hash(result.get("physical_coefficient_interval_vector", []))
            == result.get("coefficient_vector_hash")
        and ledger_hash == result.get("collocation_ledger_hash")
    )
    endpoint_checks = {
        "component": bool(
            residual_component
            <= arb(result["physical_residual_component_upper"]).upper()
        ),
        "l2": bool(
            residual_l2 <= arb(result["physical_residual_l2_upper"]).upper()
        ),
        "normalization": bool(
            abs(norm_residual).upper()
            <= arb(result["normalization_residual_upper"]).upper()
        ),
    }
    endpoint_agreement = all(endpoint_checks.values())
    verified = bool(
        exact_once and exact_order and all_blocks and glue_exact
        and transform_agreement and coefficient_count and residual_finite
        and normalization_finite and selected_contains_zero and scaling_roundtrip
        and deterministic_hashes and endpoint_agreement
        and result.get("two_cusp_assembly_certified", False)
        and result.get("verified_solve_certified", False)
        and result.get("rung4_certified") is False
    )
    return {
        "schema": "track-b-two-cusp-verification/v1",
        "verified": verified,
        "independently_reassembled_matrix": True,
        "serialized_matrix_entries_used": False,
        "expected_row_count": expected,
        "row_count": len(ledger),
        "exact_once": exact_once,
        "exact_row_order": exact_order,
        "all_four_blocks": all_blocks,
        "glue_permutations_exact": glue_exact,
        "transformed_coordinates_agree": transform_agreement,
        "coefficient_count_correct": coefficient_count,
        "selected_system_contains_zero": selected_contains_zero,
        "scaling_roundtrip_certified": scaling_roundtrip,
        "deterministic_hashes_match": deterministic_hashes,
        "endpoint_agreement": endpoint_agreement,
        "endpoint_checks": endpoint_checks,
        "maximum_component_residual_upper": str(residual_component),
        "euclidean_residual_norm_upper": str(residual_l2),
        "normalization_residual_upper": str(abs(norm_residual).upper()),
        "worst_row_index": worst,
        "worst_row_type": None if worst is None else ledger[worst]["group_element"],
        "worst_cusp_transition": None if worst is None else ledger[worst]["matrix_block"],
        "physical_residual_certified": verified,
        "rung4_certified": False,
    }


def verify_from_paths(result_path: Path, ledger_path: Path) -> dict[str, Any]:
    result = json.loads(result_path.read_text(encoding="utf-8"))
    ledger = load_jsonl(ledger_path)
    ledger_hash = hashlib.sha256(ledger_path.read_bytes()).hexdigest()
    return verify(result, ledger, ledger_hash)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--collocation-ledger", type=Path, required=True)
    parser.add_argument("--json-out", type=Path,
                        default=Path("track_b_two_cusp_verification.json"))
    ns = parser.parse_args()
    output = verify_from_paths(ns.result, ns.collocation_ledger)
    ns.json_out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))
    return 0 if output.get("verified", False) else 2


if __name__ == "__main__":
    raise SystemExit(main())
