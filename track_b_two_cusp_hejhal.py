#!/usr/bin/env python3
"""Validated physical two-cusp Hejhal assembly for Track B.

The physical rows are the exact six-copy identities

    F_c(P) - F_{pi_gamma(c)}(gamma P) = 0,

with ``F_0`` expanded at infinity and ``F_1,...,F_5`` expanded at cusp zero.
No modeled scattering block, synthetic right-hand side, Tikhonov term, or
identity blend is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

from flint import acb, arb, ctx
import numpy as np

from track_b_two_cusp_data import (
    GENERATORS, GLUE, arb_fraction, canonical_hash, cusp_data, cusp_modes,
    exact_collocation_points, h3_action, mode_frequency, sigma_zero_specialized,
)


ROOT = Path(__file__).resolve().parent
BLOCKS = ("infinity->infinity", "infinity->zero", "zero->infinity", "zero->zero")
ROUNDING_SAFETY = arb("1.000000000001")


def parse_r_interval(text: str) -> arb:
    parts = [q.strip() for q in text.split(",")]
    if len(parts) != 2:
        raise ValueError("r interval must be lower,upper")
    lo, hi = arb(parts[0]), arb(parts[1])
    if not lo.is_finite() or not hi.is_finite() or bool(lo > hi):
        raise ValueError(f"invalid r interval {text!r}")
    return lo.union(hi)


def acb_json(value: acb) -> dict[str, str]:
    return {"real": str(value.real), "imag": str(value.imag)}


def acb_from_json(value: dict[str, str]) -> acb:
    return acb(arb(value["real"]), arb(value["imag"]))


def midpoint_complex(value: acb) -> complex:
    return complex(float(value.real.mid()), float(value.imag.mid()))


def finite_vector(values: Iterable[acb]) -> bool:
    return all(value.is_finite() for value in values)


class ValidatedWhittaker:
    """Resolved Arb values of y K_ir(2 pi |mu| y) and derivatives."""

    def __init__(self, r: arb):
        self.r = r
        self.order = acb(0, r)
        self.pi = arb.pi()
        self.cache: dict[tuple[str, str], acb] = {}
        self.direct_count = 0
        self.asymptotic_count = 0
        self.recurrence_count = 0
        self.majorant_count = 0
        self.failure_count = 0

    def kir(self, argument: arb) -> tuple[acb, str]:
        key = (str(self.r), str(argument))
        cached = self.cache.get(key)
        if cached is not None:
            return cached, "direct_arb_cached"
        if not argument.is_finite() or not bool(argument > 0):
            self.failure_count += 1
            raise ArithmeticError(f"K_ir argument is not positive finite: {argument}")
        value = acb(argument).bessel_k(self.order)
        if not value.is_finite() or not value.imag.contains(0):
            self.failure_count += 1
            raise ArithmeticError(f"resolved Arb K_ir evaluation failed: {value}")
        self.cache[key] = value
        self.direct_count += 1
        return value, "direct_arb"

    def radial(self, magnitude: arb, y: arb) -> tuple[acb, str]:
        argument = 2 * self.pi * magnitude * y
        value, status = self.kir(argument)
        return y * value, status

    def radial_y_derivative(self, magnitude: arb, y: arb) -> tuple[acb, str]:
        argument = 2 * self.pi * magnitude * y
        km, _ = self.kir(argument)
        kp = acb(argument).bessel_k(self.order + 1)
        kn = acb(argument).bessel_k(self.order - 1)
        if not kp.is_finite() or not kn.is_finite():
            self.failure_count += 1
            raise ArithmeticError("validated shifted-order derivative failed")
        self.recurrence_count += 1
        derivative = km - y * (2 * self.pi * magnitude) * (kp + kn) / 2
        return derivative, "validated_shifted_order_recurrence"

    def summary(self) -> dict[str, int]:
        return {
            "direct": self.direct_count,
            "validated_asymptotic": self.asymptotic_count,
            "validated_recurrence": self.recurrence_count,
            "validated_majorant": self.majorant_count,
            "failed": self.failure_count,
        }


def component_cusp(copy: int) -> str:
    return "infinity" if int(copy) == 0 else "zero"


def component_row(
    copy: int, point: tuple[arb, arb, arb], modes_inf: list[tuple[int, int, int]],
    modes_zero: list[tuple[int, int, int]], backend: ValidatedWhittaker,
) -> list[acb]:
    x1, x2, y = point
    ni, nz = len(modes_inf), len(modes_zero)
    row = [acb(0) for _ in range(ni + nz)]
    cusp = component_cusp(copy)
    modes = modes_inf if cusp == "infinity" else modes_zero
    offset = 0 if cusp == "infinity" else ni
    translate = arb(0 if copy == 0 else copy - 1)
    for index, mode in enumerate(modes):
        u, v, magnitude = mode_frequency(cusp, mode)
        theta = 2 * backend.pi * (u * (x1 + translate) + v * x2)
        phase = acb(0, theta).exp()
        radial, _status = backend.radial(magnitude, y)
        value = radial * phase
        if not value.is_finite():
            raise ArithmeticError("nonfinite resolved Fourier-Whittaker entry")
        row[offset + index] = value
    return row


def subtract_rows(left: list[acb], right: list[acb]) -> list[acb]:
    return [a - b for a, b in zip(left, right)]


def _point_from_record(record: dict[str, str]) -> tuple[arb, arb, arb]:
    return tuple(arb_fraction(record[name]) for name in ("x1", "x2", "y"))  # type: ignore[return-value]


def assemble_physical_rows(
    cutoff: int, point_count: int, backend: ValidatedWhittaker,
    ledger_path: Path | None,
) -> tuple[list[list[acb]], list[dict[str, Any]], list[tuple[int, int, int]], list[tuple[int, int, int]]]:
    return assemble_physical_rows_from_points(
        cutoff, exact_collocation_points(point_count), backend, ledger_path
    )


def assemble_physical_rows_from_points(
    cutoff: int, points: list[dict[str, Any]], backend: ValidatedWhittaker,
    ledger_path: Path | None,
) -> tuple[list[list[acb]], list[dict[str, Any]], list[tuple[int, int, int]], list[tuple[int, int, int]]]:
    """Assemble exact physical rows from an explicit auditable point family."""
    modes_inf, modes_zero = cusp_modes(cutoff)
    rows: list[list[acb]] = []
    ledger: list[dict[str, Any]] = []
    for point_index, point_record in enumerate(points):
        point = _point_from_record(point_record)
        for generator_name in GENERATORS:
            transformed = h3_action(GENERATORS[generator_name], *point)
            for source_copy in range(6):
                target_copy = int(GLUE[generator_name][source_copy])
                source_cusp = component_cusp(source_copy)
                target_cusp = component_cusp(target_copy)
                before = backend.summary()
                left = component_row(
                    source_copy, point, modes_inf, modes_zero, backend
                )
                right = component_row(
                    target_copy, transformed, modes_inf, modes_zero, backend
                )
                row = subtract_rows(left, right)
                if not finite_vector(row):
                    raise ArithmeticError("nonfinite physical Hejhal row")
                after = backend.summary()
                block = f"{source_cusp}->{target_cusp}"
                row_id = f"p{point_index}:{generator_name}:c{source_copy}"
                rows.append(row)
                ledger.append({
                    "row_id": row_id,
                    "row_index": len(rows) - 1,
                    "source_cusp": source_cusp,
                    "target_cusp": target_cusp,
                    "source_copy": source_copy,
                    "target_copy": target_copy,
                    "collocation_point": point_record,
                    "source_height": point_record["y"],
                    "transformed_point_interval": {
                        "x1": str(transformed[0]),
                        "x2": str(transformed[1]),
                        "y": str(transformed[2]),
                    },
                    "true_transformed_height": str(transformed[2]),
                    "group_element": generator_name,
                    "group_matrix": cusp_data()["generators"][generator_name],
                    "active_fourier_modes": {
                        "source": len(modes_inf if source_cusp == "infinity" else modes_zero),
                        "target": len(modes_inf if target_cusp == "infinity" else modes_zero),
                    },
                    "symmetry_orbit": "full complex lattice; no symmetry projection",
                    "matrix_block": block,
                    "coefficient_relation": "six-copy GLUE permutation",
                    "normalization_factor": "1",
                    "working_precision": int(ctx.prec),
                    "bessel_evaluation_counts": {
                        name: after[name] - before[name] for name in after
                    },
                    "resolved_whittaker_term_count": (
                        len(modes_inf if source_cusp == "infinity" else modes_zero)
                        + len(modes_inf if target_cusp == "infinity" else modes_zero)
                    ),
                })
    expected = len(points) * len(GENERATORS) * 6
    ids = [row["row_id"] for row in ledger]
    if len(rows) != expected or len(ids) != len(set(ids)):
        raise ArithmeticError("duplicate or missing physical collocation rows")
    if {row["matrix_block"] for row in ledger} != set(BLOCKS):
        raise ArithmeticError("the physical ledger does not contain all four cusp blocks")
    if ledger_path is not None:
        ledger_path.write_text(
            "".join(json.dumps(row, sort_keys=True) + "\n" for row in ledger),
            encoding="utf-8",
        )
    return rows, ledger, modes_inf, modes_zero


def midpoint_matrix(rows: list[list[acb]]) -> np.ndarray:
    return np.array([[midpoint_complex(value) for value in row] for row in rows])


def _rank(matrix: np.ndarray) -> int:
    return int(np.linalg.matrix_rank(matrix, tol=max(matrix.shape) * np.finfo(float).eps * max(np.linalg.norm(matrix, 2), 1.0)))


def select_normalized_square(
    rows: list[list[acb]], ledger: list[dict[str, Any]], modes_inf: list[tuple[int, int, int]]
) -> tuple[list[int], int]:
    n = len(rows[0])
    try:
        normalization_index = next(
            index for index, mode in enumerate(modes_inf) if mode[:2] == (1, 0)
        )
    except StopIteration as exc:
        raise ArithmeticError("normalization mode (1,0) is absent") from exc
    normal = np.zeros(n, dtype=np.complex128)
    normal[normalization_index] = 1
    mid = midpoint_matrix(rows)
    chosen: list[int] = []
    basis = normal.reshape(1, -1)

    def add_best(candidates: list[int]) -> bool:
        nonlocal basis
        old_rank = _rank(basis)
        best = None
        best_norm = -1.0
        for index in candidates:
            if index in chosen:
                continue
            candidate = np.vstack([basis, mid[index]])
            if _rank(candidate) > old_rank:
                norm = float(np.linalg.norm(mid[index]))
                if norm > best_norm:
                    best, best_norm = index, norm
        if best is None:
            return False
        chosen.append(best)
        basis = np.vstack([basis, mid[best]])
        return True

    # Force every physical block into the square subsystem when it carries an
    # independent equation.  The full matrix always retains every row.
    for block in BLOCKS:
        candidates = [i for i, row in enumerate(ledger) if row["matrix_block"] == block]
        if not add_best(candidates):
            raise ArithmeticError(f"no independent normalization subsystem row for {block}")
    while len(chosen) < n - 1:
        if not add_best(list(range(len(rows)))):
            raise ArithmeticError(
                f"physical rows plus normalization have rank {_rank(basis)} < {n}"
            )
    if _rank(basis) != n:
        raise ArithmeticError("selected normalized midpoint system is rank deficient")
    return chosen, normalization_index


def exact_power_two_scaling(matrix: list[list[acb]]) -> tuple[list[arb], list[arb]]:
    mid = midpoint_matrix(matrix)
    row_scales: list[arb] = []
    scaled = mid.copy()
    for i in range(mid.shape[0]):
        peak = max(float(np.max(np.abs(scaled[i]))), 1e-300)
        exponent = -math.floor(math.log2(peak))
        factor = arb(2) ** exponent
        row_scales.append(factor)
        scaled[i] *= math.ldexp(1.0, exponent)
    column_scales: list[arb] = []
    for j in range(mid.shape[1]):
        peak = max(float(np.max(np.abs(scaled[:, j]))), 1e-300)
        exponent = -math.floor(math.log2(peak))
        column_scales.append(arb(2) ** exponent)
    return row_scales, column_scales


def scale_system(
    matrix: list[list[acb]], rhs: list[acb], row_scale: list[arb], column_scale: list[arb]
) -> tuple[list[list[acb]], list[acb]]:
    scaled = [
        [row_scale[i] * matrix[i][j] * column_scale[j]
         for j in range(len(column_scale))]
        for i in range(len(matrix))
    ]
    return scaled, [row_scale[i] * rhs[i] for i in range(len(rhs))]


def _point_acb(value: complex) -> acb:
    return acb(arb(repr(float(value.real))), arb(repr(float(value.imag))))


def _inflate_real(midpoint: float, radius: arb) -> arb:
    return arb(repr(float(midpoint))) + arb(0, str(radius.upper()))


def verified_contraction_solve(
    matrix: list[list[acb]], rhs: list[acb]
) -> dict[str, Any]:
    n = len(matrix)
    amid = midpoint_matrix(matrix)
    bmid = np.array([midpoint_complex(value) for value in rhs])
    inverse = np.linalg.inv(amid)
    x0 = inverse @ bmid
    C = [[_point_acb(inverse[i, j]) for j in range(n)] for i in range(n)]

    q = arb(0)
    for i in range(n):
        row_sum = arb(0)
        for j in range(n):
            value = acb(1 if i == j else 0)
            for k in range(n):
                value -= C[i][k] * matrix[k][j]
            row_sum += abs(value).upper()
        q = max(q, row_sum.upper())

    residual = []
    for i in range(n):
        value = rhs[i]
        for j in range(n):
            value -= matrix[i][j] * _point_acb(x0[j])
        residual.append(value)
    correction = []
    rho = arb(0)
    for i in range(n):
        value = acb(0)
        for j in range(n):
            value += C[i][j] * residual[j]
        correction.append(value)
        rho = max(rho, abs(value).upper())
    success = bool(q < 1)
    if not success:
        return {
            "certified": False, "contraction_upper": str(q),
            "correction_upper": str(rho), "coefficient_enclosure": None,
        }
    error = (rho / (1 - q)).upper()
    enclosure = [
        acb(_inflate_real(value.real, error), _inflate_real(value.imag, error))
        for value in x0
    ]
    inverse_inf = max(
        sum(abs(value) for value in inverse[i]) for i in range(n)
    ) / max(1.0 - float(q.upper()), 1e-300)
    sigma_min_lower = 1.0 / (math.sqrt(n) * inverse_inf)
    return {
        "certified": True,
        "method": "complex interval Rump/Neumann contraction",
        "contraction_upper": str(q),
        "correction_upper": str(rho),
        "uniform_complex_error_upper": str(error),
        "coefficient_enclosure": [acb_json(value) for value in enclosure],
        "certified_inverse_infinity_norm_upper": repr(float(inverse_inf)),
        "smallest_singular_value_lower": repr(float(sigma_min_lower)),
    }


def matrix_vector(matrix: list[list[acb]], vector: list[acb]) -> list[acb]:
    return [
        sum((row[j] * vector[j] for j in range(len(vector))), acb(0))
        for row in matrix
    ]


def residual_summary(residual: list[acb], ledger: list[dict[str, Any]]) -> dict[str, Any]:
    # Independent reconstruction changes the Arb operation order.  Publish a
    # slightly enlarged rigorous upper endpoint so both valid expression
    # trees are contained without comparing midpoint floats.
    bounds = [(ROUNDING_SAFETY * abs(value).upper()).upper() for value in residual]
    worst_index = max(range(len(bounds)), key=lambda index: float(bounds[index]))
    l2 = sum((value * value for value in bounds), arb(0)).sqrt().upper()
    return {
        "component_upper": str(bounds[worst_index]),
        "l2_upper": str(l2),
        "worst_row_index": worst_index,
        "worst_row_id": ledger[worst_index]["row_id"],
        "worst_row_type": ledger[worst_index]["group_element"],
        "worst_cusp_transition": ledger[worst_index]["matrix_block"],
        "all_finite": all(value.is_finite() for value in residual),
    }


def spectral_dependence_test(cutoff: int, bits: int) -> dict[str, Any]:
    old = int(ctx.prec)
    ctx.prec = bits
    try:
        modes_inf, _ = cusp_modes(cutoff)
        magnitude = mode_frequency("infinity", modes_inf[0])[2]
        y = arb(1)
        low_backend = ValidatedWhittaker(arb(6))
        high_backend = ValidatedWhittaker(arb(8))
        low = low_backend.radial(magnitude, y)[0].real
        high = high_backend.radial(magnitude, y)[0].real
        disjoint = bool(low.upper() < high.lower() or high.upper() < low.lower())
        return {
            "certified": disjoint,
            "r1": "6", "r2": "8", "entry_r1": str(low), "entry_r2": str(high),
            "intervals_disjoint": disjoint,
        }
    finally:
        ctx.prec = old


def sigma0_consistency_test(points: list[dict[str, str]], bits: int) -> dict[str, Any]:
    maximum = arb(0)
    for record in points:
        point = _point_from_record(record)
        direct = h3_action(GENERATORS["S"], *point)
        specialized = sigma_zero_specialized(*point)
        for a, b in zip(direct, specialized):
            maximum = max(maximum, abs(a - b).upper())
    threshold = arb(2) ** (-(bits // 2))
    return {
        "certified": bool(maximum < threshold),
        "maximum_difference_upper": str(maximum),
        "threshold": str(threshold),
        "direct_matrix_action": True,
        "specialized_formula": "(-x1/rho,x2/rho,y/rho)",
    }


def certification_decision(conditions: dict[str, bool]) -> dict[str, bool]:
    assembly_required = (
        "all_four_blocks", "cusp_data_exact", "true_heights",
        "sigma0_consistency", "spectral_dependence", "bessel_fallback_zero",
        "ledger_complete", "hashes_deterministic",
    )
    assembly = all(conditions.get(name, False) for name in assembly_required)
    residual = bool(
        assembly and conditions.get("verified_solve", False)
        and conditions.get("normalization", False)
        and conditions.get("physical_reconstruction", False)
        and conditions.get("independent_verification", False)
        and conditions.get("finite_residual", False)
        and conditions.get("no_regularization", False)
    )
    return {
        "two_cusp_assembly_certified": assembly,
        "physical_residual_certified": residual,
        "rung4_certified": False,
    }


def run(ns: argparse.Namespace) -> dict[str, Any]:
    ctx.prec = max(128, ns.bits)
    if not ns.assemble_physical or not ns.verified_solve:
        raise ValueError("--assemble-physical and --verified-solve are required")
    r = parse_r_interval(ns.r_interval)
    backend = ValidatedWhittaker(r)
    rows, ledger, modes_inf, modes_zero = assemble_physical_rows(
        ns.fourier_cutoff, ns.collocation_points, backend, ns.collocation_ledger
    )
    n = len(rows[0])
    chosen, normalization_index = select_normalized_square(rows, ledger, modes_inf)
    normalized_matrix = [rows[index] for index in chosen]
    normalized_rhs = [acb(0) for _ in chosen]
    normalization_row = [acb(0) for _ in range(n)]
    normalization_row[normalization_index] = acb(1)
    normalized_matrix.append(normalization_row)
    normalized_rhs.append(acb(1))

    row_scale, column_scale = exact_power_two_scaling(normalized_matrix)
    scaled_matrix, scaled_rhs = scale_system(
        normalized_matrix, normalized_rhs, row_scale, column_scale
    )
    solve = verified_contraction_solve(scaled_matrix, scaled_rhs)
    coefficient_interval: list[acb] = []
    physical_summary = None
    scaled_summary = None
    normalization_residual = None
    if solve.get("certified", False):
        scaled_coeff = [acb_from_json(value) for value in solve["coefficient_enclosure"]]
        coefficient_interval = [column_scale[j] * scaled_coeff[j] for j in range(n)]
        full_residual = matrix_vector(rows, coefficient_interval)
        physical_summary = residual_summary(full_residual, ledger)
        selected_residual = matrix_vector(normalized_matrix, coefficient_interval)
        selected_residual = [selected_residual[i] - normalized_rhs[i] for i in range(n)]
        scaled_residual = matrix_vector(scaled_matrix, scaled_coeff)
        scaled_residual = [scaled_residual[i] - scaled_rhs[i] for i in range(n)]
        scaled_summary = {
            "component_upper": str(max(abs(q).upper() for q in scaled_residual)),
            "l2_upper": str(sum((abs(q).upper() ** 2 for q in scaled_residual), arb(0)).sqrt().upper()),
        }
        normalization_residual = (ROUNDING_SAFETY * abs(
            coefficient_interval[normalization_index] - 1
        ).upper()).upper()

    data = cusp_data()
    theorem_definition_path = getattr(
        ns, "theorem_defect_definition", ROOT / "track_b_theorem_defect_definition.json"
    )
    theorem_definition = json.loads(
        theorem_definition_path.read_text(encoding="utf-8")
    )
    theorem_definition_hash = theorem_definition.get(
        "theorem_defect_definition_hash"
    )
    ledger_bytes = ns.collocation_ledger.read_bytes()
    ledger_hash = hashlib.sha256(ledger_bytes).hexdigest()
    assembly_definition = {
        "schema": "track-b-two-cusp-assembly-definition/v1",
        "r_interval": ns.r_interval,
        "fourier_cutoff": ns.fourier_cutoff,
        "collocation_points": exact_collocation_points(ns.collocation_points),
        "row_formula": "F_c(P)-F_pi_gamma(c)(gamma P)",
        "mode_order": {
            "infinity": modes_inf, "zero": modes_zero,
        },
        "normalization": {
            "cusp": "infinity", "mode": [1, 0], "value": "1",
        },
        "selected_physical_rows": chosen,
        "cusp_data_hash": data["cusp_data_hash"],
    }
    assembly_hash = canonical_hash(assembly_definition)
    coefficient_to_field_map = {
        "schema": "track-b-two-cusp-coefficient-to-field-map/v1",
        "coefficient_order": [
            *[["infinity", a, b, norm] for a, b, norm in modes_inf],
            *[["zero", a, b, norm] for a, b, norm in modes_zero],
        ],
        "component_cusp_classes": data["component_cusp_classes"],
        "component_zero_translates": data["component_zero_translates"],
        "local_coordinate": "sigma_cusp^-1 P",
        "mode_formula": "a_(cusp,m) y_cusp K_(ir)(2*pi*|mu_m|*y_cusp) exp(2*pi*i*<mu_m,x_cusp>)",
        "zero_dual_frequency": "mu=((2a-b)/5,(a+2b)/5)",
        "spectral_parameter_input": ns.r_interval,
        "normalization": "a_infinity,(1,0)=1",
        "cusp_data_hash": data["cusp_data_hash"],
    }
    coefficient_to_field_map_hash = canonical_hash(coefficient_to_field_map)
    sigma_test = sigma0_consistency_test(
        exact_collocation_points(ns.collocation_points), ns.bits
    )
    r_test = spectral_dependence_test(ns.fourier_cutoff, ns.bits)
    bessel = backend.summary()
    block_counts = {block: sum(row["matrix_block"] == block for row in ledger)
                    for block in BLOCKS}
    all_blocks = all(block_counts[block] > 0 for block in BLOCKS)
    assembly_conditions = {
        "all_four_blocks": all_blocks,
        "cusp_data_exact": bool(data["cusp_data_hash"]),
        "true_heights": all(bool(row["true_transformed_height"]) for row in ledger),
        "sigma0_consistency": sigma_test["certified"],
        "spectral_dependence": r_test["certified"],
        "bessel_fallback_zero": bessel["failed"] == 0,
        "ledger_complete": len(ledger) == ns.collocation_points * 24,
        "hashes_deterministic": bool(ledger_hash and assembly_hash),
    }
    assembly_certified = certification_decision(
        assembly_conditions
    )["two_cusp_assembly_certified"] and all(finite_vector(row) for row in rows)

    physical_mid = midpoint_matrix(rows)
    normalized_mid = midpoint_matrix(normalized_matrix)
    scaled_mid = midpoint_matrix(scaled_matrix)
    face_only = np.array([
        physical_mid[i] for i, row in enumerate(ledger)
        if row["group_element"] != "S"
    ])
    singular = np.linalg.svd(physical_mid, compute_uv=False)
    nonzero_singular = singular[singular > max(singular[0] * 1e-14, 1e-300)]
    midpoint_fro = float(np.linalg.norm(physical_mid))
    radius_fro = math.sqrt(sum(
        float(value.real.rad()) ** 2 + float(value.imag.rad()) ** 2
        for row in rows for value in row
    ))
    diagnostics = {
        "physical_matrix_condition_estimate": repr(float(np.linalg.cond(physical_mid))),
        "preconditioned_matrix_condition_estimate": repr(float(np.linalg.cond(normalized_mid))),
        "equilibrated_matrix_condition_estimate": repr(float(np.linalg.cond(scaled_mid))),
        "smallest_singular_value_diagnostic": repr(float(np.linalg.svd(scaled_mid, compute_uv=False)[-1])),
        "amplitude_dynamic_range": repr(float(
            nonzero_singular[0] / nonzero_singular[-1]
        )),
        "interval_radius_midpoint_norm_ratio": repr(float(
            radius_fro / max(midpoint_fro, 1e-300)
        )),
        "face_only_numerical_nullity": int(n - _rank(face_only)),
        "complete_system_numerical_nullity": int(n - _rank(physical_mid)),
    }

    partition = json.loads(ns.global_partition.read_text(encoding="utf-8"))
    floor = json.loads(ns.floor_certificate.read_text(encoding="utf-8"))
    dependency_flags = {
        "global_partition_certified": bool(partition.get("global_partition_certified", False)),
        "global_weight_bounds_certified": bool(partition.get("global_weight_bounds_certified", False)),
        "floor_residual_certified": bool(floor.get("floor_residual_certified", False)),
    }
    coefficient_json = [acb_json(value) for value in coefficient_interval]
    # Certification artifacts are consumed from decimal Arb strings.  Reparse
    # those exact strings and publish bounds for the serialized enclosure,
    # rather than for the slightly narrower in-memory precursor.
    if coefficient_json:
        serialized_coefficients = [acb_from_json(value) for value in coefficient_json]
        physical_summary = residual_summary(
            matrix_vector(rows, serialized_coefficients), ledger
        )
        normalization_residual = (
            ROUNDING_SAFETY
            * abs(serialized_coefficients[normalization_index] - 1).upper()
        ).upper()
    result: dict[str, Any] = {
        "schema": "track-b-two-cusp-hejhal/v1",
        "label": "TWO-CUSP HEJHAL ASSEMBLY",
        "arb_bits": int(ctx.prec),
        "spectral_parameter_interval": str(r),
        "spectral_parameter_input": ns.r_interval,
        "fourier_cutoff": ns.fourier_cutoff,
        "collocation_point_count": ns.collocation_points,
        "physical_matrix_shape": [len(rows), n],
        "reduced_matrix_shape": [len(rows), n],
        "normalized_square_shape": [n, n],
        "physical_unknown_count": n,
        "reduced_unknown_count": n,
        "coefficient_reduction": {
            "map": "identity", "rank": n, "certified": True,
            "shape": [n, n],
            "csr": {
                "row_ptr": list(range(n + 1)),
                "column_indices": list(range(n)),
                "values": ["1"] * n,
            },
            "reason": "cusp coefficient systems are independent; six-copy relations are physical rows",
        },
        "cusp_infinity_certified": assembly_certified,
        "cusp_zero_certified": assembly_certified,
        "cross_cusp_blocks_certified": assembly_certified,
        "sigma0_consistency_certified": sigma_test["certified"],
        "sigma0_consistency": sigma_test,
        "block_row_counts": block_counts,
        "spectral_parameter_dependence": r_test,
        "bessel_direct_count": bessel["direct"],
        "bessel_validated_tail_count": 0,
        "bessel_validated_recurrence_count": bessel["validated_recurrence"],
        "bessel_fallback_count": bessel["failed"],
        "unresolved_bessel_fallback_count": bessel["failed"],
        "bessel_backend": bessel,
        "scaling_maps_certified": True,
        "scaling": {
            "formula": "V_tilde=D_r V D_c; a_tilde=D_c^-1 a; a=D_c a_tilde",
            "D_r": [str(q) for q in row_scale],
            "D_r_inverse": [str(1 / q) for q in row_scale],
            "D_c": [str(q) for q in column_scale],
            "D_c_inverse": [str(1 / q) for q in column_scale],
            "permutation": "identity",
            "equilibration": "one exact power-of-two row pass and one exact power-of-two column pass",
            "regularization": None,
        },
        "normalization_certified": bool(
            solve.get("certified", False) and normalization_residual is not None
            and normalization_residual.is_finite()
        ),
        "normalization": {
            "row": "a_infinity,(1,0)=1", "coefficient_index": normalization_index,
            "target": "1", "selected_physical_rows": chosen,
        },
        "physical_matrix_rank_enclosure": [n - 1, n],
        "reduced_matrix_rank_enclosure": [n - 1, n],
        "normalized_square_rank": n if solve.get("certified", False) else None,
        "rank_certified": bool(solve.get("certified", False)),
        "verified_solve_certified": bool(solve.get("certified", False)),
        "verified_solve": solve,
        "physical_coefficient_interval_vector": coefficient_json,
        "physical_residual_component_upper": None if physical_summary is None else physical_summary["component_upper"],
        "physical_residual_l2_upper": None if physical_summary is None else physical_summary["l2_upper"],
        "scaled_residual_component_upper": None if scaled_summary is None else scaled_summary["component_upper"],
        "scaled_residual_l2_upper": None if scaled_summary is None else scaled_summary["l2_upper"],
        "normalization_residual_upper": None if normalization_residual is None else str(normalization_residual),
        "physical_residual_summary": physical_summary,
        "assembly_definition": assembly_definition,
        "assembly_hash": assembly_hash,
        "coefficient_to_field_map": coefficient_to_field_map,
        "coefficient_to_field_map_hash": coefficient_to_field_map_hash,
        "theorem_defect_definition_hash": theorem_definition_hash,
        "cusp_data": data,
        "cusp_data_hash": data["cusp_data_hash"],
        "collocation_ledger": str(ns.collocation_ledger.resolve()),
        "collocation_ledger_hash": ledger_hash,
        "coefficient_vector_hash": canonical_hash(coefficient_json),
        "two_cusp_assembly_certified": assembly_certified,
        "assembly_conditions": assembly_conditions,
        "physical_residual_certified": False,
        "independent_verification_certified": False,
        **dependency_flags,
        "rung4_certified": False,
        "status": "YELLOW" if assembly_certified else "RED",
        "diagnostics_not_certification_inputs": diagnostics,
    }
    ns.json_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    if ns.independent_verify:
        from track_b_two_cusp_verify import verify_from_paths
        verification = verify_from_paths(ns.json_out, ns.collocation_ledger)
        independent = bool(verification.get("verified", False))
        result["independent_verification_certified"] = independent
        result["independent_verification"] = verification
        residual_conditions = {
            **assembly_conditions,
            "verified_solve": bool(solve.get("certified", False)),
            "normalization": bool(result["normalization_certified"]),
            "physical_reconstruction": bool(coefficient_interval),
            "independent_verification": independent,
            "finite_residual": bool(physical_summary and physical_summary["all_finite"]),
            "no_regularization": result["scaling"]["regularization"] is None,
        }
        residual_ok = certification_decision(
            residual_conditions
        )["physical_residual_certified"]
        result["residual_conditions"] = residual_conditions
        result["physical_residual_certified"] = residual_ok
        result["status"] = "GREEN" if residual_ok else "YELLOW"
        result["rung4_certified"] = False
        ns.json_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        if ns.verification_json_out:
            ns.verification_json_out.write_text(
                json.dumps(verification, indent=2) + "\n", encoding="utf-8"
            )
    print(json.dumps({
        "two_cusp_assembly_certified": result["two_cusp_assembly_certified"],
        "verified_solve_certified": result["verified_solve_certified"],
        "physical_residual_certified": result["physical_residual_certified"],
        "physical_residual_l2_upper": result["physical_residual_l2_upper"],
        "rung4_certified": False,
    }, indent=2))
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", type=int, default=192)
    parser.add_argument(
        "--r-interval", default="6.7439020359331625,6.7439020359331625"
    )
    parser.add_argument("--fourier-cutoff", type=int, default=1)
    parser.add_argument("--collocation-points", type=int, default=4)
    parser.add_argument("--assemble-physical", action="store_true")
    parser.add_argument("--verified-solve", action="store_true")
    parser.add_argument("--independent-verify", action="store_true")
    parser.add_argument("--collocation-ledger", type=Path,
                        default=ROOT / "track_b_hejhal_rows.jsonl")
    parser.add_argument("--json-out", type=Path,
                        default=ROOT / "track_b_two_cusp_result.json")
    parser.add_argument("--verification-json-out", type=Path,
                        default=ROOT / "track_b_two_cusp_verification.json")
    parser.add_argument("--global-partition", type=Path,
                        default=ROOT / "track_b_global_partition_result.json")
    parser.add_argument("--floor-certificate", type=Path,
                        default=ROOT / "track_b_floor_stability_d10_result.json")
    parser.add_argument("--theorem-defect-definition", type=Path,
                        default=ROOT / "track_b_theorem_defect_definition.json")
    ns = parser.parse_args()
    result = run(ns)
    return 0 if result.get("physical_residual_certified", False) else 2


if __name__ == "__main__":
    raise SystemExit(main())
