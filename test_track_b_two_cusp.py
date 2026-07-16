from __future__ import annotations

import json
from pathlib import Path
import unittest

from flint import acb, arb, ctx

from track_b_two_cusp_data import (
    GENERATORS, arb_fraction, cusp_modes, exact_collocation_points, h3_action,
    sigma_zero_specialized,
)
from track_b_two_cusp_hejhal import (
    ValidatedWhittaker, assemble_physical_rows, certification_decision,
    component_row, matrix_vector, scale_system, spectral_dependence_test,
)
from track_b_two_cusp_verify import slow_row, verify_from_paths


class TrackBTwoCuspTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        ctx.prec = 192

    def test_sigma0_direct_and_specialized_actions_agree(self) -> None:
        for record in exact_collocation_points(6):
            point = tuple(arb_fraction(record[name]) for name in ("x1", "x2", "y"))
            direct = h3_action(GENERATORS["S"], *point)
            specialized = sigma_zero_specialized(*point)
            for actual, expected in zip(direct, specialized):
                self.assertTrue((actual - expected).contains(0))
                self.assertLess(abs(actual - expected).upper(), arb("1e-40"))

    def test_all_four_blocks_and_slow_direct_rows_agree(self) -> None:
        backend = ValidatedWhittaker(arb("6.7439020359331625"))
        rows, ledger, modes_inf, modes_zero = assemble_physical_rows(
            1, 1, backend, None
        )
        self.assertEqual(
            {record["matrix_block"] for record in ledger},
            {"infinity->infinity", "infinity->zero", "zero->infinity", "zero->zero"},
        )
        for actual, record in zip(rows, ledger):
            expected, transformed = slow_row(
                record, modes_inf, modes_zero, arb("6.7439020359331625")
            )
            self.assertTrue(all(
                (a - b).real.contains(0) and (a - b).imag.contains(0)
                for a, b in zip(actual, expected)
            ))
            self.assertTrue(
                arb(record["true_transformed_height"]).contains(transformed[2])
            )

    def test_spectral_parameter_dependence(self) -> None:
        result = spectral_dependence_test(1, 192)
        self.assertTrue(result["certified"])
        self.assertTrue(result["intervals_disjoint"])

    def test_scaling_and_physical_residual_round_trip(self) -> None:
        matrix = [[acb(2), acb(1)], [acb(1), acb(3)]]
        vector = [acb(4), acb(-2)]
        rhs = [acb(1), acb(5)]
        dr = [arb(2), arb(4)]
        dc = [arb(8), arb(16)]
        scaled_matrix, scaled_rhs = scale_system(matrix, rhs, dr, dc)
        scaled_vector = [vector[j] / dc[j] for j in range(2)]
        physical_residual = [a - b for a, b in zip(matrix_vector(matrix, vector), rhs)]
        scaled_residual = [a - b for a, b in zip(
            matrix_vector(scaled_matrix, scaled_vector), scaled_rhs
        )]
        for i in range(2):
            difference = scaled_residual[i] - dr[i] * physical_residual[i]
            self.assertTrue(difference.real.contains(0))
            self.assertTrue(difference.imag.contains(0))
            roundtrip = dc[i] * scaled_vector[i] - vector[i]
            self.assertTrue(roundtrip.real.contains(0))
            self.assertTrue(roundtrip.imag.contains(0))

    def test_identity_reduction_commutes_with_matrix(self) -> None:
        matrix = [[acb(1), acb(2)], [acb(3), acb(5)]]
        reduced = [acb(7), acb(-4)]
        physical = list(reduced)  # R=I exactly.
        left = matrix_vector(matrix, physical)
        right = matrix_vector(matrix, reduced)
        self.assertTrue(all(
            (a - b).real.contains(0) and (a - b).imag.contains(0)
            for a, b in zip(left, right)
        ))

    def test_each_certification_gate_fails_closed(self) -> None:
        conditions = {
            "all_four_blocks": True,
            "cusp_data_exact": True,
            "true_heights": True,
            "sigma0_consistency": True,
            "spectral_dependence": True,
            "bessel_fallback_zero": True,
            "ledger_complete": True,
            "hashes_deterministic": True,
            "verified_solve": True,
            "normalization": True,
            "physical_reconstruction": True,
            "independent_verification": True,
            "finite_residual": True,
            "no_regularization": True,
        }
        self.assertTrue(certification_decision(conditions)["physical_residual_certified"])
        for name in conditions:
            broken = dict(conditions)
            broken[name] = False
            decision = certification_decision(broken)
            self.assertFalse(decision["physical_residual_certified"], name)
            self.assertFalse(decision["rung4_certified"], name)
            if name in {
                "all_four_blocks", "cusp_data_exact", "true_heights",
                "sigma0_consistency", "spectral_dependence",
                "bessel_fallback_zero", "ledger_complete", "hashes_deterministic",
            }:
                self.assertFalse(decision["two_cusp_assembly_certified"], name)

    def test_forced_bessel_failure_is_fail_closed(self) -> None:
        backend = ValidatedWhittaker(arb(6))
        with self.assertRaises(ArithmeticError):
            backend.kir(arb(-1))
        self.assertEqual(backend.failure_count, 1)
        conditions = {
            "all_four_blocks": True, "cusp_data_exact": True,
            "true_heights": True, "sigma0_consistency": True,
            "spectral_dependence": True, "bessel_fallback_zero": False,
            "ledger_complete": True, "hashes_deterministic": True,
        }
        decision = certification_decision(conditions)
        self.assertFalse(decision["two_cusp_assembly_certified"])
        self.assertFalse(decision["physical_residual_certified"])
        self.assertFalse(decision["rung4_certified"])

    def test_closing_artifacts_verify_independently(self) -> None:
        result_path = Path("track_b_two_cusp_result.json")
        ledger_path = Path("track_b_hejhal_rows.jsonl")
        if not result_path.is_file() or not ledger_path.is_file():
            self.skipTest("closing artifacts have not been generated")
        result = json.loads(result_path.read_text(encoding="utf-8"))
        self.assertTrue(result["two_cusp_assembly_certified"])
        self.assertTrue(result["verified_solve_certified"])
        self.assertTrue(result["physical_residual_certified"])
        self.assertFalse(result["rung4_certified"])
        self.assertTrue(verify_from_paths(result_path, ledger_path)["verified"])


if __name__ == "__main__":
    unittest.main()

