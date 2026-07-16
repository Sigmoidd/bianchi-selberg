from __future__ import annotations

import json
from pathlib import Path
import unittest

from track_b_global_partition_arb import apply_stability, certify_partition
from track_b_partition_geometry import build
from track_b_rung4_integrator import integrate


class TrackBIntegratorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.mass = json.loads(
            Path("track_b_projected_mass_arb_result.json").read_text(encoding="utf-8")
        )
        cls.partition_audit = Path("track_b_integrator_partition_test.jsonl")
        geometry = build(6)
        closing = certify_partition(
            geometry, (1, 1, 1), 8, 192, cls.partition_audit, 1, True
        )
        check = certify_partition(
            geometry, (1, 1, 1), 8, 192, None, 1, False
        )
        cls.partition = apply_stability(closing, check)
        cls.overlap = {
            "certified": True,
            "all_transitions_covered": True,
            "common_fiber_transport_certified": True,
            "stabilizer_averaging_certified": True,
            "two_cusp_coordinates_certified": True,
            "first_gradients_certified": True,
            "theorem_DK_compatible": True,
            "active_transition_ids": list(cls.partition["active_transition_ids"]),
            "tau_upper": "0",
            "delta0_upper": "1e-3",
            "delta1_upper": "1e-3",
            "weighted_residual_upper": "1e-6",
        }
        cls.paths = {"mass": "m", "partition": "p", "overlap": "o"}

    @classmethod
    def tearDownClass(cls) -> None:
        cls.partition_audit.unlink(missing_ok=True)

    def test_certified_interval_but_counting_fails_closed(self) -> None:
        out = integrate(
            self.mass,
            self.partition,
            self.overlap,
            "0.1",
            False,
            self.paths,
        )
        self.assertTrue(out["track_b_interval_certified"])
        self.assertTrue(out["hard"]["width_lt_tol"])
        self.assertFalse(out["rung4_certified"])
        self.assertIn("counting_certified", out["blockers"])

    def test_transition_mismatch_fails_closed(self) -> None:
        bad = dict(self.overlap)
        bad["active_transition_ids"] = ["different"]
        out = integrate(
            self.mass, self.partition, bad, "0.1", False, self.paths
        )
        self.assertFalse(out["track_b_interval_certified"])
        self.assertFalse(out["rung4_certified"])
        self.assertIn("transition-set mismatch", out["reason"])

    def valid_floor(self) -> dict:
        return {
            "floor_residual_certified": True,
            "continuum_remainder_certified": True,
            "stability_check_passed": True,
            "geometry_incidence_certified": True,
            "projected_symmetries_certified": True,
            "bessel_fallback_count": 0,
            "floor_l2_upper": "0.005",
            "allowed_budget_lower": "0.010",
            "floor_width": self.partition["floor_width"],
            "floor_weight_formula_id": self.partition[
                "floor_weight_formula_id"
            ],
            "floor_geometry_incidence_hash": self.partition[
                "floor_geometry_incidence_hash"
            ],
        }

    def test_floor_endpoint_comparison_is_consumed_as_interval(self) -> None:
        out = integrate(
            self.mass, self.partition, self.overlap, "0.1", False,
            self.paths, self.valid_floor(),
        )
        self.assertTrue(out["floor_comparison"]["endpoint_comparison_certified"])
        self.assertTrue(out["hard"]["floor_residual_certified"])
        self.assertFalse(out["rung4_certified"])

    def test_each_floor_dependency_fails_closed(self) -> None:
        for key in (
            "floor_residual_certified",
            "continuum_remainder_certified",
            "stability_check_passed",
            "geometry_incidence_certified",
            "projected_symmetries_certified",
        ):
            with self.subTest(key=key):
                floor = self.valid_floor()
                floor[key] = False
                out = integrate(
                    self.mass, self.partition, self.overlap, "0.1", False,
                    self.paths, floor,
                )
                self.assertFalse(out["rung4_certified"])
                self.assertIn("floor certificate flags", out["reason"])

    def test_floor_fallback_and_budget_overlap_fail_closed(self) -> None:
        fallback = self.valid_floor()
        fallback["bessel_fallback_count"] = 1
        out = integrate(
            self.mass, self.partition, self.overlap, "0.1", False,
            self.paths, fallback,
        )
        self.assertFalse(out["rung4_certified"])
        self.assertIn("Bessel fallback", out["reason"])

        overlap = self.valid_floor()
        overlap["floor_l2_upper"] = "0.010"
        out = integrate(
            self.mass, self.partition, self.overlap, "0.1", False,
            self.paths, overlap,
        )
        self.assertFalse(out["rung4_certified"])
        self.assertIn("endpoint comparison", out["reason"])

    def test_missing_global_weight_proof_fails_closed(self) -> None:
        partition = dict(self.partition)
        partition["partition_constants_certified"] = False
        out = integrate(
            self.mass, partition, self.overlap, "0.1", False,
            self.paths, self.valid_floor(),
        )
        self.assertFalse(out["rung4_certified"])
        self.assertIn("partition flags", out["reason"])

    def test_partition_audit_and_compatibility_fail_closed(self) -> None:
        for key, value, expected in (
            ("provisional", True, "provisional"),
            ("partition_definition_hash", "0" * 64, "definition hash"),
            ("audit_sha256", "0" * 64, "audit ledger hash"),
        ):
            with self.subTest(key=key):
                partition = dict(self.partition)
                partition[key] = value
                out = integrate(
                    self.mass, partition, self.overlap, "0.1", False,
                    self.paths, self.valid_floor(),
                )
                self.assertFalse(out["rung4_certified"])
                self.assertIn(expected, out["reason"])

        floor = self.valid_floor()
        floor["floor_weight_formula_id"] = "different"
        out = integrate(
            self.mass, self.partition, self.overlap, "0.1", False,
            self.paths, floor,
        )
        self.assertIn("formula mismatch", out["reason"])

        floor = self.valid_floor()
        floor["floor_width"] = "0.25"
        out = integrate(
            self.mass, self.partition, self.overlap, "0.1", False,
            self.paths, floor,
        )
        self.assertIn("width mismatch", out["reason"])


if __name__ == "__main__":
    unittest.main()
