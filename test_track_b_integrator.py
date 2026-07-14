from __future__ import annotations

import json
from pathlib import Path
import unittest

from track_b_rung4_integrator import integrate


class TrackBIntegratorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.mass = json.loads(
            Path("track_b_projected_mass_arb_result.json").read_text(encoding="utf-8")
        )
        cls.partition = {
            "certified": True,
            "coverage_certified": True,
            "local_finiteness_certified": True,
            "transitions_complete": True,
            "stabilizers_certified": True,
            "theorem_DK_compatible": True,
            "active_transition_ids": ["g0"],
            "b0_upper": "1e0",
            "b1_upper": "1e0",
        }
        cls.overlap = {
            "certified": True,
            "all_transitions_covered": True,
            "common_fiber_transport_certified": True,
            "stabilizer_averaging_certified": True,
            "two_cusp_coordinates_certified": True,
            "first_gradients_certified": True,
            "theorem_DK_compatible": True,
            "active_transition_ids": ["g0"],
            "tau_upper": "0",
            "delta0_upper": "1e-3",
            "delta1_upper": "1e-3",
        }
        cls.paths = {"mass": "m", "partition": "p", "overlap": "o"}

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


if __name__ == "__main__":
    unittest.main()
