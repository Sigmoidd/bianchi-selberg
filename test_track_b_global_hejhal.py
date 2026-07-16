from __future__ import annotations

import copy
import unittest

from track_b_global_hejhal_verify import EXPECTED_CHANNELS, structural_checks


def closing_payload() -> dict:
    channels = [
        {"id": name, "upper": "0.001", "certified": True}
        for name in sorted(EXPECTED_CHANNELS)
    ]
    return {
        "schema": "track-b-global-hejhal-defect/v1",
        "channel_ledger": channels,
        "aggregation_plan": [{
            "rule": "triangle",
            "channels": sorted(EXPECTED_CHANNELS),
        }],
        "global_defect_upper": "[0.007 +/- 1e-12]",
        "allowed_defect_lower": "0.01",
        "gradient_pullback_certified": True,
        "fourier_tails_certified": True,
        "tail_majorant_retained_mode_count": 0,
        "reprojection_hashes_match": True,
        "floor_certificate_compatible": True,
        "physical_collocation_is_theorem_channel": False,
        "global_threshold_certified": True,
        "threshold_trial_compatibility_certified": True,
        "cutoff_stability_passed": True,
        "resolved_bessel_fallback_count": 0,
        "global_hejhal_defect_certified": False,
        "rung4_certified": False,
        "hejhal_existence_certified": False,
        "dual_certification": False,
    }


class TrackBGlobalHejhalFailClosedTests(unittest.TestCase):
    def assertRejected(self, payload: dict) -> None:
        checks = structural_checks(payload)
        self.assertFalse(all(checks.values()))
        self.assertFalse(payload["global_hejhal_defect_certified"])
        self.assertFalse(payload["rung4_certified"])
        self.assertFalse(payload["hejhal_existence_certified"])
        self.assertFalse(payload["dual_certification"])

    def test_synthetic_complete_structure_is_internally_consistent(self) -> None:
        self.assertTrue(all(structural_checks(closing_payload()).values()))

    def test_missing_cusp_zero_continuum_channel(self) -> None:
        p = closing_payload()
        p["channel_ledger"] = [q for q in p["channel_ledger"] if q["id"] != "cusp_zero_continuum"]
        self.assertRejected(p)

    def test_missing_nonfloor_face(self) -> None:
        p = closing_payload()
        p["channel_ledger"] = [q for q in p["channel_ledger"] if q["id"] != "nonfloor_face_value"]
        self.assertRejected(p)

    def test_duplicated_face_contribution(self) -> None:
        p = closing_payload()
        p["channel_ledger"].append(copy.deepcopy(next(q for q in p["channel_ledger"] if q["id"] == "nonfloor_face_value")))
        self.assertRejected(p)

    def test_wrong_gradient_pullback(self) -> None:
        p = closing_payload(); p["gradient_pullback_certified"] = False
        self.assertRejected(p)

    def test_uncertified_fourier_tail(self) -> None:
        p = closing_payload(); p["fourier_tails_certified"] = False
        self.assertRejected(p)

    def test_tail_majorant_applied_to_retained_modes(self) -> None:
        p = closing_payload(); p["tail_majorant_retained_mode_count"] = 1
        self.assertRejected(p)

    def test_reprojection_hash_mismatch(self) -> None:
        p = closing_payload(); p["reprojection_hashes_match"] = False
        self.assertRejected(p)

    def test_floor_certificate_incompatible(self) -> None:
        p = closing_payload(); p["floor_certificate_compatible"] = False
        self.assertRejected(p)

    def test_collocation_substituted_for_continuum(self) -> None:
        p = closing_payload(); p["physical_collocation_is_theorem_channel"] = True
        self.assertRejected(p)

    def test_quadratic_aggregation_without_disjointness(self) -> None:
        p = closing_payload(); p["aggregation_plan"][0]["rule"] = "disjoint_l2"
        p["aggregation_plan"][0]["justification_certified"] = False
        self.assertRejected(p)

    def test_global_threshold_replaced_by_local_floor_budget(self) -> None:
        p = closing_payload(); p["threshold_trial_compatibility_certified"] = False
        self.assertRejected(p)

    def test_cutoff_stability_failure(self) -> None:
        p = closing_payload(); p["cutoff_stability_passed"] = False
        self.assertRejected(p)

    def test_independent_verifier_total_mismatch(self) -> None:
        p = closing_payload(); p["global_defect_upper"] = "0.008"
        self.assertRejected(p)


if __name__ == "__main__":
    unittest.main()
