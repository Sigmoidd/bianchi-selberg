from __future__ import annotations

import copy
import json
from pathlib import Path
import unittest

from track_b_collocation_family import (
    candidate_point_checks, exact_candidate_points, phase_audit,
    transformed_height_audit,
)
from track_b_observable_stability import compare, freeze_gate


ROOT = Path(__file__).resolve().parent


class TrackBObservableCollocationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = {
            m: json.loads((ROOT / f"track_b_observable_M{m}_result.json").read_text(encoding="utf-8"))
            for m in (2, 3, 4, 5)
        }
        cls.stability = json.loads(
            (ROOT / "track_b_observable_stability_result.json").read_text(encoding="utf-8")
        )

    def assertFreezeClosed(self, hard: dict[str, bool]) -> None:
        self.assertFalse(all(hard.values()))
        self.assertFalse(self.stability["trial_frozen"])
        self.assertFalse(self.stability["global_hejhal_defect_certified"])
        self.assertFalse(self.stability["rung4_certified"])
        self.assertFalse(self.stability["dual_certification"])

    def test_exact_phase_aliases_are_detected(self) -> None:
        bad = [{
            "point_id": "bad", "x1": "0", "x2": "0", "y": "6/5",
            "design": "exact rational non-symmetric Humbert-core point",
        }]
        self.assertGreater(phase_audit(4, bad)["exact_alias_count"], 0)

    def test_new_family_has_no_exact_phase_alias(self) -> None:
        self.assertEqual(phase_audit(5, exact_candidate_points(20))["exact_alias_count"], 0)

    def test_identical_transformed_height_signatures_are_detected(self) -> None:
        bad = [
            {"point_id": "a", "x1": "1/5", "x2": "2/7", "y": "6/5"},
            {"point_id": "b", "x1": "-1/5", "x2": "-2/7", "y": "6/5"},
        ]
        self.assertFalse(
            transformed_height_audit(bad)["different_source_points_have_distinct_height_signatures"]
        )

    def test_old_M4_nullspace_is_high_shell_dominated(self) -> None:
        audit = json.loads((ROOT / "track_b_M4_nullspace_audit.json").read_text(encoding="utf-8"))
        self.assertEqual(audit["observed_numerical_nullity"], 7)
        self.assertGreater(min(q["fraction_in_new_M4_modes"] for q in audit["null_vectors"]), 0.7)

    def test_deficient_full_candidate_matrix_rejects_freeze(self) -> None:
        results = copy.deepcopy(list(self.results.values()))
        results[2]["observability_diagnostics"]["candidate_normalized_numerical_rank"] -= 1
        hard = freeze_gate(results, self.stability["comparisons"])
        self.assertFreezeClosed(hard)

    def test_selected_rows_fit_but_omitted_rows_fail(self) -> None:
        results = copy.deepcopy(list(self.results.values()))
        results[2]["independent_oversampled_verification"] = False
        hard = freeze_gate(results, self.stability["comparisons"])
        self.assertFreezeClosed(hard)

    def test_normalization_change_between_cutoffs_is_rejected(self) -> None:
        high = copy.deepcopy(self.results[4]); high["normalization_hash"] = "changed"
        comparison = compare(self.results[3], high)
        self.assertFalse(comparison["certified"])

    def test_unrecorded_complex_phase_alignment_is_rejected(self) -> None:
        comparisons = copy.deepcopy(self.stability["comparisons"])
        comparisons[1]["unrecorded_phase_alignment"] = True
        hard = freeze_gate(list(self.results.values()), comparisons)
        self.assertFreezeClosed(hard)

    def test_duplicate_candidate_points_are_rejected(self) -> None:
        points = exact_candidate_points(2)
        points.append(copy.deepcopy(points[0]))
        points[-1]["point_id"] = "new-id-same-coordinate"
        self.assertFalse(candidate_point_checks(points)["candidate_coordinates_exact_once"])

    def test_nonphysical_synthetic_enrichment_row_is_rejected(self) -> None:
        results = copy.deepcopy(list(self.results.values()))
        results[2]["enrichment_rows_all_exact_physical"] = False
        hard = freeze_gate(results, self.stability["comparisons"])
        self.assertFreezeClosed(hard)

    def test_bessel_fallback_is_rejected(self) -> None:
        results = copy.deepcopy(list(self.results.values()))
        results[2]["resolved_bessel_fallback_count"] = 1
        hard = freeze_gate(results, self.stability["comparisons"])
        self.assertFreezeClosed(hard)

    def test_delta_failure_keeps_trial_unfrozen(self) -> None:
        self.assertGreaterEqual(float(self.stability["delta_3_4_upper"].split()[0].lstrip("[")), 0.25)
        self.assertFalse(self.stability["trial_frozen"])

    def test_old_floor_certificate_is_not_imported(self) -> None:
        compatibility = json.loads(
            (ROOT / "track_b_legacy_trial_compatibility.json").read_text(encoding="utf-8")
        )
        floor = next(q for q in compatibility["artifacts"] if q["artifact"] == "floor")
        self.assertTrue(floor["valid_for_original_trial"])
        self.assertFalse(floor["compatible_with_current_trial"])
        self.assertFalse(compatibility["legacy_artifacts_imported_into_current_trial"])


if __name__ == "__main__":
    unittest.main()
