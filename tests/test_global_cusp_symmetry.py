"""Acceptance tests for the global cusp-symmetry work (branch research/global-cusp-symmetry).

Tiers: tests in classes ``Exact*`` check PROVED exact identities; ``Experimental*`` check
sample-point / floating diagnostics; ``FailClosed*`` check that no certification flag moved.
"""
from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import random
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import global_cusp_character_algebra as ca                       # noqa: E402
import global_cusp_gaussian as gg                                 # noqa: E402
import global_quarter_turn as qt                                  # noqa: E402
import global_conjugate_prime as cp                               # noqa: E402
import global_oldspace as oldsp                                   # noqa: E402
import global_trial_compatibility as tc                           # noqa: E402

try:
    import flint  # noqa: F401
    HAVE_FLINT = True
except ImportError:                                               # pragma: no cover
    HAVE_FLINT = False


class ExactGaussianPrimeTests(unittest.TestCase):
    def test_prime_reduction_and_conjugation(self):
        self.assertEqual(gg.gmul(gg.PI, gg.PI_BAR), gg.g(5, 0))
        self.assertTrue(gg.divisible_by(gg.g(-3, 1), gg.PI))          # i = 3 mod (2+i)
        self.assertTrue(gg.divisible_by(gg.g(-2, 1), gg.PI_BAR))      # i = 2 mod (2-i)
        self.assertFalse(gg.divisible_by(gg.PI_BAR, gg.PI))
        res = cp.residue_field_isomorphism()
        self.assertTrue(all(res.values()), res)


class ExactCharacterAlgebraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = ca.compute_certificate()

    def test_everything_passes(self):
        self.assertTrue(self.cert["all_checks_pass"])

    def test_conventions_match_repository(self):
        data = json.loads((ROOT / "six_copy_DK_algebra_result.json").read_text())
        for name, perm in ca.FWD.items():
            self.assertEqual(list(perm), data["forward_actions"][name])
        self.assertTrue(all(self.cert["first_principles_perms_match_repository"].values()))

    def test_order_60_group_is_A5_doubly_transitive(self):
        self.assertEqual(self.cert["group_order"], 60)
        self.assertEqual(self.cert["class_sizes"], [1, 12, 12, 15, 20])
        self.assertTrue(self.cert["doubly_transitive"])
        self.assertEqual(self.cert["inner_product_chi_chi"], "2")

    def test_zeta_sums_exactly_zero(self):
        self.assertEqual(ca.zeta_sums(), {0: True, 1: True, 2: True, 3: True, 4: True})

    def test_rank_two_projector_and_noncommutation(self):
        self.assertEqual(self.cert["A_T_rank"], 2)
        self.assertEqual(self.cert["I_minus_A_T_rank"], 4)
        self.assertTrue(self.cert["A_T_idempotent"] and self.cert["A_T_symmetric"])
        self.assertFalse(self.cert["A_T_commutes_with"]["S"])
        self.assertTrue(self.cert["commutator_A_T_rho_S_nonzero"])
        self.assertFalse(self.cert["V4_group_invariant"])

    def test_commutant_is_span_I_J_and_obstruction(self):
        self.assertEqual(self.cert["commutant_dimension"], 2)
        self.assertTrue(self.cert["commutant_equals_span_I_J"])
        self.assertEqual(self.cert["constant_equivariant_annihilators_of_singular_subspace"], 0)

    def test_exact_arithmetic_has_teeth(self):
        # a wrong k must fail: sum zeta^{2j} over j but with wrong exponent pattern is not zero
        s = ca.Zeta5.const(0)
        for j in range(5):
            s = s + ca.Zeta5.zeta_pow(j * 0)
        self.assertFalse(s.is_zero())


class ExactQuarterTurnTests(unittest.TestCase):
    def test_normalization_and_square(self):
        res = qt.check_normalization()
        self.assertTrue(all(res.values()), res)

    def test_point_actions_and_cusp_rotations(self):
        res = qt.check_point_actions()
        self.assertTrue(all(res.values()), res)

    def test_cusp_lattices_and_fourier_rotation(self):
        self.assertTrue(all(v for k, v in qt.check_cusp_lattices().items() if isinstance(v, bool)))
        self.assertTrue(all(qt.check_fourier_rotation_identity().values()))

    def test_bundle_level_map_not_just_fiber_permutation(self):
        res = qt.check_bundle_action()
        for key in ("intertwining_on_random_words", "rho_q_squared_times_rho_R_is_identity",
                    "finite_q_squared_is_nontrivial_so_fiber_map_alone_is_not_an_involution"):
            self.assertTrue(res[key], key)
        self.assertTrue(all(res["intertwining_on_generators"].values()))
        self.assertTrue(all(res["generator_images"].values()))

    def test_matrix_level_convention_is_T1_squared_not_cubed(self):
        c = ca.compute_certificate()
        self.assertTrue(c["q_T1_qinv_equals_T1_cubed_as_permutations"])
        self.assertTrue(c["matrix_level_rho_q_rho_T1_rho_qinv"]["equals_rho_T1_squared"])
        self.assertFalse(c["matrix_level_rho_q_rho_T1_rho_qinv"]["equals_rho_T1_cubed"])

    def test_wrong_permutation_is_rejected(self):
        """Negative control: the inverse of q does NOT intertwine; neither does the identity."""
        rng = random.Random(5)
        words = [gg.random_word(rng, 6) for _ in range(40)]
        q = ca.quarter_turn_perm("pi")
        for bad in (ca.inverse(q), tuple(range(6))):
            rq, rqi = ca.rho_matrix(bad), ca.rho_matrix(ca.inverse(bad))
            ok = all(
                ca.matmul(ca.matmul(rq, qt.rho_of(qt.conj_by_Q(m))), rqi) == qt.rho_of(m) for m in words
            )
            self.assertFalse(ok)

    def test_mode_orbit_accounting(self):
        a4 = qt.accounting(4)
        self.assertEqual((a4["unknowns_full"], a4["unknowns_after_R_only"], a4["unknowns_per_Q_parity_sector"]),
                         (80, 40, 20))
        a1 = qt.accounting(1)
        self.assertEqual((a1["unknowns_full"], a1["unknowns_per_Q_parity_sector"]), (24, 6))
        for m in range(1, 7):
            self.assertTrue(qt.accounting(m)["orbits_all_size_4"])

    def test_sector_columns_encode_sign_pattern(self):
        cols = qt.sector_basis(1, -1)
        mi, _ = qt.cusp_modes(1)
        col = cols[0]                              # infinity orbit of norm 1
        self.assertEqual(sorted(col.values()), [-1, -1, 1, 1])
        # a_{i nu} = -a_nu : pair nu=(1,0) with i nu=(0,1)
        idx = {(m[0], m[1]): k for k, m in enumerate(mi)}
        self.assertEqual(col[idx[(1, 0)]], -col[idx[(0, 1)]])
        self.assertEqual(col[idx[(1, 0)]], col[idx[(-1, 0)]])        # R-even


class ExactConjugatePrimeTests(unittest.TestCase):
    def test_conjugate_prime_correspondence(self):
        cert = cp.compute_certificate()
        self.assertTrue(cert["all_checks_pass"])
        fl = cert["fiber_level"]
        self.assertTrue(fl["Gaussian_conjugation_is_not_a_cyclotomic_Galois_relabeling"])
        self.assertTrue(fl["quarter_turn_perm_pibar_is_inverse_of_pi"])


class ExactOldspaceTests(unittest.TestCase):
    def test_both_lifts_and_commutation(self):
        cert = oldsp.compute_certificate()
        self.assertTrue(cert["all_checks_pass"], cert["checks"])
        self.assertIn("newness of any quasimode", cert["not_certified"])


@unittest.skipUnless(HAVE_FLINT, "python-flint not installed")
class ExperimentalFieldLevelTests(unittest.TestCase):
    def test_sector_relations_on_repository_field_map(self):
        res = qt.field_level_intertwining_check(cutoff=2, points=2)
        self.assertTrue(res["passes"], res)

    def test_repository_mode_lists_agree(self):
        from track_b_two_cusp_data import cusp_modes
        for m in range(1, 7):
            self.assertEqual(cusp_modes(m), qt.cusp_modes(m))

    def test_exact_gaussian_action_agrees_with_arb_action(self):
        from flint import arb, ctx
        from track_b_two_cusp_data import h3_action as arb_action, GENERATORS as ARB_GENS
        ctx.prec = 128
        for name, gen in gg.GENERATORS.items():
            for (z, y) in gg.random_points(4, 3):
                ex = gg.h3_action(gen, (z, y))
                ab = arb_action(ARB_GENS[name], arb(z[0].numerator) / z[0].denominator,
                                arb(z[1].numerator) / z[1].denominator, arb(y.numerator) / y.denominator)
                self.assertTrue(ab[0].contains(arb(ex[0][0].numerator) / ex[0][0].denominator) or
                                abs(ab[0] - arb(ex[0][0].numerator) / ex[0][0].denominator) < arb(10) ** -30)
                self.assertTrue(abs(ab[2] - arb(ex[1].numerator) / ex[1].denominator) < arb(10) ** -30)


class FailClosedTests(unittest.TestCase):
    """No certification flag may move, and incompatible trials must be reported as such."""

    def test_trial_compatibility_fails_closed(self):
        out = tc.audit()
        self.assertFalse(out["compatible"])
        self.assertFalse(out["same_trial"])
        self.assertFalse(out["same_spectral_parameter"])
        self.assertTrue(any("floor certificate carries no trial hash" in b for b in out["blockers"]))

    def test_existing_false_flags_preserved(self):
        d = json.loads((ROOT / "track_b_global_hejhal_defect_result.json").read_text())
        for flag in ("global_hejhal_defect_certified", "rung4_certified", "hejhal_existence_certified",
                     "both_cusp_continuum_certified", "all_nonfloor_faces_certified",
                     "fourier_tails_certified", "reprojection_certified", "global_threshold_certified"):
            self.assertIs(d[flag], False, flag)

    def test_recorded_mass_hash_is_crlf_hash(self):
        out = tc.audit()
        if out["hash_note"] is not None:
            self.assertTrue(out["hash_note"]["recorded_mass_trial_sha_matches_CRLF_bytes"])

    def test_no_new_module_writes_a_certification_flag(self):
        for name in ("global_cusp_character_algebra", "global_quarter_turn", "global_conjugate_prime",
                     "global_oldspace", "global_trial_compatibility", "global_cusp_symmetry_experiment"):
            text = (ROOT / f"{name}.py").read_text()
            self.assertNotIn('"rung4_certified": True', text)
            self.assertNotIn("global_hejhal_defect_certified\": True", text)


if __name__ == "__main__":
    unittest.main()
