from fractions import Fraction
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from flint import arb, acb, ctx
from core.assemble import evaluate
from core.bspline import g0, gpp0, second_derivative_coefficient
from core.certificate import certificate_payload, report_payload
from core.terms import prime_term
from core.testfunctions import sinc2k
from fields.quadratic import get_field, Lval1
from groups import get_group
from groups.systoles import verify_all

ROOT = Path(__file__).resolve().parents[1]


class TraceCoreTests(unittest.TestCase):
    def setUp(self):
        ctx.prec = 100

    def test_exact_second_derivative_known_values(self):
        # Independent low-degree spline identities, including the quartic and
        # sextic corrections that the previous cubic stencil missed.
        for k, c in [(2, Fraction(-1, 4)), (3, Fraction(-1, 8)), (4, Fraction(-1, 12))]:
            self.assertEqual(second_derivative_coefficient(k), c)
            expected = arb(c.numerator)/c.denominator/arb("0.2")**3
            self.assertTrue(gpp0(k, "0.2").overlaps(expected))

    def test_fourier_transform_at_zero(self):
        # integral sinc^4(delta*r) dr/(2*pi)=1/(3*delta).
        self.assertTrue(g0(2, "0.2").overlaps(arb(1)/(3*arb("0.2"))))

    def test_entire_sinc_at_zero_and_complex_arguments(self):
        self.assertTrue(sinc2k(0, 2, "0.2").contains(1))
        z = acb(1, 2)*arb("0.2")
        self.assertTrue(sinc2k(acb(1, 2), 3, "0.2").overlaps((z.sin()/z)**6))

    def test_all_systoles_and_class_numbers(self):
        rows = verify_all()
        self.assertEqual(len(rows), 9)
        self.assertEqual([r["checked"] for r in rows], [24,18,32,20,14,10,2,2,2])

    def test_closed_form_eisenstein_systole(self):
        v = (1+arb(21).sqrt())/4
        self.assertTrue(get_group("omega").systole().overlaps((v+(v*v-1).sqrt()).log()))

    def test_character_and_splitting(self):
        self.assertTrue(Lval1("i").overlaps(arb.pi()/4))
        self.assertTrue(Lval1("omega").overlaps(arb.pi()/(3*arb(3).sqrt())))
        self.assertEqual(get_field(2).prime_norms(2), (2, 1))
        self.assertEqual(get_field(7).prime_norms(2), (2, 2))
        self.assertEqual(get_field(11).prime_norms(2), (4, 1))
        self.assertEqual(get_field(19).prime_norms(3), (9, 1))

    def test_prime_term_baseline_and_multiplicity(self):
        from core.bspline import gbspline
        delta = arb("0.24")
        expected = arb(2).log()/2*gbspline(arb(2).log(), 2, delta)
        self.assertTrue(prime_term(1, 2, delta).overlaps(expected))
        self.assertTrue(prime_term(3, 2, delta).is_zero())
        self.assertTrue(prime_term(7, 2, delta).overlaps(2*expected))

    def test_inert_prime_and_prime_powers(self):
        from core.bspline import gbspline
        delta = arb("0.4")  # support 1.6: norm 4 survives, norm 9 does not.
        expected = arb(4).log()/4*gbspline(arb(4).log(), 2, delta)
        self.assertTrue(prime_term(19, 2, delta).overlaps(expected))
        # In Q(sqrt(-7)), two norm-2 primes contribute at norms 2 and 4.
        expected = arb(2).log()*gbspline(arb(2).log(), 2, delta)
        expected += arb(2).log()/2*gbspline(arb(4).log(), 2, delta)
        self.assertTrue(prime_term(7, 2, delta).overlaps(expected))

    def test_invalid_parameters(self):
        for k in [1, 0, -1, 2.5, True]:
            with self.assertRaises(ValueError):
                evaluate("i", k=k, verbose=False)
        for frac in [0, -1, 1.01, float("nan")]:
            with self.assertRaises(ValueError):
                evaluate("i", frac=frac, verbose=False)
        for R in [6.99, float("inf")]:
            with self.assertRaises(ValueError):
                evaluate("i", R=R, verbose=False)

    def test_unknown_fields_do_not_fall_back_to_omega(self):
        with self.assertRaises(ValueError):
            get_group("typo")

    def test_baseline_enclosures_and_proof_gate(self):
        # Historical six-decimal endpoints expanded outwards. Quadrature
        # changes may tighten enclosures; containment in this envelope is the
        # regression criterion, not identical numerical radii.
        for kind, lo, hi in [("i", "0.30455", "0.31784"),
                             ("omega", "0.52519", "0.54371")]:
            E = evaluate(kind, verbose=False)
            self.assertTrue(E.bound.lower() > arb(lo))
            self.assertTrue(E.bound.upper() < arb(hi))
            self.assertTrue(2*E.k*arb(E.delta) < E.group.systole())
            self.assertFalse(report_payload(E)["spectral_certificate"])
            with self.assertRaises(ValueError):
                certificate_payload(E)

    def test_new_fields_block_full_assembly(self):
        for d in [2, 7, 11, 19]:
            with self.assertRaisesRegex(ValueError, "inventory is incomplete"):
                evaluate(d, verbose=False)

    def test_mechanical_screen_cannot_be_a_certificate(self):
        E = evaluate(2, verbose=False, include_elliptic=False)
        self.assertTrue(E.bound.lower() > arb("-0.44"))
        self.assertTrue(E.bound.upper() < arb("-0.40"))
        self.assertTrue(E.terms["NCE"].is_zero())
        self.assertFalse(report_payload(E)["spectral_certificate"])
        with self.assertRaises(ValueError):
            certificate_payload(E)

    def test_corrected_k3_omega_bound(self):
        E = evaluate("omega", k=3, verbose=False)
        self.assertTrue(E.bound.lower() > arb("0.75"))
        self.assertTrue(E.bound.upper() < arb("0.82"))

    def test_flip_script_from_unrelated_working_directory(self):
        with tempfile.TemporaryDirectory() as cwd:
            script = ROOT/"scripts/feasibility/flip_check.py"
            copied = Path(cwd)/"flip_check.py"
            copied.write_text(script.read_text())
            for source in [script, copied]:
                result = subprocess.run([sys.executable, str(source), "--repo", str(ROOT)],
                    cwd=cwd, capture_output=True, text=True, check=True)
                self.assertIn("closed Klein-four subgroup of order 4", result.stdout)
                self.assertIn("NOT PROVED: maximality", result.stdout)


if __name__ == "__main__":
    unittest.main()
