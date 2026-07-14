from __future__ import annotations

import math
import unittest

from flint import arb, ctx

from verified_kir import kir_enclosure
from six_copy_hejhal import DUAL_SCALE, gaussian_modes, mode_frequencies
from continuum_defect_arb import act


class VerifiedKirTests(unittest.TestCase):
    def test_reference_value(self):
        q = kir_enclosure(6.7439020359331625, 5.026548245743669, 160)
        ref = 2.8917910552406475872622790302334687893494002649448557e-5
        self.assertLessEqual(q.lower, ref)
        self.assertGreaterEqual(q.upper, ref)
        self.assertLess(q.radius, 1e-40)

    def test_r_dependence(self):
        a = kir_enclosure(6.0, 2 * math.pi * 0.8, 128)
        b = kir_enclosure(8.0, 2 * math.pi * 0.8, 128)
        self.assertNotEqual(a.midpoint, b.midpoint)

    def test_fail_closed_domain(self):
        with self.assertRaises(ValueError):
            kir_enclosure(6.0, 0.0, 128)


class CuspLatticeTests(unittest.TestCase):
    def test_dual_scale(self):
        self.assertAlmostEqual(DUAL_SCALE.real, 0.4)
        self.assertAlmostEqual(DUAL_SCALE.imag, 0.2)
        # mu * conj(2+i) is a Gaussian integer for every indexed beta.
        for mu in mode_frequencies(gaussian_modes(12), cusp=1):
            beta = mu * complex(2, -1)
            self.assertAlmostEqual(beta.real, round(beta.real), places=12)
            self.assertAlmostEqual(beta.imag, round(beta.imag), places=12)

    def test_balanced_physical_cutoff(self):
        M = 12
        inf = mode_frequencies(gaussian_modes(M), cusp=0)
        zero = mode_frequencies(gaussian_modes(5 * M), cusp=1)
        self.assertLessEqual(max(abs(x) ** 2 for x in inf), M + 1e-12)
        self.assertLessEqual(max(abs(x) ** 2 for x in zero), M + 1e-12)
        # Every infinity-lattice frequency is present in the zero-cusp dual set.
        zset = {(round(z.real, 12), round(z.imag, 12)) for z in zero}
        for z in inf:
            self.assertIn((round(z.real, 12), round(z.imag, 12)), zset)


class IntervalActionTests(unittest.TestCase):
    def setUp(self):
        ctx.prec = 160

    def assert_contains_float(self, ball, value):
        self.assertLessEqual(float(ball.lower()), value)
        self.assertGreaterEqual(float(ball.upper()), value)

    def test_exact_generator_actions(self):
        x1, x2, y = arb("0.1"), arb("0.2"), arb("0.8")
        expected = {
            "T1": (1.1, 0.2, 0.8),
            "R": (-0.1, -0.2, 0.8),
            "TiR": (-0.1, 0.8, 0.8),
            "S": (-0.1 / 0.69, 0.2 / 0.69, 0.8 / 0.69),
        }
        for name, want in expected.items():
            got = act(name, x1, x2, y)
            for ball, value in zip(got, want):
                self.assertAlmostEqual(float(ball.mid()), value, places=14)

    def test_S_box_contains_corner_images(self):
        xb = arb("0.15", "0.05")
        zb = arb("-0.25", "0.05")
        yb = arb("0.9", "0.1")
        gx, gz, gy = act("S", xb, zb, yb)
        for x in (0.1, 0.2):
            for z in (-0.3, -0.2):
                for y in (0.8, 1.0):
                    den = x * x + z * z + y * y
                    self.assert_contains_float(gx, -x / den)
                    self.assert_contains_float(gz, z / den)
                    self.assert_contains_float(gy, y / den)


if __name__ == "__main__":
    unittest.main()
