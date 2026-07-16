from __future__ import annotations

from pathlib import Path
import unittest

from flint import acb, arb, ctx

from continuum_defect_arb import parse_trial
from track_b_direct_weighted_arb import (
    floor_cell_taylor_model,
    floor_residual_param_jet,
    projected_symmetry_algebra_certificate,
)
from track_b_floor_taylor import TaylorModel, polynomial_l2_squared_upper
from track_b_overlap_arb import DifferentialTrial, action_jacobian, parse_matrix
from track_b_projected_mass_arb import projected_coefficients


class TrackBFloorTaylorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        ctx.prec = 192
        data, row = parse_trial(Path("six_copy_hejhal_balanced_coeffs.json"))
        cls.ev = DifferentialTrial(
            data["parameters"]["M"], str(row["r"]), row["coefficients"]
        )
        cls.odd = projected_coefficients(cls.ev)["odd"]
        cls.width = arb("0.30")
        cls.point = (arb("-0.4375"), arb("0.0625"), arb("0.01875"))
        x, t, s = cls.point
        cls.point_model = floor_cell_taylor_model(
            cls.ev, cls.odd, (x, x), (t, t), (s, s), cls.width, 4
        )
        cls.cell_bounds = (
            (-arb(1) / 2, -arb(3) / 8),
            (arb(0), arb(1) / 8),
            (arb(0), cls.width / 8),
        )
        cls.cell_model = floor_cell_taylor_model(
            cls.ev, cls.odd, *cls.cell_bounds, cls.width, 4
        )

    def test_coordinate_inversion_identity(self) -> None:
        x, t, s = arb("-0.2"), arb("0.3"), arb("0.1")
        rho = s.exp()
        y = (rho - x * x - t * t).sqrt()
        direct, _jac = action_jacobian(
            parse_matrix([[[0, 0], [-1, 0]], [[1, 0], [0, 0]]]), x, t, y
        )
        formula = (-x / rho, t / rho, y / rho)
        for actual, expected in zip(direct, formula):
            difference = actual - expected
            self.assertTrue(difference.contains(0))
            self.assertLess(abs(difference).upper(), arb("1e-40"))
        self.assertTrue((y * y - (rho - x * x - t * t)).contains(0))

    def test_measure_identity(self) -> None:
        x, t, s = arb("-0.2"), arb("0.3"), arb("0.1")
        rho = s.exp()
        y2 = rho - x * x - t * t
        y = y2.sqrt()
        dy_ds = rho / (2 * y)
        transformed = dy_ds / (y ** 3)
        formula = rho / (2 * y2 * y2)
        difference = transformed - formula
        self.assertTrue(difference.contains(0))
        self.assertLess(abs(difference).upper(), arb("1e-40"))

    def test_point_commutator_agrees_with_legacy_jet(self) -> None:
        x, t, s = self.point
        legacy, _gradient = floor_residual_param_jet(
            self.ev, self.odd, x, t, s, self.width
        )
        modeled = self.point_model["residual_model"].evaluate(
            (arb(0), arb(0), arb(0))
        )
        # The requested convention uses D=W o S-W, the negative of legacy D.
        difference = modeled + legacy
        self.assertTrue(difference.real.contains(0))
        self.assertTrue(difference.imag.contains(0))
        self.assertLess(abs(difference).upper(), arb("1e-35"))

    def test_projected_zero_relations_are_algebraic(self) -> None:
        result = projected_symmetry_algebra_certificate(self.ev, self.odd)
        self.assertTrue(result["certified"])
        self.assertTrue(result["laplace_eigen_equation"]["certified"])
        for name in ("T1", "R", "TiR"):
            self.assertEqual(result["relations"][name]["value_defect"], "0")
            self.assertEqual(
                result["relations"][name]["first_gradient_defect"], "0"
            )

    def test_legendre_polynomial_integral(self) -> None:
        constant = TaylorModel.constant(4, 1)
        q, _count = polynomial_l2_squared_upper(
            constant, arb(1), arb(1), arb(1)
        )
        self.assertGreaterEqual(q.lower(), arb(8))

        x = TaylorModel.variable(4, arb(0), arb(1), 0)
        q, _count = polynomial_l2_squared_upper(x, arb(1), arb(1), arb(1))
        self.assertGreaterEqual(q.lower(), arb(8) / 3)

        model = TaylorModel(4, {
            (0, 0, 0): acb(arb(1) / 3),
            (1, 0, 0): acb(arb(2) / 5),
            (0, 2, 1): acb(-arb(1) / 7),
        }, arb(0))
        certified, _count = polynomial_l2_squared_upper(
            model, arb(1), arb(1), arb(1)
        )
        exact = acb(0)
        for alpha, ca in model.coeff.items():
            for beta, cb in model.coeff.items():
                weight = arb(1)
                for axis in range(3):
                    power = alpha[axis] + beta[axis]
                    if power % 2:
                        weight = arb(0)
                        break
                    weight *= arb(2) / (power + 1)
                exact += ca * cb.conjugate() * weight
        self.assertTrue(exact.imag.contains(0))
        self.assertLessEqual(exact.real.upper(), certified.upper())

    def test_remainder_contains_direct_interior_values(self) -> None:
        xi = (arb("0.2"), arb("-0.3"), arb("0.4"))
        coordinates = []
        for bounds, normalized in zip(self.cell_bounds, xi):
            center = (bounds[0] + bounds[1]) / 2
            halfwidth = (bounds[1] - bounds[0]) / 2
            coordinates.append(center + halfwidth * normalized)
        x, t, s = coordinates
        legacy, _gradient = floor_residual_param_jet(
            self.ev, self.odd, x, t, s, self.width
        )
        rho = s.exp()
        y2 = rho - x * x - t * t
        sqrt_j = (s / 2).exp() / (arb(2).sqrt() * y2)
        direct = -arb(6).sqrt() * sqrt_j * legacy
        enclosure = self.cell_model["model"].evaluate(xi)
        difference = enclosure - direct
        self.assertTrue(difference.real.contains(0))
        self.assertTrue(difference.imag.contains(0))

    def test_no_bessel_fallback(self) -> None:
        self.assertEqual(self.cell_model["bessel_fallback_count"], 0)

    def test_tensor_degree_retains_mixed_index(self) -> None:
        x = TaylorModel.variable((4, 4, 1), arb(0), arb(1), 0)
        t = TaylorModel.variable((4, 4, 1), arb(0), arb(1), 1)
        s = TaylorModel.variable((4, 4, 1), arb(0), arb(1), 2)
        model = (x ** 4) * (t ** 4) * s
        self.assertIn((4, 4, 1), model.coeff)
        self.assertTrue(model.error.contains(0))


if __name__ == "__main__":
    unittest.main()
