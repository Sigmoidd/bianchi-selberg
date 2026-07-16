from __future__ import annotations

import copy
import json
from pathlib import Path
import unittest

from flint import arb, ctx

from track_b_global_partition_arb import (
    FLOOR_WIDTH,
    RealJet2,
    _stabilizer_certificate,
    bounds_ball,
    certification_decision,
    certify_partition,
    hyperbolic_gradient_norm,
    hyperbolic_gradient_upper,
    hyperbolic_laplacian,
    log_rho_jet_from_parameter_box,
    normalized_gradient_formula,
    normalized_laplacian_formula,
    raw_partition_jets,
    smoothstep_clamped,
)
from track_b_global_partition_verify import verify


class TrackBGlobalPartitionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        ctx.prec = 192
        cls.geometry = json.loads(Path("track_b_partition_result.json").read_text())

    def test_bounds_ball_preserves_nonzero_endpoints(self) -> None:
        box = bounds_ball(arb("0.7"), arb("1.2"))
        self.assertGreater(box.lower(), arb(0))
        self.assertLessEqual(box.lower(), arb("0.7"))
        self.assertGreaterEqual(box.upper(), arb("1.2"))

    def test_exact_tensor_partition_identity(self) -> None:
        x, t, s = arb("0.2"), arb("0.2"), arb("0.1")
        y = (s.exp() - x * x - t * t).sqrt()
        raw, _gates = raw_partition_jets(x, t, y, s)
        total = sum((jet for _name, jet in raw), RealJet2.constant(0))
        self.assertTrue(total.v.contains(1))
        self.assertTrue(all(q.contains(0) for q in total.g))
        self.assertTrue(all(q.contains(0) for row in total.h for q in row))

    def test_project_laplacian_quotient_rule(self) -> None:
        x = RealJet2.variable(arb("0.2"), 0)
        t = RealJet2.variable(arb("0.1"), 1)
        y = RealJet2.variable(arb("1.1"), 2)
        phi = 1 + x * x + t * y
        Phi = 2 + x + y * y
        quotient = phi / Phi
        gradient_formula = normalized_gradient_formula(phi, Phi)
        for actual, expected in zip(quotient.g, gradient_formula):
            self.assertTrue((actual - expected).contains(0))
            self.assertLess(abs(actual - expected).upper(), arb("1e-45"))
        direct = hyperbolic_laplacian(quotient, y.v)
        formula = normalized_laplacian_formula(phi, Phi, y.v)
        self.assertTrue((direct - formula).contains(0))
        self.assertLess(abs(direct - formula).upper(), arb("1e-45"))

    def test_metric_conversion_from_s_coordinate(self) -> None:
        x, t, y = arb("0.2"), arb("0.1"), arb("1.0")
        xj = RealJet2.variable(x, 0)
        tj = RealJet2.variable(t, 1)
        yj = RealJet2.variable(y, 2)
        direct = (xj * xj + tj * tj + yj * yj).log()
        s = (x * x + t * t + y * y).log()
        converted = log_rho_jet_from_parameter_box(x, t, y, s)
        for actual, expected in zip(direct.g, converted.g):
            self.assertTrue((actual - expected).contains(0))
        self.assertTrue((hyperbolic_gradient_norm(direct, y)
                         - hyperbolic_gradient_norm(converted, y)).contains(0))
        self.assertTrue((hyperbolic_laplacian(direct, y)
                         - hyperbolic_laplacian(converted, y)).contains(0))

    def test_exact_stabilizer_averaging_ledger(self) -> None:
        certificate = _stabilizer_certificate(self.geometry)
        self.assertTrue(certificate["certified"])
        self.assertEqual(certificate["edge_orders"], [3, 3, 2, 3])
        self.assertTrue(certificate["global_stratum_group_reindexing_exact"])
        self.assertTrue(all(certificate["vertex_group_reindexing_exact"].values()))
        self.assertFalse(certificate["sampled_equality_used"])

    def test_floor_weight_formula_consistency(self) -> None:
        x, t, s = arb("-0.4"), arb("0.45"), arb("0.12")
        y = (s.exp() - x * x - t * t).sqrt()
        _raw, gates = raw_partition_jets(x, t, y, s)
        global_core = 1 - gates["floor"]
        local = smoothstep_clamped(RealJet2.constant(s / arb(FLOOR_WIDTH)))
        self.assertTrue((global_core.v - local.v).contains(0))

    def test_each_global_gate_fails_closed(self) -> None:
        conditions = {
            "coverage_certified": True,
            "denominator_positive_certified": True,
            "partition_sum_certified": True,
            "stabilizer_averages_certified": True,
            "floor_weight_consistency_certified": True,
            "weight_gradients_certified": True,
            "weight_laplacians_certified": True,
            "fallback_zero": True,
        }
        self.assertTrue(certification_decision(conditions)["global_weight_bounds_certified"])
        for name in conditions:
            broken = dict(conditions)
            broken[name] = False
            decision = certification_decision(broken)
            self.assertFalse(decision["global_weight_bounds_certified"], name)
            self.assertFalse(decision["rung4_certified"], name)

    def test_independent_verifier_rejects_ledger_faults(self) -> None:
        result = certify_partition(self.geometry, (2, 1, 2), 8, 192, None)
        # Reconstruct the same four records deterministically through a small
        # temporary audit generated by the public certificate function.
        audit = Path("track_b_partition_verify_test.jsonl")
        result = certify_partition(self.geometry, (2, 1, 2), 8, 192, audit)
        records = [json.loads(line) for line in audit.read_text().splitlines()]
        result["stability_check_passed"] = True
        result["provisional"] = False
        self.assertTrue(verify(result, records)["verified"])

        mutations = []
        mutations.append(records[:-1])                         # missing region
        mutations.append(records + [copy.deepcopy(records[0])])  # duplicate
        denominator = copy.deepcopy(records)
        denominator[0]["Phi_lower"] = "0"
        mutations.append(denominator)
        fallback = copy.deepcopy(records)
        fallback[0]["fallback_count"] = 1
        mutations.append(fallback)
        floor = copy.deepcopy(records)
        floor[0]["floor_weight_consistency"] = False
        mutations.append(floor)
        unbounded = copy.deepcopy(records)
        unbounded[0]["maximum_grad_chi_upper"] = "+inf"
        mutations.append(unbounded)
        for bad in mutations:
            self.assertFalse(verify(result, bad)["verified"])

        missing_stabilizer = copy.deepcopy(result)
        missing_stabilizer["stabilizer_certificate"][
            "vertex_group_reindexing_exact"
        ]["v_00"] = False
        self.assertFalse(verify(missing_stabilizer, records)["verified"])
        audit.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
