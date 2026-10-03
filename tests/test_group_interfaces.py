"""Extension contracts, exact levels and isolation from full-group formulas."""
from dataclasses import replace
import json
from pathlib import Path
import unittest

from flint import arb, ctx
from core.assemble import evaluate
from core.backends import AnalyticRegistry, Geometry
from core.backends.level_one import LevelOneBackend
from core.certificate import certificate_payload, report_payload
from groups import GroupKey, LevelIdeal, GroupRegistry, get_group
from groups.builtins import d2_backend
from groups.exact import Radical
from groups.arithmetic import QuadraticRing
from groups.matrix import IDENTITY, MatrixOps
from groups.relative_orders import QuadraticOrder, ROOT, power

ROOT_PATH = Path(__file__).resolve().parents[1]


class SyntheticCongruenceBackend:
    """Interface fixture only; supplies no congruence-group mathematical proof."""
    backend_id = "synthetic-test-only"
    documentation = ("tests/test_group_interfaces.py (synthetic fixture, not a proof)",)

    def validate(self, group):
        if group.key != GroupKey.congruence("gamma0", LevelIdeal.rational(2, 3)):
            raise ValueError("synthetic fixture supports only its exact test key")

    def geometry(self, group):
        self.validate(group)
        return Geometry(group.key, arb(1), arb(1))

    def terms(self, group, geometry, k, delta, R, include_elliptic):
        return {name: arb(0) for name in ("CE", "Ch0", "PARg0", "PSI", "PHIINT")}


class GroupInterfaceTests(unittest.TestCase):
    def setUp(self):
        ctx.prec = 100

    def test_distinct_prime_ideals_with_equal_norm(self):
        # The two primes above 2 in Q(sqrt(-7)) have distinct HNF and keys.
        P = LevelIdeal.principal(7, (0, 1))
        Q = LevelIdeal.principal(7, (1, -1))
        self.assertEqual(P.norm, Q.norm)
        self.assertEqual(P.norm, 2)
        self.assertNotEqual(P.hnf, Q.hnf)
        self.assertTrue(P.contains((0, 1)))
        self.assertFalse(Q.contains((0, 1)))
        self.assertNotEqual(GroupKey.congruence("gamma0", P), GroupKey.congruence("gamma0", Q))

    def test_associates_have_the_same_level_key(self):
        for d, generator in [(1, (2, 1)), (2, (0, 1)), (7, (1, -1)), (11, (3, 2)), (19, (-2, 3))]:
            I = LevelIdeal.principal(d, generator)
            J = LevelIdeal.principal(d, tuple(-c for c in generator))
            self.assertEqual(I, J)
            self.assertEqual(I.norm, QuadraticRing(d).norm(generator))
            self.assertTrue(I.contains(generator))
            self.assertTrue(I.contains(QuadraticRing(d).mul(generator, (0, 1))))

    def test_ideal_and_field_validation(self):
        for args in [(2, (3, 0, 1)), (2, (0, 0, 1)), (2, (2, 2, 1)), (2, (1, 0, True))]:
            with self.assertRaises(ValueError):
                LevelIdeal(*args)
        for d in [True, 2.0, "i"]:
            with self.assertRaises(ValueError):
                GroupKey(d)
        with self.assertRaises(ValueError):
            LevelIdeal.principal(2, (0, 0))
        with self.assertRaises(ValueError):
            GroupKey(2, "full", (3, 0, 3))
        self.assertEqual(GroupKey(2, "principal"), GroupKey(2))
        with self.assertRaisesRegex(ValueError, "key and field"):
            replace(get_group(2), key=GroupKey(7))

    def test_congruence_membership_respects_psl_sign(self):
        I = LevelIdeal.rational(2, 3)
        ops = MatrixOps(QuadraticRing(2))
        principal = GroupKey.congruence("principal", I)
        gamma0 = GroupKey.congruence("gamma0", I)
        gamma1 = GroupKey.congruence("gamma1", I)
        U = ((1, 0), (1, 0), (0, 0), (1, 0))
        V = ((1, 0), (0, 0), (3, 0), (1, 0))
        W = ((1, 0), (0, 0), (1, 0), (1, 0))
        for key in [principal, gamma0, gamma1]:
            self.assertTrue(key.contains(IDENTITY))
            self.assertTrue(key.contains(ops.neg(IDENTITY)))
            self.assertTrue(key.contains(V))
            self.assertTrue(key.contains(ops.neg(V)))
            self.assertFalse(key.contains(W))
        self.assertFalse(principal.contains(U))
        self.assertTrue(gamma0.contains(U))
        self.assertTrue(gamma1.contains(U))
        det_minus_one = ((1, 0), (0, 0), (0, 0), (-1, 0))
        self.assertFalse(gamma0.contains(det_minus_one))
        for key in [principal, gamma0, gamma1]:
            for A in [V, ops.neg(V)]:
                for B in [V, ops.neg(V)]:
                    self.assertTrue(key.contains(ops.mul(A, B)))

    def test_relative_order_arithmetic_is_reusable_for_other_fields(self):
        # Norm and multiplication determinant agree over d=7 and 11 too.
        for d in [7, 11, 19]:
            ring = QuadraticRing(d)
            order = QuadraticOrder(ring, (1, 0), (1, 0))
            for unit in [(1, 2, -1, 3), (0, 1, 2, 0)]:
                self.assertEqual(order.absolute_norm(unit), ring.norm(order.norm(unit)))
                self.assertEqual(MatrixOps(ring).det(order.matrix(unit)), order.norm(unit))
            self.assertEqual(power(order, ROOT, 3), (-1, 0, 0, 0))
        self.assertEqual(Radical(55, -12, 21).sign(), 1)
        self.assertEqual(Radical(3, -1, 10).sign(), -1)

    def test_proof_registry_is_bound_to_exact_group_and_id(self):
        G = get_group(2)
        registry = GroupRegistry()
        registry.register(G.key, lambda: G, d2_backend())
        self.assertEqual(get_group(G.key, registry).key, G.key)
        self.assertEqual(G.verify_inventory(registry)["class_count"], 4)
        with self.assertRaisesRegex(ValueError, "proof id"):
            replace(G, inventory_proof_id="different-proof").verify_inventory(registry)
        level = GroupKey.congruence("gamma0", LevelIdeal.rational(2, 3))
        with self.assertRaisesRegex(ValueError, "proof id"):
            replace(G, key=level).verify_inventory(registry)
        with self.assertRaises(ValueError):
            registry.register(GroupKey(7), lambda: get_group(7), d2_backend())
        with self.assertRaises(ValueError):
            registry.register(G.key, lambda: G)
        registry.register(G.key, lambda: replace(G, inventory_status="incomplete"), replace=True)
        with self.assertRaisesRegex(ValueError, "no self-contained"):
            registry.verify(G)

    def test_new_analytic_backend_dispatch_needs_no_core_field_branch(self):
        class DelegatingBackend(LevelOneBackend):
            backend_id = "level-one-dispatch-test"
        analytic = AnalyticRegistry()
        analytic.register(DelegatingBackend())
        G = replace(get_group(2), trace_backend_id=DelegatingBackend.backend_id)
        E = evaluate(G, verbose=False, analytic_registry=analytic)
        payload = certificate_payload(E, analytic_registry=analytic)
        self.assertTrue(payload["spectral_certificate"])
        self.assertEqual(payload["analytic_backend"], DelegatingBackend.backend_id)
        self.assertEqual(payload["group_identity"], GroupKey(2).payload())
        self.assertTrue(E.bound.upper() < arb("0.431"))

    def test_full_group_formulas_cannot_run_on_a_congruence_clone(self):
        key = GroupKey.congruence("gamma0", LevelIdeal.rational(2, 3))
        G = replace(get_group(2), key=key)
        for include in [False, True]:
            with self.assertRaisesRegex(ValueError, "full group"):
                evaluate(G, include_elliptic=include, verbose=False)
        with self.assertRaisesRegex(ValueError, "full group"):
            G.analytic_data(False)

    def test_congruence_backend_uses_common_pipeline_without_claiming_a_proof(self):
        key = GroupKey.congruence("gamma0", LevelIdeal.rational(2, 3))
        G = replace(get_group(2), key=key, inventory_status="incomplete",
                    inventory_proof_id=None, systole_trace=None,
                    trace_backend_id=SyntheticCongruenceBackend.backend_id)
        registry = GroupRegistry()
        registry.register(key, lambda: G)
        analytic = AnalyticRegistry()
        analytic.register(SyntheticCongruenceBackend())
        E = evaluate(key, include_elliptic=False, verbose=False,
                     group_registry=registry, analytic_registry=analytic)
        report = report_payload(E, analytic_registry=analytic)
        self.assertEqual(report["group_identity"], key.payload())
        self.assertIsNone(report["systole_trace"])
        self.assertFalse(report["spectral_certificate"])
        with self.assertRaises(ValueError):
            certificate_payload(E, group_registry=registry, analytic_registry=analytic)
        # A mismatched geometry cannot borrow the full group's systole.
        class WrongGeometry(SyntheticCongruenceBackend):
            def geometry(self, group):
                return Geometry(GroupKey(2), arb(1), arb(1))
        wrong = AnalyticRegistry()
        wrong.register(WrongGeometry())
        with self.assertRaisesRegex(ValueError, "different group identity"):
            evaluate(G, verbose=False, include_elliptic=False, analytic_registry=wrong)

    def test_v1_certificate_still_has_identical_numerical_endpoints(self):
        E = evaluate(2, verbose=False)
        frozen = json.loads((ROOT_PATH/"certificates/d2-k2.json").read_text())
        self.assertEqual(E.bound.lower().str(40), frozen["lower_endpoint_ball"])
        self.assertEqual(E.bound.upper().str(40), frozen["upper_endpoint_ball"])

    def test_v2_certificate_preserves_v1_math_and_records_backend_scope(self):
        E = evaluate(2, verbose=False)
        payload = certificate_payload(E)
        v1 = json.loads((ROOT_PATH/"certificates/d2-k2.json").read_text())
        v2 = json.loads((ROOT_PATH/"certificates/d2-k2-v2.json").read_text())
        replayed = json.loads(json.dumps(payload))
        # Python version is provenance, rather than a numerical requirement.
        replayed["arithmetic"]["python"] = v2["arithmetic"]["python"]
        self.assertEqual(replayed, v2)
        for key in ("inventory_proof", "parameters", "lower_endpoint_ball", "upper_endpoint_ball"):
            self.assertEqual(v1[key], v2[key])
        self.assertEqual(v2["schema_version"], 2)
        self.assertEqual(v2["group_identity"], GroupKey(2).payload())
        self.assertEqual(v2["inventory_backend"], "d2-arithmetic-v1")
        self.assertEqual(v2["analytic_backend"], "level-one-v1")

    def test_analytic_backend_rejects_incomplete_terms_and_nonfinite_geometry(self):
        key = GroupKey.congruence("gamma0", LevelIdeal.rational(2, 3))
        G = replace(get_group(2), key=key, inventory_status="incomplete",
                    trace_backend_id=SyntheticCongruenceBackend.backend_id)
        class InfiniteGeometry(SyntheticCongruenceBackend):
            def geometry(self, group):
                return Geometry(group.key, arb("inf"), arb(1))
        class MissingTerms(SyntheticCongruenceBackend):
            def terms(self, *args):
                return {"PHIINT": arb(0)}
        class FloatTerms(SyntheticCongruenceBackend):
            def terms(self, *args):
                return {name: 0.0 for name in ("CE", "Ch0", "PARg0", "PSI", "PHIINT")}
        for backend in [InfiniteGeometry(), MissingTerms(), FloatTerms()]:
            registry = AnalyticRegistry()
            registry.register(backend)
            with self.assertRaises(ValueError):
                evaluate(G, include_elliptic=False, verbose=False, analytic_registry=registry)


if __name__ == "__main__":
    unittest.main()
