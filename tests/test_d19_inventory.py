"""d=19 completeness gates, independent translation traces, and frozen replay."""
from dataclasses import replace
import json
from pathlib import Path
import unittest
from flint import arb, ctx
from core.assemble import evaluate
from core.certificate import certificate_payload
from groups import get_group
from groups.d19_inventory import RING, verify_inventory, verify_units

ROOT = Path(__file__).resolve().parents[1]


def trace_coefficient(C):
    """Use the matrix trace, without the stored radical norm or coefficient method."""
    T = C.primitive_translation
    t = RING.add(T[0], T[3])
    discriminant = RING.sub(RING.mul(t, t), (4, 0))
    length = ((arb(RING.norm(t)) + arb(RING.norm(discriminant)).sqrt())/4).acosh()
    sine_square = arb(1) if C.m == 2 else arb(3)/4
    return length/(4*C.finite_centralizer_order*sine_square)


class D19InventoryTests(unittest.TestCase):
    def setUp(self):
        ctx.prec = 100

    def test_complete_inventory_and_inverse_obstruction(self):
        proof = verify_inventory()
        self.assertEqual(proof['number_theory']['maximal_class_numbers'], {'A2': 1, 'A3': 1})
        self.assertEqual(proof['class_count'], 2)
        self.assertEqual(verify_units('A2')['determinant_image'], [-1, 1])
        self.assertEqual(verify_units('A3')['determinant_image'], [-1, 1])
        self.assertEqual([C['finite_centralizer_order'] for C in proof['classes']], [4, 3])
        self.assertTrue(proof['splitting']['order_three_inverse_classes_merge'])

    def test_record_changes_are_rejected_before_evaluation(self):
        G = get_group(19)
        for index, C in enumerate(G.elliptic_classes):
            for change in (dict(norm_denominator=2), dict(finite_centralizer_order=2*C.finite_centralizer_order),
                           dict(primitive_translation=C.representative), dict(cuspidal=True)):
                classes = list(G.elliptic_classes)
                classes[index] = replace(C, **change)
                with self.subTest(index=index, change=change), self.assertRaises(ArithmeticError):
                    evaluate(replace(G, elliptic_classes=tuple(classes)), verbose=False)
            with self.subTest(drop=index), self.assertRaises(ArithmeticError):
                evaluate(replace(G, elliptic_classes=tuple(c for i,c in enumerate(G.elliptic_classes)
                                                        if i != index)), verbose=False)
        C = G.elliptic_classes[0]
        with self.assertRaises(ArithmeticError):
            evaluate(replace(G, elliptic_classes=(replace(C, flip=None), *G.elliptic_classes[1:])), verbose=False)

    def test_coefficients_from_independent_matrix_trace_path(self):
        G = get_group(19)
        production = [C.coefficient() for C in G.elliptic_classes]
        ctx.prec = 200
        try:
            independent = [trace_coefficient(C) for C in G.elliptic_classes]
            for p, v in zip(production, independent):
                self.assertTrue(p.contains(v), (p, v))
            total = sum(independent, arb(0))
            closed = (170+39*arb(19).sqrt()).log()/8 + 2*(151+20*arb(57).sqrt()).log()/9
            self.assertTrue(total.overlaps(closed))
        finally:
            ctx.prec = 100

    def test_frozen_full_certificate_exact_endpoints_and_support(self):
        E = evaluate(19, frac=1., R=256, verbose=False)
        payload = certificate_payload(E)
        frozen = json.loads((ROOT/'certificates/d19-k2.json').read_text())
        for key in ('parameters', 'bound_ball', 'lower_endpoint_ball', 'upper_endpoint_ball'):
            self.assertEqual(payload[key], frozen[key])
        self.assertEqual(json.loads(json.dumps(payload['inventory_proof'])), frozen['inventory_proof'])
        self.assertTrue(E.bound.upper() < arb('0.9'))
        self.assertTrue(E.bound.lower() > arb('0.89'))
        self.assertTrue(E.terms['CE'].is_zero())
        with self.assertRaisesRegex(ValueError, 'support'):
            certificate_payload(replace(E, delta=E.delta*1.001))

    def test_search_claim_is_best_of_55_only(self):
        search = json.loads((ROOT/'certificates/d19-search.json').read_text())
        self.assertFalse(search['global_optimality_proved'])
        self.assertEqual(len(search['evaluations']), 55)
        winner = min(search['evaluations'], key=lambda row: row['upper'])
        for key in ('k', 'frac', 'R', 'delta_hex'):
            self.assertEqual(winner[key], search['selected'][key])
