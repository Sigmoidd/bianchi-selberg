from dataclasses import replace
import json
from pathlib import Path
import unittest
from flint import arb, ctx
from core.assemble import evaluate
from core.certificate import certificate_payload
from groups import get_group
from groups.d7_inventory import verify_inventory, verify_units
from groups.exact import Radical

ROOT = Path(__file__).resolve().parents[1]


class D7InventoryTests(unittest.TestCase):
    def setUp(self):
        ctx.prec = 100

    def test_complete_inventory_and_determinant_splitting(self):
        proof = verify_inventory()
        self.assertEqual(proof['class_count'], 2)
        self.assertEqual(proof['number_theory']['maximal_class_numbers'], {'A2': 1, 'A3': 1})
        self.assertEqual(verify_units('A2')['determinant_image'], [1])
        self.assertEqual(verify_units('A3')['determinant_image'], [-1, 1])
        self.assertEqual([C['finite_centralizer_order'] for C in proof['classes']], [2, 3])
        self.assertEqual(proof['splitting']['trace_zero_GL_to_SL'], 2)
        self.assertEqual(proof['splitting']['trace_zero_PSL_classes'], 1)
        self.assertEqual(proof['splitting']['trace_one_PSL_classes'], 1)

    def test_rational_radical_square_and_reciprocal(self):
        r = Radical(5, 1, 21, 2)
        self.assertEqual(r.square(), Radical(23, 5, 21, 2))
        self.assertEqual(r.unit_reciprocal(), Radical(5, -1, 21, 2))
        self.assertGreater(r.compare(1), 0)
        self.assertLess(r.unit_reciprocal().compare(1), 0)

    def test_record_tampering_cannot_lower_bound(self):
        G = get_group(7)
        C = G.elliptic_classes[1]
        for change in (dict(norm_denominator=4), dict(finite_centralizer_order=6),
                       dict(primitive_translation=C.representative), dict(cuspidal=True)):
            with self.assertRaises(ArithmeticError):
                evaluate(replace(G, elliptic_classes=(G.elliptic_classes[0], replace(C, **change))), verbose=False)
        with self.assertRaises(ArithmeticError):
            evaluate(replace(G, elliptic_classes=G.elliptic_classes[:1]), verbose=False)

    def test_full_selected_certificate_replays_and_support_is_checked(self):
        E = evaluate(7, frac=1., R=256, verbose=False)
        payload = certificate_payload(E)
        frozen = json.loads((ROOT/'certificates/d7-k2.json').read_text())
        self.assertTrue(E.bound.upper() < arb('0.3641'))
        self.assertTrue(E.bound.lower() > arb('0.3639'))
        self.assertTrue(E.terms['CE'].is_zero())
        coefficient = ((8+3*arb(7).sqrt()).log()/8
                       +((23+5*arb(21).sqrt())/2).log()/9)
        self.assertTrue(coefficient.overlaps(arb(payload['elliptic_coefficient_ball'])))
        self.assertEqual(payload['parameters'], frozen['parameters'])
        self.assertEqual(json.loads(json.dumps(payload['inventory_proof'])), frozen['inventory_proof'])
        self.assertTrue(E.bound.overlaps(arb(frozen['bound_ball'])))
        with self.assertRaisesRegex(ValueError, 'support'):
            certificate_payload(replace(E, delta=E.delta*1.001))

    def test_search_manifest_has_only_finite_optimality_claim(self):
        search = json.loads((ROOT/'certificates/d7-search.json').read_text())
        self.assertFalse(search['global_optimality_proved'])
        winner = min(search['evaluations'], key=lambda row: row['upper'])
        for key in ('k', 'frac', 'R', 'delta_hex'):
            self.assertEqual(winner[key], search['selected'][key])
        self.assertEqual(len(search['evaluations']), 55)
