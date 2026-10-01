from dataclasses import replace
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from flint import arb, ctx
from core.assemble import evaluate
from core.certificate import certificate_payload
from groups import get_group
from groups.d2_inventory import Radical, verify_inventory, verify_unit_reduction

ROOT = Path(__file__).resolve().parents[1]


class D2InventoryTests(unittest.TestCase):
    def setUp(self):
        ctx.prec = 100

    def test_complete_arithmetic_replay(self):
        proof = verify_inventory()
        self.assertEqual(proof["class_count"], 4)
        self.assertEqual([c["m"] for c in proof["classes"]], [2, 2, 3, 3])
        self.assertEqual(proof["lattices"]["trace_zero_GL_types"], ["A2", "S2"])
        self.assertEqual(proof["number_theory"]["discriminants"],
                         {"A2": 1024, "S2": 256, "A3": 576})

    def test_determinant_obstruction_is_not_discarded(self):
        self.assertEqual(verify_unit_reduction("A3")["determinant_image"], [1])
        for name in ["A2", "S2"]:
            self.assertEqual(verify_unit_reduction(name)["determinant_image"], [-1, 1])
        G = get_group(2)
        self.assertNotEqual(G.elliptic_classes[2].representative,
                            G.elliptic_classes[3].representative)

    def test_radical_comparisons_on_both_sides_of_zero(self):
        # Pell approximants exercise cancellation without floating point.
        self.assertEqual(Radical(665857, -470832, 2).sign(), 1)
        self.assertEqual(Radical(-665857, 470832, 2).sign(), -1)
        self.assertEqual(Radical(19601, -13860, 2).sign(), 1)
        self.assertEqual(Radical(5, -2, 6).sign(), 1)
        self.assertEqual(Radical(0, 0, 2).sign(), 0)

    def test_inventory_gate_rejects_missing_or_altered_classes(self):
        G = get_group(2)
        variants = [replace(G, elliptic_classes=G.elliptic_classes[:-1])]
        C = G.elliptic_classes[0]
        for change in [dict(finite_centralizer_order=2), dict(norm=(3, 2, 2)),
                       dict(primitive_translation=C.representative),
                       dict(flip=None), dict(normalization_status="historical")]:
            variants.append(replace(G, elliptic_classes=(replace(C, **change),)+G.elliptic_classes[1:]))
        variants.append(replace(G, ce_integral=1))
        for variant in variants:
            with self.assertRaises(ArithmeticError):
                evaluate(variant, verbose=False)

    def test_status_string_cannot_certify_another_field(self):
        forged = replace(get_group(7), inventory_status="self-contained")
        with self.assertRaisesRegex(ValueError, "no self-contained inventory verifier"):
            evaluate(forged, verbose=False)

    def test_full_coefficient_and_frozen_certificate(self):
        E = evaluate(2, verbose=False)
        payload = certificate_payload(E)
        coefficient = (arb(3)/16*(3+2*arb(2).sqrt()).log()
                       +arb(2)/9*(5+2*arb(6).sqrt()).log())
        self.assertTrue(E.group.analytic_data()["C_ell"].overlaps(coefficient))
        self.assertTrue(E.bound.lower() > arb("0.424"))
        self.assertTrue(E.bound.upper() < arb("0.431"))
        self.assertTrue(payload["spectral_certificate"])
        self.assertTrue(E.terms["CE"].is_zero())
        frozen = json.loads((ROOT/"certificates/d2-k2.json").read_text())
        self.assertEqual(json.loads(json.dumps(payload["inventory_proof"])), frozen["inventory_proof"])
        self.assertEqual(payload["parameters"], frozen["parameters"])
        self.assertTrue(E.bound.overlaps(arb(frozen["bound_ball"])))
        self.assertTrue(arb(frozen["upper_endpoint_ball"]) < 1)

    def test_export_rechecks_inventory_and_support(self):
        E = evaluate(2, verbose=False)
        incomplete = replace(E.group, elliptic_classes=E.group.elliptic_classes[:-1])
        with self.assertRaises(ArithmeticError):
            certificate_payload(replace(E, group=incomplete))
        with self.assertRaisesRegex(ValueError, "support"):
            certificate_payload(replace(E, delta=1.0))
        with self.assertRaisesRegex(ValueError, "B < 1"):
            certificate_payload(replace(E, bound=arb(2)))

    def test_certificate_cli_from_unrelated_directory(self):
        with tempfile.TemporaryDirectory() as cwd:
            output = Path(cwd)/"certificate.json"
            subprocess.run([sys.executable, str(ROOT/"examples/d2_certificate.py"),
                            "--output", str(output)], cwd=cwd, check=True,
                           capture_output=True, text=True)
            payload = json.loads(output.read_text())
            self.assertTrue(payload["spectral_certificate"])
            self.assertTrue(arb(payload["upper_endpoint_ball"]) < 1)


if __name__ == "__main__":
    unittest.main()
