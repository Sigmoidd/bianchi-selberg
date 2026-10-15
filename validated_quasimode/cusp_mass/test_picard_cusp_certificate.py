import hashlib
import json
import subprocess
import unittest
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parent
EXPECTED='b04f90ca943e5d6d574f5f7c47d7649c8fdadd6b40607f6cc9064a46cdbb671b'

class CuspMassCertificateTests(unittest.TestCase):
    def test_exact_frozen_source(self):
        coeff=ROOT/'analytic_quasimode_coefficients_50.json'
        if not coeff.exists(): coeff=ROOT.parent.parent/'analytic_quasimode_coefficients_50.json'
        raw=coeff.read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(),EXPECTED)
        data=json.loads(raw)
        self.assertEqual(Fraction(data['r']),Fraction(6622119348586011,10**15))
        self.assertEqual(data['M'],50)
        self.assertEqual(len(data['orbits']),40)
        orb=data['orbits'][0]
        self.assertEqual(Fraction(orb['coefficient_real']),1)
        self.assertEqual(Fraction(orb['coefficient_imag']),0)
        self.assertEqual(orb['modes'],[[-1,0],[0,-1],[1,0],[0,1]])

    def test_exact_quadrature_error_audit(self):
        p=subprocess.run(['python','audit_cusp_simpson_bounds.py'],cwd=ROOT,capture_output=True,text=True,check=True)
        self.assertIn('1e-10',p.stdout)

    def test_mpfr_coarse_and_fine(self):
        for f in ['cusp_mass_mpfr.log','cusp_mass_mpfr_refined.log']:
            s=(ROOT/f).read_text()
            self.assertIn('PASS: mu_tail > 4e-7',s)
            self.assertIn('certified tail mass lower squared [1.912747',s)

if __name__=='__main__':unittest.main()
