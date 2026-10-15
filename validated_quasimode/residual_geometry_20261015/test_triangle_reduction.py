#!/usr/bin/env python3
"""Exact geometry plus numerical (not proof) symmetries for projected trial."""
import unittest,json
from pathlib import Path
import numpy as np
from certify_triangle_symmetry import TRIANGLE_INDICES,proof
from level_one_quasimode_analytic import Quasimode,representatives
from fast_scalar_residual import ScalarResidual
from ford_active_cosets import exact_prune_proof

DIR=Path(__file__).parent
class ExactGeometry(unittest.TestCase):
    def test_exact_pruning(self):
        self.assertEqual(len(proof()),14)
        self.assertEqual(TRIANGLE_INDICES,(0,1,2,4,5,13))

class Regression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.q=Quasimode(M=50,degree=45,dps=18)
        o=json.loads((DIR/'analytic_quasimode_coefficients_50_reflection_even.json').read_text())['orbits']
        cls.a=np.array([complex(float(row['coefficient_real']),float(row['coefficient_imag'])) for row in o])
        cls.s=ScalarResidual(cls.q,cls.a)
    def evaluate(self,coords,reps=None):
        self.q.gammas=reps or representatives()
        v=np.array(coords,dtype=float)
        return self.s.value_and_residual(v[:,0],v[:,1],v[:,2])
    def test_triangle_six_cosets_match_full(self):
        pts=[(.25,.03,1.14),(.41,.31,.92),(.45,.43,.81),(.46,.13,1.1),(.34,.34,1.5)]
        a,b=self.evaluate(pts)
        c,d=self.evaluate(pts,[representatives()[i] for i in TRIANGLE_INDICES])
        self.assertLess(max(abs(a-c)),1e-17)
        self.assertLess(max(abs(b-d)),1e-17)
    def test_symmetries(self):
        pts=[(.11,.41,1.22),(.42,.15,1.06),(.3,.3,1.42)]
        w,r=self.evaluate(pts)
        wx,rx=self.evaluate([(-x,t,y) for x,t,y in pts])
        ws,rs=self.evaluate([(t,x,y) for x,t,y in pts])
        self.assertLess(max(abs(w-wx)),1e-16)
        self.assertLess(max(abs(r-rx)),1e-16)
        self.assertLess(max(abs(w+ws)),1e-16)
        self.assertLess(max(abs(r+rs)),1e-16)
if __name__=='__main__':unittest.main(verbosity=2)
