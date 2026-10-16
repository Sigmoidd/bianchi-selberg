#!/usr/bin/env python3
import unittest, math, json
from pathlib import Path
import numpy as np
from level_one_quasimode_analytic import Quasimode,representatives
from ford_active_cosets import exact_prune_proof,all_vertical_splits,vertical_cut_roots,PRUNED
from fast_scalar_residual import ScalarResidual
from reflection_even_trial import real_reflection_coefficients
ROOT=Path(__file__).resolve().parent

class PrunedGeometry(unittest.TestCase):
    def test_exact_coset_elimination(self):
        self.assertEqual(len(exact_prune_proof()),9)
        self.assertEqual(PRUNED,{3,6,9,11,12})
    def test_exact_reflection_projected_fourier_basis(self):
        src,modes,avg,coeff=real_reflection_coefficients()
        self.assertEqual(avg[(1,0)],1)
        for a,b,n in modes:
            self.assertEqual(avg[(a,b)],avg[(a,-b)])
            self.assertEqual(avg[(-b,a)],-avg[(a,b)])
    def test_vertical_roots_are_transition_boundaries(self):
        for x,t in [(-.42,.11),(.3,.39),(0.,0.)]:
            for gamma in exact_prune_proof()[1:]:
                _,_,c,d=gamma
                for cutoff in (.6,.7):
                    for y in vertical_cut_roots(gamma,x,t,cutoff):
                        z=complex(x,t)
                        h=y/(abs(c*z+d)**2+abs(c)**2*y*y)
                        self.assertAlmostEqual(h,cutoff,places=11)
        for x,t in [(0.,0.),(-.25,.2),(.48,.45)]:
            cuts=all_vertical_splits(x,t)
            self.assertGreaterEqual(len(cuts),2)
            self.assertEqual(cuts,sorted(cuts))

class NumericRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.q=Quasimode(50,degree=70,dps=20)
        coeffs=json.loads((ROOT/'analytic_quasimode_coefficients_50.json').read_text())['orbits']
        cls.a=np.array([complex(float(o['coefficient_real']),float(o['coefficient_imag'])) for o in coeffs])
        cls.x=np.array([-.45,-.22,.09,.44]);cls.t=np.array([.16,.3,.07,.43]);cls.y=np.array([1.08,1.15,1.19,1.11])
    def test_nine_and_fourteen_cosets_are_identical(self):
        full=ScalarResidual(self.q,self.a)
        w0,r0=full.value_and_residual(self.x,self.t,self.y)
        self.q.gammas=exact_prune_proof()
        w1,r1=full.value_and_residual(self.x,self.t,self.y)
        self.assertLess(max(abs(w1-w0)),1e-17)
        self.assertLess(max(abs(r1-r0)),1e-17)
    def test_scalar_against_orbit_evaluator(self):
        self.q.gammas=exact_prune_proof()
        W,R,_=self.q.W_residual_cols(self.x,self.t,self.y)
        w,r=ScalarResidual(self.q,self.a).value_and_residual(self.x,self.t,self.y)
        self.assertLess(max(abs(W@self.a-w)),1e-17)
        self.assertLess(max(abs(R@self.a-r)),1e-17)
    def test_real_projection_agrees_with_real_part(self):
        self.q.gammas=exact_prune_proof()
        w0,r0=ScalarResidual(self.q,self.a).value_and_residual(self.x,self.t,self.y)
        w1,r1=ScalarResidual(self.q,self.a.real.astype(complex)).value_and_residual(self.x,self.t,self.y)
        self.assertLess(max(abs(w1-w0.real)),1e-17)
        self.assertLess(max(abs(r1-r0.real)),1e-17)
if __name__=='__main__':unittest.main(verbosity=2)
