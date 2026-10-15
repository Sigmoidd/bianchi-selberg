import unittest
from fractions import Fraction as F
from boundary_weak_residual import (build_basis,radial_tables,numerical_boundary,exact_conditions,direct_gauss,cheb_coefficients)

class TestBoundary(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.terms=build_basis();cls.pol,cls.der=radial_tables()
 def test_basis(self):
  self.assertEqual(len(self.terms),21)
  self.assertTrue(sum(abs(c) for *_,c in self.terms)<37)
 def test_exact_trace_and_tail(self):
  self.assertIn('conditional_L2_gate',exact_conditions())
 def test_boundary_convergence(self):
  a=direct_gauss(32,28,self.terms,self.pol,self.der)
  b=direct_gauss(90,70,self.terms,self.pol,self.der)
  self.assertLess(abs(a['J_L2_numerical']-b['J_L2_numerical']),2e-19)
  self.assertLess(b['J_L2_numerical'],1e-12)
 def test_cheb_diagnostic(self):
  row,_=cheb_coefficients(96,self.terms,self.pol,self.der)
  self.assertLess(row['abs_coefficient_sum_numerical'],3.5e-12)
 def test_flux_sign_and_magnitude(self):
  r=direct_gauss(60,50,self.terms,self.pol,self.der)
  self.assertLess(r['one_sided_floor_flux_numerical'],0)
  self.assertLess(abs(r['one_sided_floor_flux_numerical']+4.0193e-17),1e-21)

if __name__=='__main__':unittest.main()
