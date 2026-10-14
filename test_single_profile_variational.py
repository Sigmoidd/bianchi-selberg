"""Independent local checks of the one-profile quotient trial and flux algebra."""
import json, math, unittest
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
import certify_single_profile_floor as c

ROOT=Path(__file__).resolve().parent

class SingleProfileTests(unittest.TestCase):
    def test_parameter_below_45_exact(self):
        r=F(6633,1000)
        self.assertEqual(1+r*r,F(44996689,1000000))
        self.assertLess(1+r*r,45)

    def test_analytic_error_envelopes(self):
        self.assertLess(F(160000)*(F(10,13)**159)+F(1,10**100),F(1,10**12))
        self.assertLess(F(10000)*(F(5,7)**111),F(1,10**12))

    def test_vertical_pairing_and_odd_symmetry(self):
        mp.mp.dps=60
        x=mp.mpf('0.13');t=mp.mpf('0.27')
        u=lambda x,t:mp.cos(2*mp.pi*x)-mp.cos(2*mp.pi*t)
        self.assertTrue(mp.almosteq(u(x+1,t),u(x,t)))
        self.assertTrue(mp.almosteq(u(-x,-t),u(x,t)))
        self.assertTrue(mp.almosteq(u(-x,t),u(x,t)))
        self.assertTrue(mp.almosteq(u(-t,x),-u(x,t)))
        self.assertTrue(mp.almosteq(u(-x,0),u(x,0)))
        self.assertTrue(mp.almosteq(u(-x,mp.mpf('0.5')),u(x,mp.mpf('0.5'))))

    def test_floor_flux_formula_with_direct_bessel_derivative(self):
        mp.mp.dps=70
        x=mp.mpf('0.21');t=mp.mpf('0.13');y=mp.sqrt(1-x*x-t*t)
        r=mp.mpf(6633)/1000
        z=2*mp.pi*y;v=1j*r
        K=mp.re(mp.besselk(v,z));Kp=-mp.re((mp.besselk(v-1,z)+mp.besselk(v+1,z))/2)
        U=mp.cos(2*mp.pi*x)-mp.cos(2*mp.pi*t)
        Ux=-2*mp.pi*mp.sin(2*mp.pi*x)
        Ut= 2*mp.pi*mp.sin(2*mp.pi*t)
        f=2*y*K*U
        fy=2*(K+y*2*mp.pi*Kp)*U
        fx=2*y*K*Ux
        ft=2*y*K*Ut
        direct= -f*(x*fx+t*ft+y*fy)/y**2
        formula=-4*(K*K*(U*U+U*(x*Ux+t*Ut))+z*K*Kp*U*U)
        self.assertLess(abs(direct-formula),mp.mpf('1e-65'))

    def test_replay_flags(self):
        for name in ('single_profile_floor_interval_result.json','single_profile_floor_interval_independent_result.json'):
            with self.subTest(name=name):
                d=json.loads((ROOT/name).read_text())
                self.assertTrue(d['sign_negative'])
                self.assertEqual(d['lambda0'],'44996689/1000000')
                self.assertLess(float(d['floor_flux_upper'].strip('[]').split(',')[0]),-3e-11)

if __name__=='__main__':unittest.main()
