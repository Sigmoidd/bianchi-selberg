"""Fast algebra/geometry regressions; none of these certify an L2 spectral bound."""
import math
import unittest
import numpy as np
from level_one_quasimode_analytic import Quasimode,cutoff,modes_and_orbits,representatives

class AlgebraAndGeometryTests(unittest.TestCase):
    def test_14_exact_candidate_representatives(self):
        reps=representatives()
        self.assertEqual(len(reps),14)
        self.assertEqual(len(set(reps)),14)
        for a,b,c,d in reps:
            self.assertEqual(a*d-b*c,1)
        # All excluded Gaussian denominators have height below the cutoff.
        self.assertLess(math.sqrt(2)/3, .6)  # |c|^2 >=3
        self.assertLess(1/(2*(2-1/math.sqrt(2))),.6) # c unit, |d|>=2
        self.assertLess(1/(2*math.sqrt(2)*(math.sqrt(5)-1)),.6) # N(c)=2, N(d)>=5

    def test_symmetry_is_exact_in_mode_basis(self):
        modes,orbits,S=modes_and_orbits(50)
        self.assertEqual(len(orbits),40)
        self.assertEqual(len(modes),160)
        index={(a,b):i for i,(a,b,_) in enumerate(modes)}
        for row in S:
            for i,(a,b,_) in enumerate(modes):
                self.assertEqual(row[index[(-b,a)]],-row[i])
                self.assertEqual(row[index[(-a,-b)]],row[i])

    def test_cutoff_is_c2(self):
        y=np.array([.58,.6,.65,.7,.72])
        v,v1,v2=cutoff(y,.6,.7)
        self.assertEqual(v[0],0)
        self.assertEqual(v[-1],1)
        for j in [0,1,3,4]:
            self.assertEqual(v1[j],0)
            self.assertEqual(v2[j],0)
        self.assertAlmostEqual(v[2],.5)

    def test_hyperbolic_height_identity(self):
        x=np.array([.17,-.31]);t=np.array([.19,.41]);y=np.array([.9,1.18])
        for gamma in representatives():
            xp,tp,h,zj,gradh=Quasimode.transform(gamma,x,t,y)
            self.assertTrue(np.allclose((y*y*sum(v*v for v in gradh))/(h*h),1,rtol=1e-12,atol=1e-12))

if __name__=='__main__':unittest.main()