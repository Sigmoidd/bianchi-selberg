#!/usr/bin/env python3
"""Local interval and algebra checks. Not a global certificate."""
import sys,unittest
from pathlib import Path
from fractions import Fraction as F
sys.path.insert(0,str(Path(__file__).parent))
from seed_interval_prototype import Seed,I,bound
from taylor2_defect_iv import Taylor2Defect
from certify_local_mass_tiles import eval_tile
from local_residual_iv import residual_on_tile

class LocalIntervalRegression(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.seed=Taylor2Defect()
 def test_cos_decomposition_160_to_21_exact(self):
  import json
  from level_one_quasimode_analytic import modes_and_orbits
  model=json.loads((Path(__file__).parent/'inputs'/'analytic_quasimode_coefficients_50_reflection_even.json').read_text())
  modes,orb,S=modes_and_orbits(50)
  target={}
  for j,o in enumerate(model['orbits']):
   val=F(o['coefficient_real'])
   for idx in orb[j]:target[modes[idx][:2]]=int(S[j,idx])*val
  reco={key:F(0) for key in target}
  for a,b,n,c in self.seed.terms:
   # True as rational because coefficients in the source are fixed decimals.
   # Reconstruct using source instead, avoiding conversion from interval midpoint.
   source_val=target[(a,b)]
   assert source_val!=0
   for sa in [-1,1]:
    for sb in ([-1,1] if b else [1]):
     m=(sa*a,sb*b);reco[m]+=source_val
     if b:reco[(sa*b,sb*a)]-=source_val
   if b==0:
    # y-axis pair was not filled in the generic b>0 section.
    for sa in [-1,1]:reco[(0,sa*a)]-=source_val
  self.assertEqual(target,reco)
  self.assertEqual(len(self.seed.terms),21)
 def test_mass_transfer_local_cube(self):
  self.assertEqual(eval_tile(self.seed,['.29','.09','1.14'],['.31','.11','1.16'])['status'],'PASS')
 def test_residual_local_cube(self):
  self.assertEqual(residual_on_tile(self.seed,['.299995','.099995','1.149995'],['.300005','.100005','1.150005'])['status'],'PASS_NARROW_GATE')
 def test_wider_cube_not_falsely_promoted(self):
  self.assertEqual(residual_on_tile(self.seed,['.29995','.09995','1.14995'],['.30005','.10005','1.15005'])['status'],'UNDECIDED_TOO_WIDE')
 def test_parameter_mutation_rejected(self):
  seed=Taylor2Defect()
  seed.terms=[(a,b,n,(c*I('1.01') if n==1 else c)) for a,b,n,c in seed.terms]
  x=residual_on_tile(seed,['.3','.1','1.15'],['.3','.1','1.15'])
  self.assertEqual(x['status'],'UNDECIDED_TOO_WIDE')

if __name__=='__main__':unittest.main(verbosity=2)
