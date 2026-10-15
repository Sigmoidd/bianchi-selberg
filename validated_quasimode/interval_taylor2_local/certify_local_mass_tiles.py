#!/usr/bin/env python3
"""Strictly local, mpmath.iv mass-transfer proof primitive; no global completeness claim."""
from __future__ import annotations
from decimal import Decimal as D
from mpmath import iv
from seed_interval_prototype import I,bound
from taylor2_defect_iv import Taylor2Defect

def bump_on_h(h):
 """Monotone quintic cutoff, with rigorously outward interval endpoints."""
 c0=I(3)/5;c1=I(7)/10
 if h.b<=c0:return I(0)
 if h.a>=c1:return I(1)
 def at(x):
  if x<=c0:return I(0)
  if x>=c1:return I(1)
  z=10*(x-c0)
  return z*z*z*(I(10)-I(15)*z+I(6)*z*z)
 return iv.mpf([at(h.a).a,at(h.b).b])

def eval_tile(seed,lo,hi):
 boxes=[iv.mpf([str(a),str(b)]) for a,b in zip(lo,hi)]
 num=I(0);G=I(1);used=[]
 for j in range(5):
  u,v,h=seed.T(j,*[__import__('taylor2_defect_iv').Jet.variable(b,i) for i,b in enumerate(boxes)])
  # h.v encloses all heights in tile; intersection with y>=.6 required
  if h.v.b<=I(3)/5:continue
  if h.v.a<I(3)/5:
   # D enclosure may access the unvalidated Chebyshev region below .6;
   # so classify as undecided, rather than extrapolate a radial theorem.
   return {'status':'UNDECIDED_ACTIVE_CUTOFF_CROSSING','index':j,'height':bound(h.v)}
  w=bump_on_h(h.v)
  Dval,_,_=seed.taylor2(j,lo,hi)
  num+=w*Dval;G+=w
  used.append(j)
 q=num/G
 return {'status':'PASS' if max(abs(bound(q)[0]),abs(bound(q)[1]))<1e-6 else 'UNDECIDED_TOO_WIDE',
         'enclosure':[format(float(q.a),'.18e'),format(float(q.b),'.18e')],
         'active_indices':used}

if __name__=='__main__':
 from time import perf_counter
 s=Taylor2Defect();st=perf_counter()
 for d in [D('0.005'),D('0.01'),D('0.02')]:
  m=[D('.3'),D('.1'),D('1.15')];lo=[str(v-d/2) for v in m];hi=[str(v+d/2) for v in m]
  q=eval_tile(s,lo,hi)
  print('tile',d,q,'seconds',round(perf_counter()-st,2),flush=True)
