#!/usr/bin/env python3
"""LOCAL validated interval commutator with centered second-order defect jets.
No global Ford enclosure claimed; a PASS only certifies a named rational cell.
"""
from __future__ import annotations
from decimal import Decimal as D
from mpmath import iv
from seed_interval_prototype import I,bound
from taylor2_defect_iv import Taylor2Defect,Jet,Z
from certify_local_mass_tiles import bump_on_h


def cutoff_derivatives(h):
 lo,hi=h.a,h.b
 if hi<=I(3)/5 or lo>=I(7)/10:return I(0),I(0)
 if lo<I(3)/5 or hi>I(7)/10:
  return iv.mpf([0,18.75]),iv.mpf([-600,600])
 z=10*(h-I(3)/5)
 wp=300*z*z*(1-z)*(1-z)
 wpp=6000*z*(1-z)*(1-2*z)
 return wp,wpp

def centered_gradient(seed,j,lo,hi):
 _,c,h=seed.taylor2(j,lo,hi)
 ds=[(D(str(b))-D(str(a)))/2 for a,b in zip(lo,hi)]
 dev=[iv.mpf([str(-d),str(d)]) for d in ds]
 g=tuple(c.g[i]+sum((h.h[i][k]*dev[k] for k in range(3)),I(0)) for i in range(3))
 return _,g

def residual_on_tile(seed,lo,hi,verbose=False):
 xs=[Jet.variable(iv.mpf([str(a),str(b)]),i) for i,(a,b) in enumerate(zip(lo,hi))]
 G=I(1);N=I(0);gs=[I(0) for _ in range(3)];nds=[I(0) for _ in range(3)]
 cache=[]
 for j in range(5):
  u,v,h=seed.T(j,*xs)
  if h.v.b<=I(3)/5:continue
  if h.v.a<I(3)/5:
   return dict(status='UNDECIDED_ACTIVE_CUTOFF_CROSSING',index=j,height=bound(h.v))
  w=bump_on_h(h.v);wp,wpp=cutoff_derivatives(h.v)
  H=h.v
  lw=H*wp-H*H*wpp
  dw=tuple(wp*q for q in h.g)
  dv,gd=centered_gradient(seed,j,lo,hi)
  G+=w;N+=w*dv
  for i in range(3):
   gs[i]+=dw[i]
   nds[i]+=w*gd[i]+dw[i]*dv
  cache.append((j,w,lw,dw,dv,gd))
 delta=N/G
 gdelta=tuple((nds[i]-delta*gs[i])/G for i in range(3))
 y=xs[2].v
 r=I(0)
 for j,w,lw,dw,dv,gd in cache:
  r+=lw*(dv-delta)-2*y*y*sum((dw[i]*(gd[i]-gdelta[i]) for i in range(3)),I(0))
 r/=G
 maxabs=max(abs(bound(r)[0]),abs(bound(r)[1]))
 out=dict(status='PASS_NARROW_GATE' if maxabs<3.4e-10 else 'UNDECIDED_TOO_WIDE',
          enclosure=[format(float(r.a),'.18e'),format(float(r.b),'.18e')],
          active_cosets=[z[0] for z in cache],gate='3.4e-10')
 if verbose:out['correction']=[str(delta.a),str(delta.b)]
 return out

if __name__=='__main__':
 from time import perf_counter
 s=Taylor2Defect();st=perf_counter()
 for d in (D('0'),D('0.00001'),D('0.0001'),D('0.0005'),D('0.001')):
  m=[D('.3'),D('.1'),D('1.15')]
  lo=[str(v-d/2) for v in m];hi=[str(v+d/2) for v in m]
  q=residual_on_tile(s,lo,hi)
  print('side',d,q,'seconds',round(perf_counter()-st,2),flush=True)
