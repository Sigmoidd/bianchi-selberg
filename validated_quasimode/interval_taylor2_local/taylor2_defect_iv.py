#!/usr/bin/env python3
"""Second-order, cancellation-preserving interval Taylor box engine.
Everything depends only on exact rational input coefficients, validated
radial interpolation enclosures and mpmath.iv outward-rounded operations.
This program is a LOCAL certificate primitive, not a global-cover proof.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from seed_interval_prototype import Seed,I,bound
from mpmath import iv
import math
from decimal import Decimal
Z=I(0);ONE=I(1)

def dot(a,b):return sum((p*q for p,q in zip(a,b)),Z)
class Jet:
 __slots__=('v','g','h')
 def __init__(self,v,g=None,h=None):
  self.v=v if hasattr(v,'a') else I(v)
  self.g=g or (Z,Z,Z)
  self.h=h or ((Z,Z,Z),(Z,Z,Z),(Z,Z,Z))
 @staticmethod
 def variable(v,j):return Jet(v,tuple(ONE if i==j else Z for i in range(3)))
 def __add__(a,b):
  if not isinstance(b,Jet):b=Jet(b)
  return Jet(a.v+b.v,tuple(a.g[i]+b.g[i] for i in range(3)),tuple(tuple(a.h[i][j]+b.h[i][j] for j in range(3)) for i in range(3)))
 __radd__=__add__
 def __neg__(a):return Jet(-a.v,tuple(-x for x in a.g),tuple(tuple(-x for x in row) for row in a.h))
 def __sub__(a,b):return a+-b if isinstance(b,Jet) else a+Jet(-b)
 def __rsub__(a,b):return (-a)+b
 def __mul__(a,b):
  if not isinstance(b,Jet):b=Jet(b)
  return Jet(a.v*b.v,tuple(a.g[i]*b.v+a.v*b.g[i] for i in range(3)),tuple(tuple(a.h[i][j]*b.v+a.v*b.h[i][j]+a.g[i]*b.g[j]+a.g[j]*b.g[i] for j in range(3)) for i in range(3)))
 __rmul__=__mul__
 def inverse(a):
  v=a.v;vv=v*v;vvv=vv*v
  return Jet(ONE/v,tuple(-a.g[i]/vv for i in range(3)),tuple(tuple(2*a.g[i]*a.g[j]/vvv-a.h[i][j]/vv for j in range(3)) for i in range(3)))
 def __truediv__(a,b):return a*b.inverse() if isinstance(b,Jet) else a*Jet(b).inverse()
 def __rtruediv__(a,b):return a.inverse()*b
 def cos(a):
  sn=iv.sin(a.v);cn=iv.cos(a.v)
  return Jet(cn,tuple(-sn*v for v in a.g),tuple(tuple(-cn*a.g[i]*a.g[j]-sn*a.h[i][j] for j in range(3)) for i in range(3)))

class Taylor2Defect(Seed):
 def __init__(self):
  super().__init__();self.two_pi=2*self.pi
 def radialjet(self,n,y):
  tau=I(15)/8*(y.v-I(17)/15)
  v=self.cheb(self.cf[n],tau)+iv.mpf(['-1e-20','1e-20'])
  dy=self.cheb(self.cd[n],tau)+iv.mpf(['-4e-20','4e-20'])
  lam=ONE+I('6.622119348586011')**2
  k2=self.two_pi*self.two_pi*n
  ddy=k2*v+dy/y.v-lam*v/(y.v*y.v)
  return Jet(v,tuple(dy*y.g[i] for i in range(3)),tuple(tuple(ddy*y.g[i]*y.g[j]+dy*y.h[i][j] for j in range(3)) for i in range(3)))
 def f(self,x,t,y):
  vals={n:self.radialjet(n,y) for a,b,n,c in self.terms}
  cos_x={a:(x*(a*self.two_pi)).cos() for a,b,n,c in self.terms}
  cos_x.update({b:(x*(b*self.two_pi)).cos() for a,b,n,c in self.terms})
  cos_t={a:(t*(a*self.two_pi)).cos() for a,b,n,c in self.terms}
  cos_t.update({b:(t*(b*self.two_pi)).cos() for a,b,n,c in self.terms})
  s=Jet(Z)
  for a,b,n,c in self.terms:
   s+=(vals[n]*c)*(cos_x[a]*cos_t[b]-cos_x[b]*cos_t[a])
  return s
 def T(self,j,x,t,y):
  aa,bb,cc,dd=self.reps[j]
  a=(int(aa.real),int(aa.imag));b=(int(bb.real),int(bb.imag));c=(int(cc.real),int(cc.imag));d=(int(dd.real),int(dd.imag))
  def cm(p,q):return (p[0]*q[0]-p[1]*q[1],p[0]*q[1]+p[1]*q[0])
  def conj(p):return (p[0],-p[1])
  def ca(p,q):return (p[0]+q[0],p[1]+q[1])
  z=(x,t);q=ca(cm(c,z),d);p=ca(cm(a,z),b)
  den=q[0]*q[0]+q[1]*q[1]+(c[0]*c[0]+c[1]*c[1])*y*y
  ac=cm(a,conj(c))
  nm=ca(cm(p,conj(q)),(ac[0]*y*y,ac[1]*y*y))
  return nm[0]/den,nm[1]/den,y/den
 def D(self,j,x,t,y):
  u,v,h=self.T(j,x,t,y)
  return self.f(u,v,h)-self.f(x,t,y)
 def jet(self,j,lower,upper):
  xs=[Jet.variable(iv.mpf([str(a),str(b)]),i) for i,(a,b) in enumerate(zip(lower,upper))]
  return self.D(j,*xs)
 def taylor2(self,j,lower,upper):
  mids=[I((Decimal(str(a))+Decimal(str(b)))/2) for a,b in zip(lower,upper)]
  diam=[(Decimal(str(b))-Decimal(str(a)))/2 for a,b in zip(lower,upper)]
  center=self.D(j,*[Jet.variable(v,i) for i,v in enumerate(mids)])
  full=self.jet(j,lower,upper)
  dev=[iv.mpf([str(-d),str(d)]) for d in diam]
  result=center.v+dot(center.g,dev)
  for i in range(3):
   for j in range(3):result+=full.h[i][j]*dev[i]*dev[j]/2
  return result,center,full

if __name__=='__main__':
 from time import perf_counter
 s=Taylor2Defect();ts=perf_counter();mid=(.3,.1,1.15)
 for d in [0,.005,.01,.02,.04]:
  low=[v-d/2 for v in mid];hi=[v+d/2 for v in mid]
  r,c,h=s.taylor2(2,low,hi)
  print('side',d,'mean-value2',bound(r),'center',bound(c.v),'grad_center',[bound(a) for a in c.g],'Hmax',max(max(abs(float(z.a)),abs(float(z.b))) for row in h.h for z in row),'s',perf_counter()-ts,flush=True)
