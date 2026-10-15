#!/usr/bin/env python3
"""Validated-arithmetic prototype: box-wise centered automorphy defects.
Uses shipped MPFR-verified Chebyshev radial coefficient intervals.
IMPORTANT: an individual PASS certifies only that box, not the whole Ford wedge.
"""
from __future__ import annotations
import math
from mpmath import iv
from seed_interval_prototype import Seed,I,bound,ivwidth,Cadd,Cmul,Cconj

class CenteredDefect(Seed):
 def radial(self,h):
  tau=I(15)/8*(h-I(17)/15)
  return ({n:self.cheb(c,tau)+iv.mpf(['-1e-20','1e-20']) for n,c in self.cf.items()},
          {n:self.cheb(c,tau)+iv.mpf(['-4e-20','4e-20']) for n,c in self.cd.items()})
 def seed_grad(self,x,t,h):
  vals,ders=self.radial(h)
  P=2*self.pi
  ux={i:iv.cos(P*i*x) for a,b,n,c in self.terms for i in (a,b)}
  ut={i:iv.cos(P*i*t) for a,b,n,c in self.terms for i in (a,b)}
  sx={i:iv.sin(P*i*x) for a,b,n,c in self.terms for i in (a,b)}
  st={i:iv.sin(P*i*t) for a,b,n,c in self.terms for i in (a,b)}
  f=fx=ft=fh=I(0)
  for a,b,n,c in self.terms:
   A=ux[a]*ut[b]-ux[b]*ut[a]
   Ax=P*(-a*sx[a]*ut[b]+b*sx[b]*ut[a])
   At=P*(-b*ux[a]*st[b]+a*ux[b]*st[a])
   f+=c*vals[n]*A;fx+=c*vals[n]*Ax;ft+=c*vals[n]*At;fh+=c*ders[n]*A
  return f,(fx,ft,fh)
 def transform_grad(self,j,x,t,y):
  a,b,c,d=self.reps[j]
  ac=(I(int(a.real)),I(int(a.imag)));bc=(I(int(b.real)),I(int(b.imag)))
  cc=(I(int(c.real)),I(int(c.imag)));dc=(I(int(d.real)),I(int(d.imag)))
  z=(x,t);q=Cadd(Cmul(cc,z),dc);p=Cadd(Cmul(ac,z),bc)
  D=q[0]*q[0]+q[1]*q[1]+(cc[0]*cc[0]+cc[1]*cc[1])*y*y
  acu=Cmul(ac,Cconj(cc));num=Cadd(Cmul(p,Cconj(q)),(acu[0]*y*y,acu[1]*y*y))
  zp=(num[0]/D,num[1]/D);h=y/D
  # complex gradient of projective numerator
  ic=(-cc[1],cc[0]); ia=(-ac[1],ac[0]); ip=(-p[1],p[0]);
  Dx=2*Cmul(cc,Cconj(q))[0];Dt=2*Cmul(ic,Cconj(q))[0];Dy=2*(cc[0]*cc[0]+cc[1]*cc[1])*y
  nx=Cadd(Cmul(ac,Cconj(q)),Cmul(p,Cconj(cc)))
  nt=Cadd(Cmul(ia,Cconj(q)),(-Cmul(ip,Cconj(cc))[0],-Cmul(ip,Cconj(cc))[1]))
  ny=(2*acu[0]*y,2*acu[1]*y)
  zs=tuple(((nj[0]-zp[0]*Dj)/D,(nj[1]-zp[1]*Dj)/D) for nj,Dj in zip((nx,nt,ny),(Dx,Dt,Dy)))
  hs=(-h*Dx/D,-h*Dt/D,(I(1)-h*Dy)/D)
  return zp[0],zp[1],h,zs,hs
 def diff_grad(self,j,x,t,y):
  s,g=self.seed_grad(x,t,y)
  xx,tt,h,zs,hs=self.transform_grad(j,x,t,y)
  sg,gg=self.seed_grad(xx,tt,h)
  gr=tuple(gg[0]*zs[k][0]+gg[1]*zs[k][1]+gg[2]*hs[k]-g[k] for k in range(3))
  return sg-s,gr
 def box_mean_value(self,j,lo,hi):
  mid=[(aa+bb)/2 for aa,bb in zip(lo,hi)]
  center,centergrad=self.diff_grad(j,*[I(v) for v in mid])
  xs=[iv.mpf([str(a),str(b)]) for a,b in zip(lo,hi)]
  _d,g=self.diff_grad(j,*xs)
  shift=I(0)
  for grad,a,b in zip(g,lo,hi):shift+=grad*iv.mpf([str(-(b-a)/2),str((b-a)/2)])
  return center+shift,center,centergrad,g
 def correction(self,lo,hi):
  mid=[(aa+bb)/2 for aa,bb in zip(lo,hi)]
  xx=[iv.mpf([str(a),str(b)]) for a,b in zip(lo,hi)]
  num=I(0);den=I(1);unresolved=0
  for j in range(len(self.reps)):
   _u,_v,h,*_=self.transform_grad(j,*xx)
   hl,hu=bound(h)
   if hu <=0.6:continue
   if hl>=0.7:w=I(1)
   elif hl>=0.6 and hu<=0.7:
    tl=10*(hl-0.6);tu=10*(hu-0.6)
    def bump(t):return t**3*(10-15*t+6*t*t)
    w=iv.mpf([max(0,bump(tl)-2e-12),min(1,bump(tu)+2e-12)])
   else:
    unresolved+=1
    def bump(t):return t**3*(10-15*t+6*t*t)
    loq=0 if hl<=0.6 else (1 if hl>=0.7 else bump(10*(hl-.6)))
    hiq=1 if hu>=0.7 else (0 if hu<=0.6 else bump(10*(hu-.6)))
    w=iv.mpf([max(0,loq-2e-12),min(1,hiq+2e-12)])
   D,_,_,_=self.box_mean_value(j,lo,hi)
   num+=w*D;den+=w
  return num/den,unresolved

if __name__=='__main__':
 from time import perf_counter
 obj=CenteredDefect();t0=perf_counter()
 for d in [0,0.005,0.01,0.02,0.04,0.08]:
  mid=(.3,.1,1.15);lo=[v-d/2 for v in mid];hi=[v+d/2 for v in mid]
  v,amb=obj.correction(lo,hi)
  print('diam=',d,'corr',bound(v),'ambiguous',amb,'elapsed',perf_counter()-t0,flush=True)
