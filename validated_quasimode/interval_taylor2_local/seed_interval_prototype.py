#!/usr/bin/env python3
"""Experimental interval prototype (mpmath.iv) for Picard mass-transfer gate.
Proves bounds only for explicitly evaluated boxes, not for global Ford wedge.
"""
from __future__ import annotations
import sys, csv
from decimal import Decimal, getcontext
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent/'inputs'))
from level_one_quasimode_analytic import representatives,modes_and_orbits
from mpmath import iv
import json, math
getcontext().prec=90
iv.dps=30
ROOT=Path(__file__).parent/'inputs'
TAB=ROOT/'radial'/'radial_coeffs_all.csv'
COEF=ROOT/'analytic_quasimode_coefficients_50_reflection_even.json'

def I(s):return iv.mpf(str(s))
def bound(v):
    a=float(v.a);b=float(v.b)
    return (a,b)
def ivwidth(v):
    a,b=bound(v);return b-a

def Cmul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def Cconj(a):return (a[0],-a[1])
def Cadd(a,b):return (a[0]+b[0],a[1]+b[1])

class Seed:
 def __init__(self):
  self.cf={};self.cd={}
  for l in TAB.read_text().splitlines():
   if l.startswith('#') or not l.strip():continue
   n,k,fl,fh,dl,dh=l.split(',');n,k=int(n),int(k)
   self.cf.setdefault(n,[]).append(I((Decimal(fl)+Decimal(fh))/2))
   self.cd.setdefault(n,[]).append(I((Decimal(dl)+Decimal(dh))/2))
  modes,orbits,S=modes_and_orbits(50)
  data=json.loads(COEF.read_text());amap={}
  for j,o in enumerate(data['orbits']):
   a=Decimal(o['coefficient_real'])
   for idx in orbits[j]:
    u,v,n=modes[idx];amap[(u,v)]=a*int(S[j,idx])
  self.terms=[]
  for a,b,n in modes:
   if a>b>=0 and amap[(a,b)]:
    c=(2 if b==0 else 4)*amap[(a,b)]
    self.terms.append((a,b,n,I(c)))
  assert len(self.terms)==21
  self.pi=iv.pi
  self.reps=[representatives()[i] for i in (1,2,4,5,13)]
  self.rerr=I('1e-20')
 
 def cheb(self,c,t):
  # Clenshaw, exact interval computation; no float polynomial evaluations.
  b1=b2=I(0)
  for ck in c[:0:-1]:
   b0=2*t*b1-b2+ck;b2,b1=b1,b0
  return c[0]+t*b1-b2
 def seed(self,x,t,y):
  tau=I(15)/8*(y-I(17)/15)
  vals={n:self.cheb(c,tau)+iv.mpf(['-1e-20','1e-20']) for n,c in self.cf.items()}
  P=2*self.pi
  cx={a:iv.cos(P*a*x) for a,b,n,c in self.terms}
  cx.update({b:iv.cos(P*b*x) for a,b,n,c in self.terms})
  ct={a:iv.cos(P*a*t) for a,b,n,c in self.terms}
  ct.update({b:iv.cos(P*b*t) for a,b,n,c in self.terms})
  s=I(0)
  for a,b,n,c in self.terms:
   s+=c*vals[n]*(cx[a]*ct[b]-cx[b]*ct[a])
  return s
 def transform(self,j,x,t,y):
  a,b,c,d=self.reps[j]
  ac=(I(int(a.real)),I(int(a.imag)));bc=(I(int(b.real)),I(int(b.imag)))
  cc=(I(int(c.real)),I(int(c.imag)));dc=(I(int(d.real)),I(int(d.imag)))
  z=(x,t);q=Cadd(Cmul(cc,z),dc);p=Cadd(Cmul(ac,z),bc)
  D=q[0]*q[0]+q[1]*q[1]+(cc[0]*cc[0]+cc[1]*cc[1])*y*y
  num=Cadd(Cmul(p,Cconj(q)),(Cmul(ac,Cconj(cc))[0]*y*y,Cmul(ac,Cconj(cc))[1]*y*y))
  return num[0]/D,num[1]/D,y/D
 def correct(self,x,t,y):
  seed=self.seed(x,t,y)
  n=I(0);G=I(1)
  for j in range(5):
   xx,tt,h=self.transform(j,x,t,y)
   lo,hi=bound(h)
   if hi<=.6:continue
   if lo>=.7:psi=I(1)
   elif lo>=.6 and hi<=.7:
    u=(h-I(3)/5)*10
    # fifth order monotonic, at any t in [0,1]
    psi=u**3*(I(10)-I(15)*u+I(6)*u*u)
   else:
    # Safe hull of monotone cutoff at lower and upper interval endpoints.
    def clipq(q):
     if q<=.6:return 0.
     if q>=.7:return 1.
     z=10*(q-.6);return z**3*(10-15*z+6*z*z)
    a=clipq(lo);b=clipq(hi)
    psi=iv.mpf([a-1e-14,b+1e-14])
   diff=self.seed(xx,tt,h)-seed
   n+=psi*diff;G+=psi
  return n/G

if __name__=='__main__':
 from time import perf_counter
 s=Seed();t0=perf_counter()
 for d in [0,0.001,0.01,0.05]:
  x=iv.mpf([.24999-d/2,.24999+d/2]);t=iv.mpf([.11999-d/2,.11999+d/2]);y=iv.mpf([1.1-d/2,1.1+d/2])
  q=s.correct(x,t,y);print('width',d,'bound correction',bound(q),'elapsed',perf_counter()-t0,flush=True)
