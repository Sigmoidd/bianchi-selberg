#!/usr/bin/env python3
"""High-precision, non-interval, non-Chebyshev regression for interval jet engine."""
from __future__ import annotations
import sys,mpmath as mp
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from seed_interval_prototype import Seed
mp.mp.dps=45
S=Seed();r=mp.mpf('6.622119348586011');P=2*mp.pi

def f_grad(x,t,y):
 f=fx=ft=fy=mp.mpf(0)
 for a,b,n,civ in S.terms:
  # coefficients from already exactly rational decimal source; their intervals are tiny
  c=mp.mpf(str((float(civ.a)+float(civ.b))/2))
  z=P*mp.sqrt(n)*y; k=mp.re(mp.besselk(1j*r,z));
  kp=-mp.re((mp.besselk(1j*r-1,z)+mp.besselk(1j*r+1,z))/2)
  F=y*k;dF=k+z*kp
  aa=P*a*x;bb=P*b*x; cc=P*a*t;dd=P*b*t
  A=mp.cos(aa)*mp.cos(dd)-mp.cos(bb)*mp.cos(cc)
  Ax=P*(-a*mp.sin(aa)*mp.cos(dd)+b*mp.sin(bb)*mp.cos(cc))
  At=P*(-b*mp.cos(aa)*mp.sin(dd)+a*mp.cos(bb)*mp.sin(cc))
  f+=c*F*A;fx+=c*F*Ax;ft+=c*F*At;fy+=c*dF*A
 return f,(fx,ft,fy)

def trans(i,x,t,y):
 a,b,c,d=S.reps[i]
 a=mp.mpc(int(a.real),int(a.imag));b=mp.mpc(int(b.real),int(b.imag));c=mp.mpc(int(c.real),int(c.imag));d=mp.mpc(int(d.real),int(d.imag))
 z=x+1j*t;q=c*z+d;pp=a*z+b;D=abs(q)**2+abs(c)**2*y*y
 zz=(pp*mp.conj(q)+a*mp.conj(c)*y*y)/D;h=y/D
 Dx=2*mp.re(c*mp.conj(q));Dt=2*mp.re(1j*c*mp.conj(q));Dy=2*abs(c)**2*y
 nx=a*mp.conj(q)+pp*mp.conj(c);nt=1j*a*mp.conj(q)-1j*pp*mp.conj(c);ny=2*a*mp.conj(c)*y
 J=tuple((nx0-zz*dx)/D for nx0,dx in zip((nx,nt,ny),(Dx,Dt,Dy)))
 hg=(-h*Dx/D,-h*Dt/D,(1-h*Dy)/D)
 return mp.re(zz),mp.im(zz),h,J,hg

def compute(x,t,y):
 f,g=f_grad(x,t,y)
 G=mp.mpf(1);dd=mp.mpf(0);Ggrad=[mp.mpf(0) for _ in range(3)]
 deltaGrad=[mp.mpf(0) for _ in range(3)]
 terms=[]
 for i in range(5):
  xp,tp,h,jac,hj=trans(i,x,t,y)
  if h<=mp.mpf('0.6'):continue
  u=(h-mp.mpf('0.6'))*10
  if h>=mp.mpf('0.7'):w=1;wp=wpp=0
  else:
   w=u**3*(10-15*u+6*u*u)
   wp=300*u*u*(1-u)*(1-u)
   wpp=6000*u*(1-u)*(1-2*u)
  fg,grad=f_grad(xp,tp,h)
  fd=fg-f
  gradfd=[mp.re(jac[k])*grad[0]+mp.im(jac[k])*grad[1]+hj[k]*grad[2]-g[k] for k in range(3)]
  gw=[wp*hj[k] for k in range(3)];Lw=h*wp-h*h*wpp
  G+=w;dd+=w*fd
  for k in range(3):Ggrad[k]+=gw[k];deltaGrad[k]+=w*gradfd[k]+gw[k]*fd
  terms.append((Lw,gw,fd,gradfd))
 delta=dd/G
 gdelta=[(deltaGrad[k]-delta*Ggrad[k])/G for k in range(3)]
 resid=sum((Lw*(d-delta)-2*y*y*sum((gw[k]*(gd[k]-gdelta[k]) for k in range(3))) for Lw,gw,d,gd in terms))/G
 return delta,resid

if __name__=='__main__':
 a,b=compute(mp.mpf('.3'),mp.mpf('.1'),mp.mpf('1.15'))
 print('mpmath at center W-f',mp.nstr(a,35));print('mpmath at center residual',mp.nstr(b,35))
