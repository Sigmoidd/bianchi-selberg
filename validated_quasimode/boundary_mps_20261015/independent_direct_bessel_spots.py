"""Independent NON-INTERVAL spot check against mpmath direct Bessel K."""
import mpmath as mp
import numpy as np
from boundary_weak_residual import build_basis,radial_tables,numerical_boundary
mp.mp.dps=48
r=mp.mpf('6.622119348586011')
terms=build_basis();pol,der=radial_tables()
for ss,tt in [('0.1','0.12'),('0.3','0.1'),('0.49','0.48')]:
 x,t=mp.mpf(ss),mp.mpf(tt);y=mp.sqrt(1-x*x-t*t)
 direct=mp.mpf(0)
 for a,b,n,c in terms:
  C=mp.mpf(c.numerator)/mp.mpf(c.denominator)
  z=2*mp.pi*mp.sqrt(n)*y; K=mp.re(mp.besselk(1j*r,z))
  Kp=mp.re(-(mp.besselk(1j*r-1,z)+mp.besselk(1j*r+1,z))/2)
  F=y*K;Fy=K+z*Kp
  ax=2*mp.pi*a*x;bx=2*mp.pi*b*x;at=2*mp.pi*a*t;bt=2*mp.pi*b*t
  T=mp.cos(ax)*mp.cos(bt)-mp.cos(bx)*mp.cos(at)
  Tx=2*mp.pi*(-a*mp.sin(ax)*mp.cos(bt)+b*mp.sin(bx)*mp.cos(at))
  Tt=2*mp.pi*(-b*mp.cos(ax)*mp.sin(bt)+a*mp.cos(bx)*mp.sin(at))
  direct-=y*C*(F*(x*Tx+t*Tt)+y*Fy*T)
 _,j,_=numerical_boundary(np.array([[float(x)]]),np.array([[float(t)]]),terms,pol,der)
 print(ss,tt,'direct=',mp.nstr(direct,22),'poly=',format(float(j[0,0]),'.17e'),'abs_diff=',mp.nstr(abs(direct-mp.mpf(float(j[0,0]))),12))
