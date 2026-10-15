#!/usr/bin/env python3
"""Exact rational constants supporting MPFR Bessel/Chebyshev radial enclosures."""
from fractions import Fraction as F
from math import factorial
rho_u=F(13,10);rho_y=F(3,2)
uymax=3*(rho_u-rho_u**-1)/2
assert uymax==F(207,260)<F(4,5)
assert 1-F(4,5)**2/2 >F(2,3)  # cos(Im u)>2/3
assert F(18,5)*F(2,3)==F(12,5)
x=F(12,5)
assert sum(x**j/factorial(j) for j in range(6))>10
assert 3*64*73/(15*(rho_u**2-1))<3000
assert F(3000)*(rho_u**-264)<F(1,10**25)
assert 2**720>10**200
assert sum(F(6)**k/factorial(k) for k in range(17))>401
# Gauss < 1e-25, tail <1e-200, hence <<1e-24.
assert F(3000)*rho_u**-264+F(1,10**200)<F(1,10**24)
# radial analytic ellipse
lo,hi=F(3,5),F(5,3)
center=(lo+hi)/2;half=(hi-lo)/2
assert center-half*(rho_y+1/rho_y)/2==F(5,9)
assert center+half*(rho_y+1/rho_y)/2==F(77,45)<2
assert F(2)*F(3)*F(5,9)==F(10,3)
assert sum(F(10,3)**j/factorial(j) for j in range(9))>27
# pi<22/7, sqrt 50<8, |y|<2 hence |z|<101
assert 2*F(22,7)*8*2<101
assert F(1+101,27)<4
n=120
err_f=F(4)*F(1,10)/(rho_y-1)*rho_y**-n
err_d=F(4)*4/(rho_y-1)*rho_y**-n
assert err_f<F(1,10**20) and err_d<F(3,10**20)
print('PASS EXACT: Bessel complex-u Gauss + truncation error < 1e-24 for m=132, z>=18/5')
print('PASS EXACT: analytic radial value envelope <1/10, derivative envelope <4 on E_3/2')
print('PASS EXACT: degree120 Cheb-Lobatto ideal interpolation errors: f<1e-20, fprime<3e-20')
print('Value ideal interp upper',float(err_f),'derivative',float(err_d))
