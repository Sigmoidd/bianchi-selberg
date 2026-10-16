#!/usr/bin/env python3
"""Exact-rational proof-gate audit. DOES NOT certify the nodal DCT."""
from fractions import Fraction as F
from math import factorial
from boundary_weak_residual import exact_conditions,build_basis

r=exact_conditions()
assert len(build_basis())==21
# Exact power-series bounds on exponentials used in the polyellipse comparison.
s6=sum(F(39,10)**k/factorial(k) for k in range(7))
assert s6>40, 'e^(39/10) > 40'
S100=sum(F(11)**k/factorial(k) for k in range(101))
term101=F(11)**101/factorial(101)
E11=S100+term101/(1-F(11,102))
assert E11<60000, 'e^11 < 60000'
# Complex ellipse x=(1/2)z,t=(1+z')/4, rho=3/2.
# |x| <= 13/24, |t| <=25/48; Im bounds 5/24, 5/48.
U=F(13,24)**2+F(25,48)**2
assert 1-U>F(13,20)**2
assert 1+U<F(13,10)**2
assert F(44,7)*F(5,3)<11
# K_ir and its derivative via real-u integral representation satisfy |K|,|K'|<1/40
# on the polyellipse: Re(2*pi*sqrt(n)*y)>39/10 and exp(39/10)>40.
# 2*pi*sqrt(n) <51; |F|<13/400, |F'|<17/10.
M=F(13,10)*37*(F(11,10)*F(13,400)*3780000+F(13,10)*F(17,10)*120000)
assert M<2*10**7
# 2-D degree 127 Chebyshev-Gauss interpolation of J, alias error
# <= 16*M*(3/2)^(-128)/(1-2/3)^2 = 144*M*(2/3)^128.
err=144*F(2*10**7)*F(2,3)**128
assert err<F(1,10**13)
# If outward-certified coeff absolute sum <3.5e-12, uniform J<3.6e-12<4e-12.
assert F(35,10**13)+F(1,10**13)<F(4,10**12)
# Trace constant 81/4 and area <=5/6; J∞<4e-12 and mu>13e-6
# imply spectral dist < 1e-5 by weak residual inclusion using sqrt(1+lambda+eps)<7.
assert 7*F(9,2)*F(4,10**12)<F(13,10**6)*F(1,10**5)
print('PASS all exact-rational conditional inclusion and polyellipse inequalities')
print('polyellipse M rational bound:', M, '(<20,000,000)')
print('Chebyshev degree-127 tail bound:', f'{float(err):.14g}')
print('e^11 rational upper bound:', f'{float(E11):.14g}')
print('PENDING: certify 128x128 Chebyshev coefficients with MPFR directed rounding')
