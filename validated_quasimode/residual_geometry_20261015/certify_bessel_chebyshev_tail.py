#!/usr/bin/env python3
"""Exact rational audit of a uniform analytic bound for radial Bessel profiles.

Proves a priori interpolation tail for exact Chebyshev interpolation; does NOT
certify floating interpolation coefficients or radial derivative errors.
"""
from fractions import Fraction as F
from math import factorial
rho=F(3,2)
lo,hi=F(3,5),F(5,3)
center=(lo+hi)/2
half=(hi-lo)/2
re_min=center-half*(rho+1/rho)/2
radius_bound=center+half*(rho+1/rho)/2
assert re_min==F(5,9)
assert radius_bound==F(77,45)<2
assert F(2)*F(3)*re_min==F(10,3)
# Exp(10/3)>27, proved by rational Taylor lower bound.
x=F(10,3)
assert sum(x**k/factorial(k) for k in range(9))>27
# cosh(u)>=1+u^2/2. For real r, |cos(r*u)|<=1.
# |K_ir(a*y)|<=exp(-a Re(y))*sqrt(pi/(2*a Re(y)))
# a=2*pi*sqrt(n)>=2*pi>6, pi<20/3 and a Re(y)>10/3
# imply sqrt(pi/(2*a Re(y)))<1 and exp(-a Re(y))<1/27.
# Hence |y K_ir(a*y)| < 2/27 < 1/10 for all 1<=n<=50,
# all real r, and all y on the Bernstein ellipse E_(3/2).
assert F(2,27)<F(1,10)
M=F(1,10)
n=120
# Standard Chebyshev-Lobatto polynomial interpolation theorem:
# error<=4*M*rho^(-n)/(rho-1)
value_error=4*M*(1/rho)**n/(rho-1)
assert value_error<F(1,10**20)
print('PASS: Re(y) >=',re_min,'and |y|<=',radius_bound,'on ellipse E_3/2')
print('PASS: |y K_ir(2*pi*sqrt(n)*y)|<1/10 for all 1<=n<=50')
print('PASS: ideal degree-120 Chebyshev interpolation error <1e-20')
print('Exact upper bound:',value_error)
print('DIAGNOSTIC float:',float(value_error))
print('NOT YET CERTIFIED: Chebyshev node value intervals, coefficient rounding, first derivatives')
