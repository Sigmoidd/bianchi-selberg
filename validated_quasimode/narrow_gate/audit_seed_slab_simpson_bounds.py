#!/usr/bin/env python3
"""Exact rational audit of all remainders for [1,5/3] seed mass."""
from fractions import Fraction as F
from math import factorial
P=lambda z,r:(z**4+6*z**3+7*z**2+z)+4*r*(z**3+3*z**2+z)+6*r*r*(z*z+z)+4*r**3*z+r**4
assert P(11,7)<128000
assert F(128000,400)<=320
# e^6>400 from a FINITE rational Taylor lower bound.
assert sum(F(6**k,factorial(k)) for k in range(22))>400
# y derivative bound from earlier exact expansion.
assert F(3,400)*(7**4+4*7**3+12*7**2+24*7+24)<35
# u Simpson: span 3 and h=3/1024, integrated over y slab width 1/12.
w=F(1,12);err_u=w*F(3,180)*F(3,1024)**4*320
err_y=w*F(1,180)*F(1,384)**4*35
assert err_u+err_y <F(9,10000000000) / 10 # 9e-11
# u tail: cosh(3)>10 and sinh(3)>9; for z>=6,
# integral tail <= e^-60 /54. e^60 > 60^14/14! >1e20.
assert F(60**30,factorial(30))>F(10**20)
assert w*F(1,54*10**20)<F(1,10**20)
assert err_u+err_y+F(1,10**20)<F(1,10**10)
print('PASS: e^6>400; u4<320; y4<35; Simpson + tail <1e-10 for every [1,5/3] slab')
print('u error upper',float(err_u),'y error upper',float(err_y))
