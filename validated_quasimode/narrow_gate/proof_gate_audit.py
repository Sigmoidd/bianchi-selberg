#!/usr/bin/env python3
"""Exact-arithmetic deduction of pointwise gates for a narrow odd-spectrum enclosure.
No interval supremum claims are made by this file.
"""
from fractions import Fraction as F
from math import isqrt

# F = {-.5<=x<=.5,0<=t<=.5, y>=sqrt(1-x^2-t^2)};
# core y<=5/3. For u=x^2+t^2<=1/2:
# 1/(1-u) <= 1+2u (since u(1-2u)>=0).
# x,t area=1/2, integrals x^2 and t^2 both 1/24.
V_upper= F(1,2)*(F(1,2)-F(9,50))+F(1,24)+F(1,24)
assert V_upper==F(73,300)
# Uniform |R| < 2e-10 => R² < V_upper * (2e-10)² < 1e-20.
assert V_upper*F(2,10**10)**2 <F(1,10**20)
# Horizontal strip 1<=y<=5/3 has area .5 and integral y^-3 = 8/25.
V_strip=F(1,2)*F(8,25)
assert V_strip==F(4,25)
# Uniform defect |D_g| <7e-6 => ||W-f||² < (7e-6)² *4/25.
# sqrt(V_strip)=2/5 exactly, so
seed_norm=F(13,10**6)
defect_pointwise=F(7,10**6)
mu_lower=seed_norm-F(2,5)*defect_pointwise
assert mu_lower==F(102,10**7) # 1.02e-5
# With R<1e-10 and mu>1.02e-5, ratio <1e-5.
assert F(1,10**10)/mu_lower <F(1,10**5)
print('PASS exact core volume:',V_upper,'<1/4')
print('PASS exact strip volume:',V_strip)
print('PASS uniform-defect plus seed-mass gate: ||W|| >',mu_lower)
print('PASS uniform-residual gate: R <1e-10; ratio <1e-5, CONDITIONAL on two UNPROVED sup inequalities')
