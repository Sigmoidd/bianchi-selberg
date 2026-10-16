#!/usr/bin/env python3
"""Exactly check implications of two UNPROVED supremum gates.

This is NOT itself a full spectral certificate; it checks only algebraic
implications of: (1) projected-seed mass >= 1.32e-5 (existing MPFR proof),
(2) support volume < 11/80 (new exact integer proof), and hypothetical
W-f and residual supremum inequalities to be discharged by intervals.
"""
from fractions import Fraction as F
from certify_triangle_symmetry import TRIANGLE_INDICES,proof

proof()
assert TRIANGLE_INDICES==(0,1,2,4,5,13)
# Published MPFR seed-strip squared lower bound is >1.7682998e-10.
# Use weaker rational floor 1.768e-10 to avoid importing a huge interval.
seed_sq_lower=F(1768,10**13)
seed_norm_target=F(132,10**7)
assert seed_sq_lower > seed_norm_target**2
strip_volume=F(4,25)
assert strip_volume ==F(2,5)**2
# To be checked on triangular fundamental wedge: 0<=t<=x<=1/2
W_minus_f_sup_target=F(1,10**6)
residual_sup_target=F(34,10**11)
# Independent cover shows vol(supp R)<11/80.
support_volume_target=F(11,80)
mu_lower=seed_norm_target-F(2,5)*W_minus_f_sup_target
assert mu_lower==F(128,10**7)
eps_target=F(1,10**5)
assert support_volume_target*residual_sup_target**2 <mu_lower**2*eps_target**2
print('PASS EXACT: uniform mass correction <1e-6 ==> mu>1.28e-5')
print('PASS EXACT: residual supremum <3.4e-10 and support vol<11/80 ==> residual/mu<1e-5')
print('STATUS: both pointwise uniform supremum assumptions UNPROVED')
