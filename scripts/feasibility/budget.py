"""Exploratory budgets; positive budget does not prove a spectral bound."""
import mpmath as mp, sys
from screen import field, mech_terms
mp.mp.dps = 20
print(f"{'d':>4} {'vol':>8} {'l0':>8} | k  {'delta':>7} {'g(0)':>7} {'B_mech':>9} {'NCE budget':>11} {'C_ell budget':>13}")
for d in [2, 7, 11, 19, 43, 67, 163]:
    F = field(d)
    for k in (2, 3):
        T = mech_terms(F, k)
        budget = 1 - T['Bmech']; cb = budget / T['g0']
        print(f"{d:>4} {mp.nstr(F['vol'],6):>8} {mp.nstr(T['ell0'],6):>8} | {k}  {mp.nstr(T['delta'],5):>7} {mp.nstr(T['g0'],5):>7} "
              f"{mp.nstr(T['Bmech'],6):>9} {mp.nstr(budget,6):>11} {mp.nstr(cb,6):>13}")
print()
print("reference per-class costs in C_ell units (N(T0)=7+4sqrt3, log N=%s):" % mp.nstr(mp.log(7+4*mp.sqrt(3)),6))
for lab, den in [("order-3, |E|=3 (log N/9)", 9), ("order-2, |E|=2 (log N/8)", 8), ("order-2, hypothetical |E|=4 (log N/16)", 16)]:
    print(f"   {lab:38s} {mp.nstr(mp.log(7+4*mp.sqrt(3))/den,6)}")
print("known totals: Z[i] C_ell=0.29266 (+CE),  Z[omega] C_ell=0.32924 (+CE)")
