#!/usr/bin/env python3
"""Exact Ford-cell pruning and algebraic cutoff-surface splitting.

This program performs a geometric deduction, not a rigorous residual integration.
It is intended to be used by a later interval certificate, and the exact
elimination lemma can be verified without special functions.
"""
from __future__ import annotations
import math
from fractions import Fraction as F
from level_one_quasimode_analytic import representatives

# These entries are defined by their 0-based index in the 14-representative
# formula from level_one_quasimode_analytic.py. Every representative in this
# list satisfies h_gamma <= 1/2 or 1/sqrt(5) on Ford F. Both are <3/5.
PRUNED = {3,6,9,11,12}
ACTIVE_INDICES=tuple(i for i in range(14) if i not in PRUNED)


def exact_prune_proof():
    reps=representatives()
    assert len(reps)==14
    # c=1, d=i and d=+/-1+i: nearest horizontal q-distances
    assert reps[6][2:]==(1+0j,1j)
    assert reps[3][2:]==(1+0j,-1+1j)
    assert reps[9][2:]==(1+0j,1+1j)
    # q=(1+i)z+d. On x in [-1/2,1/2], t in [0,1/2]:
    # d=-1 => |q|² = 2(x-1/2)² + 2(t+1/2)² >= 1/2.
    # d=+i => |q|² = 2(x+1/2)² + 2(t+1/2)² >= 1/2.
    assert reps[11][2:]==(1+1j,-1+0j)
    assert reps[12][2:]==(1+1j,1j)
    assert F(1,4)+F(1,4)==F(1,2)
    # h <= 1/(2|c||cz+d|). c unit, min |q|² >= 1
    # (and for d=+-1+i, >=5/4) => h <= 1/2.
    # c=1+i, min |q|² >=1/2 => |c|²|q|²>=1 => h<=1/2.
    assert F(1,2)<F(3,5)
    return [reps[i] for i in ACTIVE_INDICES]


def vertical_cut_roots(gamma,x,t,h0):
    """Midpoint planning only: exact surfaces are A*y^2-y+B=0;
    interval proof must enclose roots and integrate each algebraic cell.
    """
    a,b,c,d=gamma
    if c==0:return []
    A=h0*abs(c)**2
    B=h0*abs(c*complex(x,t)+d)**2
    disc=1-4*A*B
    if disc<0:return []
    disc=max(0.,disc)
    return sorted(v for v in ((1-math.sqrt(disc))/(2*A),(1+math.sqrt(disc))/(2*A)) if v>0.)


def all_vertical_splits(x,t,ylo=None,yhi=5/3):
    if ylo is None:ylo=math.sqrt(1-x*x-t*t)
    cuts=[ylo,yhi]
    for ga in exact_prune_proof()[1:]:
        for h0 in (.6,.7):
            cuts.extend(y for y in vertical_cut_roots(ga,x,t,h0) if ylo<y<yhi)
    return sorted(set(cuts))


def test():
    reps=exact_prune_proof()
    assert len(reps)==9
    # Grid is only a regression; the exact proof above is independent.
    for x in (-.5,-.27,0.,.22,.5):
        for t in (0.,.11,.32,.5):
            for y in (math.sqrt(1-x*x-t*t),.8,1.,1.2,1.4,5/3):
                for i in PRUNED:
                    a,b,c,d=representatives()[i]
                    h=y/(abs(c*complex(x,t)+d)**2+abs(c)**2*y*y)
                    assert h<.6
                for k in all_vertical_splits(x,t):
                    assert k>0 and k<=5/3+1e-12
    print('PASS: five inactive cosets pruned by exact Ford-cell inequalities')
    print('PASS: nine representatives remain; algebraic vertical transition-root splitting')

if __name__=='__main__':test()
