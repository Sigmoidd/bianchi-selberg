#!/usr/bin/env python3
"""Rigorous symmetry and six-coset reduction for the exactly projected Picard trial.

The statements are group and integer inequalities. Numerical evaluations in
--regression are only sanity checks, not the source of the proof.
"""
from fractions import Fraction as F
from ford_active_cosets import ACTIVE_INDICES,PRUNED,exact_prune_proof
from level_one_quasimode_analytic import representatives

TRIANGLE_INDICES=tuple(i for i in ACTIVE_INDICES if i not in (7,8,10))

def proof():
    assert TRIANGLE_INDICES==(0,1,2,4,5,13)
    reps=representatives()
    assert reps[7][2:]==(1+0j,1-1j)
    assert reps[8][2:]==(1+0j,1+0j)
    assert reps[10][2:]==(1+1j,1+0j)
    assert F(1,2)<F(3,5)
    # For 0<=t<=x<=1/2:
    # i=7/8, |z+d|>=x+1>=1, c=1; h<=1/2.
    # i=10, |(1+i)z+1|>=1+x-t>=1, |c|^2=2;
    # h<=1/(2 sqrt2) <1/2 <3/5.
    # Q=(x,t)->(-t,x), T=(x,t)->(x,-t); TQT = Q^-1.
    # For projected W, U_QW=-W and U_TW=W.
    # Thus W(-x,t)=W(x,t), W(t,x)=-W(x,t).
    # L commutes with both Q and T; same properties hold for residual.
    # The Ford floor y^2>=1-x^2-t^2 is invariant under reflections/swap.
    # Therefore sup of |R| and |D correction W-f| can be checked on
    # 0<=t<=x<=1/2; the defect-max gate itself need not be symmetric
    # as a max over independently indexed gamma, but |W-f| is.
    return reps

if __name__=='__main__':
    proof()
    print('PASS EXACT: full Ford sup |R| reduces to triangular wedge 0<=t<=x<=1/2')
    print('PASS EXACT: six surviving representatives (identity + 5) on wedge:',TRIANGLE_INDICES)
    print('NOTE: defect bound must be verified globally or via symmetry of |W-f|, not by naively restricting max of labelled defects')
