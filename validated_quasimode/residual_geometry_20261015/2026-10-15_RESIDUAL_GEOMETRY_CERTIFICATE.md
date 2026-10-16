# Picard analytic quasimode: verified support and symmetry reductions

Date: 2026-10-15. Branch: `research/level-one-analytic-quasimode` (draft PR #14). **Status: PARTIAL. A rigorous two-sided spectral enclosure is not yet obtained.** This file supersedes the *numerical* geometric support speculation but does not supersede the older mass certificates.

## New exact structural theorems

The projected trial from `analytic_quasimode_coefficients_50_reflection_even.json`, SHA256 `f306a3740651575df5be1b5e5130eda956e346b5a3905579b48556eb31d78f98`, obeys `W(-x,t,y)=W(x,t,y)` and `W(t,x,y)=-W(x,t,y)`. These identities follow exactly from the quarter-turn odd action and the reflection-even projection. The same identities apply to the residual `(L-lambda)W` because `L` commutes with these isometries. Consequently, the supremum of `|R|` and `|W-f|` over the Ford cell can be taken on the wedge `0<=t<=x<=1/2`. Both functions vanish on its diagonal `t=x`.

Of the nine possibly active Ford representatives from the previous exact pruning, just **six** can have `psi_gamma>0` on this wedge: representative indices `(0,1,2,4,5,13)`, with zero representing the identity. The three additional inactive representatives `(7,8,10)` have `h_gamma<=1/2<3/5`: for `c=1,d=1` and `c=1,d=1-i`, `|cz+d|>=1+x>=1`; for `c=1+i,d=1`, `|cz+d|>=1+x-t>=1`, and `|c|=sqrt2`. The full nine-representative argument remains valid globally on the Ford domain, and we only use the six-representative formula after invoking the symmetry reduction.

The residual vanishes identically wherever all nonidentity cutoffs are locally constant (each 0 or 1). Its support is therefore contained in the union of transition bands `3/5<h_gamma<7/10`. The exact integer-only program `certify_transition_support.py` covers the Ford cell by rational cuboids, retains all cuboids that *might* intersect one such transition, and exactly sums their full hyperbolic volume, including their portions that lie outside the Ford cell.

Its 256-by-128-by-384 run gives a rigorous upper bound `0.1358388152300871...` (decimal for orientation), proving the simple **exact rational** statement

    volume(support R) < 11/80.

An independent coarser 128-by-64-by-192 cover gives `0.141766560746...`. Both replays use only bounded signed `int64` arithmetic for geometric predicates; maximum integer expressions are many orders smaller than 2^63. The final volume sum is calculated using Python `fractions.Fraction`, without floating-point error. The exact rational result is in each corresponding JSON file. All boxes crossing the curved Ford floor are retained rather than guessed to be inside or outside.

## Improved spectral closure gates, STILL UNPROVED

Previously a MPFR calculation certified `||f||²_{1<=y<=5/3}>1.7682998477e-10` for the projected 40-orbit seed. In particular `||f||_strip > 1.32e-5`, slightly improving the formerly reported rational bound `1.3e-5`.

It now suffices to verify *directly*, by validated uniform interval arithmetic on `0<=t<=x<=1/2`, the two inequalities

    sup |W-f| < 1e-6   on 1<=y<=5/3,
    sup |(L-lambda)W| < 3.4e-10   on sqrt(1-x^2-t^2)<=y<=5/3.

The direct correction gate avoids taking a maximum over all individual automorphy defects. Because `volume(strip)=4/25`, its first inequality would prove

    ||W|| > 1.32e-5 - sqrt(4/25)*1e-6 = 1.28e-5.

Because `volume(supp R)<11/80`, the second would imply

    ||R||² < (11/80)*(3.4e-10)².

The *exact* rational script `proof_gate_improved.py` verifies their combination implies `||R||/||W||<1e-5`, hence the existence of a point of the **Q-odd cuspidal spectrum** within 1e-5 of the exact `lambda=1+(6.622119348586011)^2`. It does not show first-eigenvalue status, a spectral exclusion below this location, or a trace-formula certification.

The sampled screen recorded `max |R| ~1.16e-10` and `max |W-f| ~4.38e-13`, but **those screens do not certify either supremum**. The new thresholds were chosen to balance rigorous mass transfer against support geometry. The interval implementation is the remaining blocker.

## Independent replay

```sh
python certify_triangle_symmetry.py
python certify_transition_support.py --resolution 128 --y-den 192 --output support_cover_128.json
python certify_transition_support.py --resolution 256 --y-den 384 --output support_cover_256.json
python proof_gate_improved.py
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_triangle_reduction
```

Next engineering action: build a verified Bessel/Chebyshev radial table for the finite set of norms, covering `[3/5,5/3]` with outward coefficient and derivative errors; then a certified branch-and-bound on the triangular wedge using the exact quadratic cutoff surfaces. Exploit the defect-preserving formula for `R`, and avoid naive interval expressions that separately bound large cancelling values. Log every terminal box and proof of full coverage.

**No theorem status flag should be promoted until both supremum inequalities return `STRICT PASS`.**
