# Level-one odd quasimode: analytic-residual reconstruction

**Status: EXPERIMENTAL, not an interval-certified spectral enclosure.**

This artifact reconstructs the uploaded `level-one-quasimode.diff` experiment without its sixth-order finite-difference Laplacian. It supplies a complete frozen rational coefficient list (40 C4 orbits at M=50), an analytic chain-rule evaluator, and independent Gauss and Sobol numerical checks. The earlier 40-component numerical optimizer output was missing from the user's diff; the reconstructed coefficients are therefore **new** rather than an exact replay of the prior unnamed vector.

## Fixed trial

- Group: `PSL_2(Z[i])`, odd under `Q(z,y)=(iz,y)`.
- Spectral parameter `r = 6.622119348586011` (decimal rational value for a later exact theorem), `lambda = 1+r^2`.
- Finite Fourier modes `a^2+b^2 <= 50`, grouped into **40** fourfold odd orbits, normalized `a_(1,0)=1`.
- Cutoff: quintic `psi(y)=0` for `y <= 3/5`, `psi(y)=1` for `y >= 7/10`.
- Function: `W = sum_(Gamma_inf\Gamma) (psi f) o gamma / sum_(Gamma_inf\Gamma) psi o gamma`.
- Denominator `G>=1` on the Gaussian Ford cell. `W=f` for `y>=5/3`.
- All coefficients are frozen as exact *decimal rationals* in `analytic_quasimode_coefficients_50.json`, SHA-256 `b04f90ca943e5d6d574f5f7c47d7649c8fdadd6b40607f6cc9064a46cdbb671b`.

## A better finite-coset lemma: **14 instead of 26**

On the reference Ford cell `|z|<=1/sqrt(2)`, `y>=1/sqrt(2)`, a nonparabolic transformation has height

`h=y/(|cz+d|^2+|c|^2*y^2)`.

Only the identity and **13** nontrivial coset representatives can have `h>3/5`:

- `N(c)>=3` gives `h<=sqrt(2)/3<3/5`.
- `N(c)=1`: normalize `c=1`. If `N(d)>=4`, then `|z+d|>=2-1/sqrt(2)` and `h<=1/[2(2-1/sqrt(2))]<3/5`. Thus only `d` of norm `0,1,2` (nine choices) can contribute.
- `N(c)=2`: normalize `c=1+i`. Coprimality excludes `N(d)=0,2,4`; if `d` is not a unit then `N(d)>=5` and `|cz+d|>=sqrt(5)-1`, so `h<=1/[2 sqrt(2)(sqrt(5)-1)]<3/5`. Only four units can contribute.

This gives `1 + 9 + 4 = 14` cosets. The original 26-representative evaluation had 12 always-zero weights. A direct rerun with 14 gave **identical output** to the prior 26-coset analytic evaluation on the tested grids.

## Analytic residual used

Write `h_g = y(g P)`, `w_g=psi(h_g)`, `f_g=f(g P)`, and `G=sum_g w_g`, `N=sum_g w_g f_g`, `W=N/G`.

For the positive hyperbolic Laplacian `L=-y^2 (dxx+dtt+dyy)+y dy`, each finite Fourier-Whittaker term satisfies `(L-lambda) f_g=0`.

`(L psi(h_g))= h_g psi'(h_g)-h_g^2 psi''(h_g)` and `|grad h_g|_H^2=h_g^2`.

The implemented *analytic* residual is

`(L-lambda)W = [sum_g {(Lw_g)(f_g-W)-2<grad w_g,grad f_g>_H} + 2<grad G,grad W>_H]/G`.

The derivatives of `f(gP)` use the **full explicit real 3x3 Jacobian** of the hyperbolic Möbius action; radial derivatives come from differentiated Chebyshev interpolants of `y K_(i r)(2 pi |mu| y)`. No finite differences enter the reported residual integration. The Chebyshev interpolants and Bessel values are **not interval-certified**; the reported midpoint interpolation error is only diagnostic.

The independent scalar finite-difference test shows 2nd-order truncation errors `1.5e-8 -> 3.9e-9` at h `.004 -> .002`, far larger than the true residual scale `1e-11`, exactly why old numerical differentiation was unsafe.

## Independent numeric evidence

| Test | Nodes | Residual R | Ratio R / mu |
|---|---:|---:|---:|
| Reconstructed optimizer, training | 1120 | 1.18944e-11 | 6.2346e-7 |
| New frozen coefficients, Gaussian held-out | 3696 | 1.18720e-11 | 6.2211e-7 (including cusp tail) |
| New frozen coefficients, Gaussian finer | 8960 | 1.19113e-11 | 6.2434e-7 (core mass only) |
| New frozen coefficients, scrambled Sobol seed 17 | 8192 | 1.18876e-11 | 6.2310e-7 (core mass only) |

Mass including numerically integrated tail on the 3696-point grid is about `1.90834357e-5`; the tail alone has `mu_tail≈4.46946e-7`, **not** the `>=3e-6` claimed in the uploaded plan. Thus both the cusp-only lower bound and any earlier claim of a completed numerical spectral certificate must be corrected.

For the *particular* frozen finite trial, a verified interval enclosure `R/mu < 0.05` would imply existence of an eigenvalue in the global odd cuspidal sector within `0.05` of the stated lambda. This implication is a theorem; the antecedent has **not been certified**.

## Replay

```
python -m pip install numpy scipy mpmath
python -m unittest test_analytic_quasimode -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python level_one_quasimode_analytic.py --M 50 --degree 118 --nx 10 --nt 8 --ny 14 --verify-grid 14,12,22 --output analytic_quasimode_50.json --coefficients analytic_quasimode_coefficients_50.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python level_one_quasimode_verify.py --coefficients analytic_quasimode_coefficients_50.json --degree 140 --gauss 20,16,28 --sobol 8192 --output analytic_quasimode_independent_fine.json
```

The first run is deterministic up to numeric library versions and should regenerate the same coefficients; the second script does no fitting and uses only frozen coefficients.

## The last proof obligations

1. Convert the 40 frozen coefficient decimal strings and `r` to exact rationals in the theorem record (already parseable). Do not replace them with a new fit during certificate production.
2. Prove interval majorants for the 14 cosine cosets and the analytic commutator over the **complete curved Ford core**. Use cancellation-preserving Taylor/Chebyshev models with rational box faces and separately integrated squared polynomial pieces, not direct interval bounds on each big Fourier term.
3. Prove a positive interval mass lower bound. The cusp four-mode Parseval integral can be bounded rigorously with Arb/Acb and a validated radial quadrature on a finite slab; this avoids relying on a mistaken `3e-6` tail mass assumption.
4. Independently reconstruct the same frozen trial with a separate Bessel implementation and verify the 3D Jacobian, boundary/elliptic compatibility, and global automorphy.
5. Only after those steps issue `certified=true` and a spectral interval. `rung4_certified` and `hejhal_existence_certified` remain **false**.