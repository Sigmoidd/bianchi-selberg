# A referee-style rethink of Picard quasimode certification

**Date:** 15 October 2026. **Status:** A NEW, CONDITIONAL, TWO-DIMENSIONAL EIGENVALUE INCLUSION THEOREM, supported by strong independent numerical evidence. **No certified two-sided enclosure yet.**

## Executive finding

The costly 3D sup norm of the residual of the Poincaré-average `W` is not necessary for a first narrow spectral inclusion. The Fourier–Bessel seed `f` **already solves** `(L-lambda) f=0` on the Ford interior and is continuous across its identified faces. It is in the automorphic **form domain** `H^1` but generally not the self-adjoint **operator domain** `Dom(L)` because its normal derivatives do not match across the paired floor.

Work with the *weak residual* in the dual of the quadratic-form domain. Green's formula expresses it entirely as a single boundary distribution supported on the Ford hemisphere. This is the analogue of the *method of particular solutions* / boundary-tension eigenvalue inclusion, except we provide an elementary explicit trace bound specialized to this Picard geometry instead of importing a Euclidean-domain theorem.

## Frozen data

- `Gamma=PSL_2(Z[i])`, half-cell `F={(x,t,y): -1/2<=x<=1/2, 0<=t<=1/2, y>=sqrt(1-x²-t²)}` in the reflection-even sector.
- `r=6622119348586011/10^15`, `lambda=1+r²` **exactly**.
- `f` is the real 21-term cosine-difference projection of the *unchanged* 40-orbit source coefficient vector in `inputs/`. `f` is Q-odd, periodic in x,t, and even in x and t. Reflection-even symmetry explains the half-cell t range; do not mistake the reflection for a group element of PSL_2(Z[i]).
- `f` satisfies `(L-lambda)f=0` pointwise inside the cell. Values match under the inversion `S` on the hemispherical floor because `S(x,t,y)=(-x,t,y)` there and `f` is x-even. The paired x faces are periodic, and the t=0,1/2 reflection faces have zero normal derivative. At infinity the Fourier–Bessel seed decays.
- Existing separately replayed MPFR height-slab certificate gives `||f||_{L²(F)}>13/10^6` for the same projected seed (requires an independently checked mass normalization and face-orbifold convention in a referee audit).

## Boundary formula

Let `y0=sqrt(1-x²-t²)`, `D=[-1/2,1/2] ×[0,1/2]`. Let `dS_H=dx dt / y0³` and outward hyperbolic normal derivative

`J(x,t)=∂_{n_H}f=-y0[x f_x+t f_t+y0 f_y]` at `y=y0`.

For every compatible Q-odd, reflection-even `H¹` test `v`, hyperbolic Green's identity yields

`a(f,v)-lambda(f,v) = ∫_D J(x,t) conjugate(v(x,t,y0)) dx dt/y0³`, 

with `a(u,v)=∫ y^{-1} (u_x v_x+u_t v_t+u_y v_y) dx dt dy`. This is a *single-face* integral: do **not** double J again merely because S pairs two points on the floor. The reflected quotient identifies pairs already contained in D.

## Elementary trace constant

For fixed (x,t), set `h=y0>=1/sqrt(2)` and `delta=1/3`. The 1D fundamental theorem of calculus and Cauchy–Schwarz give

`|v(h)|² <= 2/delta ∫_h^{h+delta} |v(y)|² dy + 2delta∫_h^{h+delta}|v_y|²dy`.

Converting to hyperbolic weights gives

`||v|floor||²_{L²(dS_H)} <= (81/4) (||v||²_{L²(F)} + a(v,v))`,

since `sqrt(2)<3/2` and `max{2/delta(1+sqrt(2)delta)^3, 4delta(1+sqrt(2)delta)} <81/4`. In the second term we used `h^{-2}<=2`.

Thus the dual residual satisfies `||R_weak||_{H^{-1}_form} <= (9/2)||J||_{L²(floor)}`.

## Spectral theorem implication

The Q-odd (and reflection-even) subspace has no cusp constant term and is spectrally discrete. Set `d=dist(lambda,spectrum(L)|sector)`. In an orthonormal eigenbasis, 

`||R_weak||²_{V*} = ∑ |(lambda_j-lambda)|²/(1+lambda_j) |f_j|²`.

If `d>=epsilon` and `lambda>=0`, the summand ratio is >= `epsilon²/(1+lambda+epsilon)`, implying 

`epsilon/sqrt(1+lambda+epsilon) ||f|| <= ||R_weak||_{V*}`.

For `epsilon=10^-5`, `sqrt(1+lambda+epsilon)<7`. Consequently an outward-verified bound **`||J||_{L²(floor)}<4×10^-12`**, combined with `||f||>13×10^-6`, forces `d<10^-5`. (Pure rational check: `7*(9/2)*4e-12 < 13e-6 *1e-5`.) A stronger, simpler sufficient condition is `sup_D |J|<4×10^-12` since the area satisfies `area_H(floor)<5/6<1` (convexity of `(1-u)^(-3/2)` gives `(1-u)^(-3/2)<=1+4u` on `0<=u<=1/2`).

**This is a valid conditional inclusion theorem; it is not a certified value for `||J||` yet.** It proves existence of *some* Q-odd eigenvalue in the `epsilon`-disk on the real line, not its index, multiplicity or uniqueness. Proof of the compact resolvent for the Q-odd sector and all quotient face pairings must be referee-checked.

## Independent numerical checks (NOT interval bounds)

Script `boundary_weak_residual.py` reconstructs the 21-term seed from *exact rational* frozen coefficient JSON and evaluates first derivatives using the independently generated radial Chebyshev table:

| Gauss nodes in x,t | `||J||_L²` | `∫fJ dS_H` |
|---|---|---|
| 32×28 | 8.3168436949e-13 | -4.0193118141e-17 |
| 60×50 | 8.3168436602e-13 | -4.0193117957e-17 |
| 90×70 | 8.3168436570e-13 | -4.0193117942e-17 |
| 120×96 | 8.3168436810e-13 | -4.0193118042e-17 |

The boundary flux combined with the previous numerical norm `~1.90834e-5` predicts Ritz `44.85246455655`, independently reproducing the `44.85246456` numerical review. This is a valuable identity/normalization cross-check, not an exact spectral bound.

A 128×128 Gauss-Chebyshev interpolant of `J(x,t)` on the FULL rectangle, not only the wedge, has **numerical sum of absolute Chebyshev coefficients** `3.054019e-12`. This already falls below the uniform target `4e-12`, leaving ~0.946e-12 for interval nodal rounding and analytic interpolation error.

## New constructive interval proof plan: only 2D Chebyshev coefficients

Map `x=s/2`, `t=(1+u)/4`; use first-kind Gauss-Chebyshev nodes `s_j=cos(pi(j+1/2)/128)`, same for `u`. The complex polyellipse with radii `rho=3/2` admits:

`|x|<=13/24, |t|<=25/48, |Im x|<=5/24, |Im t|<=5/48`, and `Re sqrt(1-x²-t²)>13/20`, `|sqrt(1-x²-t²)|<13/10`.

From the integral representation of the complex-order Bessel K (real r), `cosh v>=1+v²/2`, and `log cosh v <= cosh v-1`, we get `|K_ir(z)|<1/40`, `|K_ir'(z)|<1/40` when `Re(z)>39/10`. Use `|2 pi sqrt(n)|<51`, `sum|c_ab|<37`, `max(a+b)<=9`, `max a<=7`, and `exp(11)<60000` to obtain the very conservative **complex-analytic bound `sup_polyellipse |J|<20,000,000`**. File `audit_theorem.py` checks the necessary rational inequalities (including the exponential bounds by exact rational Taylor series).

For tensor Gauss-Chebyshev interpolation of degree 127 in each variable, the full function has a convergent double Chebyshev series with coefficients `|a_kl| <= 4M rho^(-k-l)`. Higher modes alias at the sampled nodes into lower modes, so `|T_k T_l - aliased_term| <=2`. Therefore

`||J-P_127,127||∞ <= 16 M rho^(-128)/(1-rho^-1)² < 8.32×10^-14 <10^-13`,

with M=20,000,000. This error bound is *analytic*, independent of cancellation in J.

**Final tiny gate:** compute the 128×128 numerical interpolation coefficients as exact outward MPFR intervals, using the earlier rigorously bounded real `F_n` and `F'_n` Chebyshev polynomials, all 128² algebraic cosine nodes and directed-rounding DCTs. Verify

`sum_(0<=k,l<128) sup |coefficient_interval_kl| < 3.5×10^-12`.

Then `sup |J|<3.5e-12+1e-13<4e-12`; the conditional spectral theorem above closes the `10^-5` inclusion. This task uses 2D quadrature/interpolation, **not** the previous 3D Taylor sweeping or Poincaré averaging. 

Important: the current 128×128 coefficient sum was computed using IEEE double and is NOT a certified outward-MPFR coefficient enclosure. Neither the numerical DCT nor its convergence to a stable answer is a theorem. The analytic polyellipse envelope has been derived here but should be independently scrutinized before promotion.

## What an expert would question immediately

1. **Face pairings and orbifold symmetry:** Do the assumed coset identifications, t-reflection sector, H¹ trace conventions, and floor S-map produce exactly the claimed weak residual? Check the operator's Friedrichs form on the quotient, including stabilizers.
2. **The inherited mass:** Re-derive normalization and Parseval factors independently for the precise *reflection-projected* 40-orbit vector; check its SHA and coordinate convention.
3. **A posteriori inclusion vs completeness:** Weak-residual estimates certify proximity to *some* eigenvalue, not that it is the first, only, or exactly located within a much smaller interval.
4. **Validated Chebyshev:** Confirm no unaccounted complex Bessel singularities inside the polyellipse, the 2D aliasing factor, and the interval Fourier/Chebyshev DCT normalization.
5. **Independent solver:** Compare to a conventional finite-element or boundary-MPS solve with a deliberately different basis. Different grids of the *same code* are not independent verification.
6. **Literature/priority:** Child (2022), Booker–Strömbergsson–Venkatesh (2006), Barnett–Hassell (2011) and Dahne–Salvy (2020) establish nearby techniques, but not automatically this precise hyperbolic-threefold theorem. Avoid claiming a new general method without comparison.

## References

- Booker, Strömbergsson & Venkatesh (2006), *Effective computation of Maass cusp forms*.
- Kieran Child (2022), *Certification of Maass cusp forms of arbitrary level and character*, arXiv:2204.11761.
- Barnett & Hassell (2011), *Boundary quasi-orthogonality and sharp inclusion bounds for large Dirichlet eigenvalues*.
- Dahne & Salvy (2020), *Computation of Tight Enclosures for Laplacian Eigenvalues*.

## Status and next gate

**Completed this turn:** exact conditional weak inclusion algebra; conservative trace constant; 2D analytic Chebyshev-error bound; independently repeated boundary tension, flux and 128² DCT diagnostics. **Still missing:** a complete outward-rounded MPFR 128² nodal evaluation and interval DCT coefficient ℓ1 certificate; formal audit of quotient trace/face gluing. No theorem flags should be changed.
