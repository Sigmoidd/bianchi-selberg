# One-profile level-one Picard variational certificate: proof and replay

Prepared 2026-10-14. **Status:** strict negative floor flux enclosed using an explicit quadrature-error theorem and `mpmath.iv` interval arithmetic, reproduced with two independently sized rules. This is a substantive computer-assisted proof candidate, **not yet independently replayed using Arb/MPFR**. Do not elevate unrelated `rung4_certified` or dual-certification flags. The theorem discussed here is **existence** of a level-one quarter-turn-odd cuspidal eigenvalue *below 45*, not its location to narrow width or primacy in the entire quotient spectrum.

## 1. The single-profile trial

Let `Gamma=PSL_2(Z[i])` and use the classical Picard/Humbert cell

\[
\mathcal F=\left\{(x,t,y): -\tfrac12\le x\le\tfrac12,\quad
0\le t\le\tfrac12,\quad x^2+t^2+y^2\ge1\right\}.
\]

Fix the exact rational parameter
\[
r=6633/1000,\qquad \lambda_0=1+r^2=44996689/1000000=44.996689<45.
\]

Define the *single* Fourier profile (the factor 2 is an arbitrary positive scaling)

\[
f(x,t,y)=2yK_{ir}(2\pi y)\,[\cos(2\pi x)-\cos(2\pi t)].
\]

It solves `(Delta-lambda0)f=0` throughout the cell for
`Delta=-y^2(partial_x^2+partial_t^2+partial_y^2)+y partial_y`.
The integral representation of `K_(ir)` ensures exponential decay, so `f` and its first derivatives are square-integrable in the cusp.

### Why it is an admissible global H¹ test without a Poincaré sum

A trial for the variational principle need only have matching **values** under every quotient face pairing. Matching normal derivatives is not required for membership in the quadratic-form domain.

The five Picard face identifications are:

* `T_1`: x=-1/2 is translated to x=+1/2; both cosines are periodic.
* `R`: z -> -z on the t=0 face; f is even in x and t.
* `T_i R`: z -> i-z on the t=1/2 face; f is even in x.
* `S`: inversion `z -> -1/z` in H³; on the Ford hemisphere `x²+t²+y²=1`, S(x,t,y)=(-x,t,y); f is even in x.

Thus the traces match and define an `H¹(Gamma\H³)` form-domain vector. Stabilizer incidences obey these same value identifications. The unique cusp constant Fourier mode vanishes. Under `Q(z,y)=(iz,y)`, `f` changes sign by swapping x and t. After reduction of QP back to the reference half-cell (using R when x<0), this remains true globally. The odd cusp sector has compact resolvent, so its Rayleigh infimum is an eigenvalue.

In particular, **there is no need for a cutoff, a normalized Poincaré average, a six-copy partition, or a second derivative of f.**

## 2. Green's identity leaves one two-dimensional integral

With hyperbolic measure `dV=dx dt dy/y³` and form `E=∫ |grad_E f|² dx dt dy/y`, Green's identity on the truncated Ford cell and exponential decay imply

\[
E(f)-45\|f\|^2=(\lambda_0-45)\|f\|^2+B_{\rm floor}.
\]

The vertical-face fluxes cancel exactly by the stronger local periodicity/reflection properties. The only remaining boundary is the self-paired inversion hemisphere. Its **outward Euclidean** normal is `(-x,-t,-y)` and its Euclidean area measure is `dx dt / y`. Therefore

\[
B_{\rm floor}=-\int_{-1/2}^{1/2}\int_0^{1/2}
\frac{f(x,t,y)\,[x f_x+t f_t+y f_y]}{y^2}\,dt\,dx,
\quad y=\sqrt{1-x^2-t^2}.
\]

Put `U=cos(2πx)-cos(2πt)`, `D U=2π[-x sin(2πx)+t sin(2πt)]`,
`z=2πsqrt(1-x²-t²)`, `K=K_(ir)(z)`, and `K'=dK_(ir)/dz`. Exact differentiation yields

\[
\boxed{B_{\rm floor}=-4\int_{-1/2}^{1/2}\int_0^{1/2}
\left[K^2(U^2+U DU)+zKK'U^2\right]dt\,dx.}
\]

Since `lambda0-45=-3311/1000000<0` and f is nonzero, **it is sufficient to prove `B_floor<0`**. No mass integration or lower bound is needed.

## 3. A validated-Bessel evaluator requiring only real elementary functions

For real `z>0` and real r,

\[
K_{ir}(z)=\int_0^\infty e^{-z\cosh u}\cos(ru)du,\qquad
K'_{ir}(z)=-\int_0^\infty \cosh u\, e^{-z\cosh u}\cos(ru)du.
\]

The evaluator truncates at `u=6`, using validated 80/88-point Gauss–Legendre rules with validated node and weight brackets. For `z≥2π/sqrt(2)>4`, the omitted tail is bounded by `1e-100` using `cosh 6>200` and the exact `s=cosh u` substitution.

For Bessel integration, the Bernstein-ellipse radius is `rho_u=13/10` on the transformed interval `[0,6]`. On that ellipse `|Im u|<4/5` and `|Re u|<611/100`. The crude bounds

* `|cos(r u)|<102`, `|cosh u|<226`,
* `Re[z cosh(u)] > 8/3`, hence `|exp(-z cosh u)|<1/14`,

imply uniform analytic modulus `<2000` for the derivative integrand, `<100` for the value integrand. The standard Chebyshev approximation argument, using exactness of n-point Gauss quadrature through degree `2n-1`, gives (in the `[0,6]` measure)

\[
\epsilon_K,\epsilon_{K'}\le
160000(10/13)^{2n-1}+10^{-100}.
\]

At `n=80` this is smaller than `1.23e-13`. Both K and K' intervals are independently widened by `1e-12` per spatial node, a deliberately conservative choice.

## 4. A validated two-dimensional floor integration

For x on `[-1/2,1/2]` and t on `[0,1/2]`, use a tensor Gauss rule with `n=56` or `n=64` in each coordinate, and Bernstein radius `rho_f=7/5` in each transformed coordinate.

On these complex product ellipses,

* `|Re x| ≤ 37/70`, `|Re t| ≤ 18/35`, hence `Re(1-x²-t²) ≥ 2235/4900 > 1/4`;
* the principal square root satisfies `Re y>1/2` and `|y|<2`;
* `|Im x|≤6/35`, `|Im t|≤3/35`, giving `|U|<4`, `|DU|<15`;
* `Re z>π`, and the integral representation plus Gaussian inequalities on cosh give `|K_(ir)(z)|<1/20`, `|K'_(ir)(z)|<1/20`, and `|z|<13`.

Thus the entire **complex-valued** analytic floor-flux integrand has modulus less than 100. Tensor Gauss quadrature is therefore bounded by the deliberately loose remainder

\[
\epsilon_{\rm floor}\le 10000(5/7)^{2n-1}<10^{-12}
\quad\text{for } n\ge56.
\]

The bound follows by truncating each coordinate's Chebyshev expansion at degree `2n-1`, bounding its tail by geometric decay, applying exact Gauss integration to the truncated polynomial, and using positivity/sum of the Gauss weights. It holds in the full analytic product ellipses, not merely at the sampled nodes.

All 192 root intervals in the first run (216 in the second) were validated by strict interval signs of the Legendre polynomial at the rational decimal endpoints; intervals were checked to be disjoint and ordered, and weights were independently enclosed by the exact root-derivative expression. The sum-of-weights identity was enclosed as a consistency check.

## 5. Two independent-resolution interval replays

| Spatial × Bessel nodes | Gauss sum enclosure (approximately) | Full safe floor-flux interval (approximately) |
|---|---|---|
| 56×56, 80 | `[-3.36792710,-3.36779499] × 10^-11` | `[-3.46792710,-3.26779499] × 10^-11` |
| 64×64, 88 | `[-3.36792763,-3.36779446] × 10^-11` | `[-3.46792763,-3.26779446] × 10^-11` |

Both independently prove, **relative to the declared interval-arithmetic and analytic-error implementation**,

\[
B_{\rm floor}<-3.2677\times10^{-11}<0.
\]

Together with the explicit `H¹` form-domain argument and `lambda0<45`, this yields the intended variational existence implication.

## 6. Remaining independent proof verification

The executable certificate uses `mpmath.iv` (high-precision directed intervals). Its endpoint checks, error terms and exact formula can be replayed today, but the library is not the established python-flint/Arb environment used by the rest of the repository. A publication-grade statement should require:

1. independent verification of the two complex ellipse modulus bounds and the `10000*(5/7)^(2n-1)` remainder coefficient;
2. a repeat of the interval calculation in python-flint Arb/Acb or MPFR with explicitly directed rounding;
3. an independently written check of the Picard fundamental-domain face identifications and Green's theorem orientation;
4. a preserved exact hash manifest and refusal to reuse these results for an unrelated trial or spectral target.

The historical trace-formula exclusions, higher-level congruence problem, newspace statements, and first-eigenvalue/counting claims are **not** consequences of this certificate.

## 7. Replay commands

From the directory containing the standalone script:

```bash
python certify_single_profile_floor.py --space-n 56 --bessel-n 80 --bits 320
python certify_single_profile_floor.py --space-n 64 --bessel-n 88 --bits 320 \
  --out single_profile_floor_interval_independent_result.json
python -m unittest test_single_profile_variational -v
```

The two interval JSON files contain the full outward-enclosed values, declared quadrature errors and `sign_negative` gates. Preserve the scripts and JSON together. Do not turn this existence proof into a claim that the numerical eigenvalue `44.85246...` has been rigorously localized.
