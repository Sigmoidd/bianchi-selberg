# Picard global Q-odd 40-orbit quasimode: rigorous cusp mass and residual next gate

**Date:** 2026-10-14  
**Target:** `PSL_2(Z[i])`, odd sector, a global smooth Poincaré-averaged quasimode near `lambda = 44.85246466691722`.  
**Proof status:** The **cusp-mass inequality** is independently replayable with MPFR and exact rational bounds. **The residual is not interval-enclosed**. No spectral enclosure is claimed.

## Fixed inputs

- Frozen finite seed: `analytic_quasimode_coefficients_50.json` (40 fourfold Fourier orbits, norms <= 50), SHA256 `b04f90ca943e5d6d574f5f7c47d7649c8fdadd6b40607f6cc9064a46cdbb671b`.
- `r=6.622119348586011` treated as an **exact decimal rational**.
- The norm-one mode coefficients are exactly `(1,-1,1,-1)` in the C4 orbit, so on the nonoverlapping cusp `y>=5/3`, the norm-one contribution equals `2y K_{ir}(2πy)(cos(2πx)-cos(2πt))`.
- Over the Ford half-square `x in [-1/2,1/2], t in [0,1/2]`, horizontal Fourier orthogonality gives `mu_tail^2 >= 2 ∫_{5/3}^∞ K_{ir}(2πy)^2/y dy`. All other 39 Fourier orbits contribute nonnegative mass, and their coefficients never enter this proof.

## Exact finite-slab lower bound

For each interval `[a,b]`, Cauchy–Schwarz states

```
2 ∫_a^b K_ir(2πy)^2/y dy >= 2 [∫_a^b K_ir(2πy)/y dy]^2 / log(b/a).
```

Set `I1=[5/3,7/4]`, `I2=[7/4,19/10]`, `I3=[19/10,21/10]`. For every `I=[a,b]`, write

```
J_I = ∫_a^b dy/y ∫_0^∞ exp(-2πy cosh u) cos(r u) du
```

and compute a directed-rounding MPFR tensor composite-Simpson sum over `y ∈ [a,b]` (32 subintervals) and `u∈[0,3]` (512 subintervals). The inner integral is exactly the standard integral representation of the real `K_{ir}(2πy)` for real `r` and positive `y`.

The program certifies a **one-sided** interval for each `J_I` by subtracting `10^-10` from its interval Simpson sum, and proves that this conservative subtraction includes all quadrature and truncation errors. The terminal print interval for `J` or the squared mass refers to this **valid lower bound**, not to the exact integral itself.

### Global fourth-derivative inequalities

Put `z=2πy`, `v=cosh u>=1` and `g(y,u)=exp(-z v)cos(r u)/y`. On the entire box `5/3<=y<=21/10`, `0<=u<=3`, we have `z>157/15>10`, `z<66/5<14`, `r<7`, and `exp(-z)<1/30000`. These rational inequalities are independently tested in `audit_cusp_simpson_bounds.py`; for the exponential inequality, the degree-14 Taylor lower sum for `exp(157/15)` exceeds 30000.

The formula for four `u`-derivatives of `exp(-z cosh u)` and Leibniz's rule show

```
|∂_u^4 g| <= (3/5)(1/30000) * P4(14,7)
P4(z,r) = (z⁴+6z³+7z²+z) +4r(z³+3z²+z)
          +6r²(z²+z)+4r³z+r⁴
P4(14,7)=233303.
```

Thus `|∂_u^4 g|<5`. This uses `|sinh u|<=cosh u` and, for `0<=k<=4` and `z>=10`, the monotonicity `exp(-zv)v^k <= exp(-z)` for `v>=1`.

Four direct `y`-derivatives of `exp(-2πyv)/y` give terms
`exp(-2πyv)[(2πv)^4/y +4(2πv)^3/y²+12(2πv)^2/y³+24(2πv)/y⁴+24/y⁵]`. The same monotonicity and `2π<7`, `1/y<=1` give

```
|∂_y^4 ∫_0^3 g(y,u) du|
 < (3/30000)*(7⁴+4*7³+12*7²+24*7+24)
 = 4553/10000 < 1/2.
```

The composite-Simpson error formula `(b-a)*h^4*sup|f''''|/180`, applied in `u` and `y`, with `max(b-a)=1/5`, `h_u=3/512`, `h_y<=(1/5)/32`, implies

```
Error_u < 1.965e-11,  Error_y < 8.48e-13.
```

The neglected `u>3` tail is below `1e-40`: `cosh(3)>10`, `sinh(3)>9`, and convexity of cosh yields `∫_3^∞exp(-z cosh u)du <= exp(-z cosh 3)/(z sinh 3)`. Therefore the **total absolute error in each `J` is strictly below `10^-10`**. All inequalities use exact rational comparisons, validated by `audit_cusp_simpson_bounds.py`. `mpfr_cos`, `mpfr_exp` and `mpfr_log` evaluations are enclosed by directed MPFR rounding; the trigonometric argument interval is widened using the global Lipschitz bound `|cos'|<=1`.

### Replayed bounds

| Quantity | Coarse 512×32 MPFR lower bound | Finer 768×40 MPFR lower bound |
|---|---:|---:|
| J([5/3,7/4]) | 5.42606319029562e-8 | 5.42606318947682e-8 |
| J([7/4,19/10]) | 4.95759796644771e-8 | 4.95759795755669e-8 |
| J([19/10,21/10]) | 2.32621740030071e-8 | 2.32621738528554e-8 |
| Cauchy tail-mass squared bound | **1.91274714714e-13** | **1.91274714324e-13** |

Each displayed decimal is truncated for readability. Full outward-rounded decimal bounds are in the logs. The two runs apply the same conservative error cap, which remains valid for both since both subdivision counts meet the proven minimum. Consequently,

```
mu >= mu_tail > 4e-7
```

holds for the fixed quasimode, assuming the already-established Ford-cusp representation and normalization. The certified squared lower bound actually exceeds `1.9127e-13`, implying `mu_tail>4.37e-7`.

## Analytic residual: independent algebraic cross-check, NOT a certificate

Let `w_g=ψ(y(gP))`, `G=Σ_g w_g`, and on the Ford cell set

```
d_g = f(gP)-f(P)
H = Σ_{g!=1} w_g d_g
W = f + H/G.
```

Every exact defect `d_g` satisfies `(L-λ)d_g=0`. For `L=-y²Δ+y∂_y` and the hyperbolic inner product,

```
(L-λ)W = [Σ_{g!=1}{(Lw_g)d_g -2<∇w_g,∇d_g>}]/G
          - H(LG)/G² -2H|∇G|²/G³ + 2<∇H,∇G>/G².
```

This identity uses no second derivative of the Bessel functions. The independent implementation `verify_picard_defect_residual.py` agrees with the prior analytic `W_residual_cols` evaluator on 42 irregular test points to `2.9e-18` absolute (the maximum sampled residual itself is `8.7e-11`). These checks are **double precision with Chebyshev-interpolated Bessels** and DO NOT bound the residual globally.

The residual vanishes identically for `y>=5/3` because all nonidentity cutoff weights vanish. The only outstanding analytic integration is over the finite curved Ford core.

## Two genuinely different certification targets

- **Coarse existence enclosure:** With the proven `mu>4e-7`, any future validated `R<1e-8` would give `R/mu<0.025`, hence a true odd cuspidal eigenvalue within 0.025 of the fixed `λ=1+r²`. This is NOT a 1e-5 enclosure.
- **Narrow ~1e-5 enclosure:** Prove `mu>1e-5` from a finite collection of compact-core boxes **and** `R<1e-10`. Sampled Gauss quadrature gives `mu_core≈1.908e-5` and `R≈1.2e-11`, but these are midpoint diagnostics. Four disjoint rectangles, `y∈[1,1.4]`, `x∈[-.5,-.25] or [.25,.5]` with `t∈[0,.25]`, and `x∈[-.25,0] or [0,.25]` with `t∈[.25,.5]`, have a combined sampled L2 norm approximately `1.25e-5`. Certifying only their mass, rather than the full Ford core, is enough for the proposed narrow gate.

Neither certified inequality on `R` nor compact-core mass is available yet. **Do not promote spectral theorem flags**, and do not interpret closeness to a Ritz upper bound as a lower bound for the first odd eigenvalue.

## Exact replay

```sh
python audit_cusp_simpson_bounds.py
# On machines without mpfr.h but with libmpfr.so.6 and libgmp installed:
gcc -O2 -std=c11 -o cusp certify_picard_cusp_mass_mpfr.c -l:libmpfr.so.6 -lgmp -lm
./cusp
gcc -O2 -std=c11 -DNT=768 -DNY=40 -o cusp_fine certify_picard_cusp_mass_mpfr.c -l:libmpfr.so.6 -lgmp -lm
./cusp_fine
OPENBLAS_NUM_THREADS=1 python verify_picard_defect_residual.py analytic_quasimode_coefficients_50.json
```

**Next:** implement a validated Bessel/Chebyshev evaluator on the 14-coset compact core, then an adaptive interval Taylor model for the **defect-only** residual identity. Keep the separate three-profile Ritz upper-bound research intact; it does not prove that the quasimode's nearby eigenvalue is first.
