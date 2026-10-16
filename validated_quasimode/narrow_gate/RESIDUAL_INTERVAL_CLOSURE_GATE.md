# Picard 40-orbit quasimode: audited narrow-enclosure gates

**Date:** 2026-10-15  
**Repository:** `Sigmoidd/bianchi-selberg`, draft PR #14 (`research/level-one-analytic-quasimode`).  
**Proof state:** **PARTIAL.** New MPFR seed-strip mass inequality is proved modulo the standard Ford-domain and Fourier conventions. The two pointwise global supremum inequalities below are **not certified**. Consequently **no two-sided spectral enclosure has been established by this work**.

## Target and exact trial

`L = -y²(∂xx+∂tt+∂yy)+y∂y`, `r=6.622119348586011` **exact decimal rational**, `λ=1+r²`, and finite 40-orbit seed

\[
f(x,t,y)=\sum_{\mu\ne0}a_\mu yK_{ir}(2\pi|\mu|y)e^{2\pi i(\mu_1x+\mu_2t)},\qquad a_{i\mu}=-a_\mu.
\]

For proof purposes **use the new exactly real and reflection-even projection** in `analytic_quasimode_coefficients_50_reflection_even.json`, SHA256 `f306a3740651575df5be1b5e5130eda956e346b5a3905579b48556eb31d78f98`. It is obtained by taking the real part of the frozen original seed and averaging with its image under `T:(x,t,y)↦(x,-t,y)`. Both operations commute with `L`; `T` normalizes `PSL₂(Z[i])` and conjugates the quarter-turn `Q` to `Q⁻¹`, so the projected trial remains automorphic after Poincaré averaging, exactly `Q`-odd and cuspidal. The coefficient of `(1,0)` remains **exactly 1**.

**Critical correction:** Fourier modes are NOT automatically orthogonal on the *half-square* `x∈[-1/2,1/2],t∈[0,1/2]` for a generic Q-odd coefficient list. The original list is only approximately reflection even. The exact reflection projection is needed for the simple half-square Parseval statement below. The prior cusp-mass argument must be stated for this projected candidate or be supplemented by a separate bound on the reflection-odd part.

The weight is the original nonnegative quintic `ψ=0` at height `h≤3/5`, `ψ=1` at `h≥7/10`, with `G(P)=Σ_γ ψ(y(γP))≥1` on Ford F and `W(P)=Σ_γψ(y(γP))f(γP)/G(P)`.

## New exact coset lemma: 9 rather than 14

Five of the 13 nonidentity representations in the 14-representative implementation have zero cutoff weight **everywhere** on Ford F. Labeling indices 0 through 13 in `representatives()`, remove `{3,6,9,11,12}`. For every nonidentity term, AM–GM gives

\[
h_\gamma=\frac{y}{|cz+d|^2+|c|^2y^2}\le\frac{1}{2|c|\,|cz+d|}.
\]

- For `c=1,d=i`, `|z+d|²≥1`, hence `h≤1/2<3/5`.
- For `c=1,d=±1+i`, `|z+d|²≥5/4`, hence `h≤1/√5<3/5`.
- For `c=1+i,d=-1` or `d=i`, completing squares on the Ford rectangle gives `|cz+d|²≥1/2` and `|c|²=2`; hence `h≤1/2<3/5`.

This **exactly** reduces the Ford-cell representation to the identity plus eight nonidentity cosets; test checks the old/new analytic evaluators agree identically on representative sample points. Vertical cutoff boundaries satisfy the exact **quadratic** equation

\[
a|c|^2y^2-y+a|cz+d|^2=0,\qquad a\in\{3/5,7/10\}.
\]

The root-split diagnostic integrates these smooth vertical branches, but is **not** itself interval-certified.

## New MPFR seed-strip lower bound (completed)

For `1≤y≤5/3`, the Ford domain contains the entire horizontal half-square. Because the **projected** `f` is exactly even in `t`, the half-square norm is one half the full periodic-torus norm and Fourier Parseval applies. Keeping only the normalized orbit of norm one,

\[
\|f\|^2_{\mathrm{strip}}\ge2\int_1^{5/3}\frac{K_{ir}(2\pi y)^2}{y}\,dy.
\]

Divide `[1,5/3]` into eight slabs `[ (12+k)/12,(13+k)/12 ]`, `k=0,…,7`. Cauchy–Schwarz on each slab `I=[a,b]` yields

\[
2\int_I\frac{K_{ir}(2\pi y)^2}{y}\,dy\ge
\frac{2}{\log(b/a)}\left(\int_I\frac{K_{ir}(2\pi y)}{y}\,dy\right)^2.
\]

`certify_picard_seed_slab_mpfr.c` independently encloses the Bessel integral representation using 192-bit **outward-rounded MPFR**, composite Simpson `NT=1024` in `u∈[0,3]` and `NY=32` in `y` per slab. `audit_seed_slab_simpson_bounds.py` uses exact rational arithmetic to check a strict absolute error below `10⁻¹⁰` **for each slab**:

- `|∂ᵤ⁴ exp(-2πy cosh u)cos(ru)/y| < 320`, using `2πy>6`, `r<7`, `exp(-6)<1/400`, and monotonicity of `exp(-z v)v^k` for `v≥1`, `k≤4`.
- `|∂ᵧ⁴ ∫₀³ exp(-2πy cosh u)cos(ru)/y du|<35`.
- Composite Simpson errors sum to less than `3.35×10⁻¹¹`, with `u>3` tail below `10⁻²⁰`.

Two completed replays:

| Run | Rigorous Cauchy lower bound for `||f||²_strip` |
|---|---:|
| 1024×32 MPFR | `>1.7682998477673873e-10` |
| 1536×40 MPFR | `>1.7682998481138830e-10` |

Hence for the exactly reflection-even seed, **`||f||_strip>1.3e-5`**. Log lines headed `certified J` are intervals for **the computed one-sided lower bound**, NOT tight enclosures of the exact integrals.

## Exact reduction of both missing bounds to pointwise gates

The identity cutoff equals one throughout the strip. Writing `Dγ=f(γP)-f(P)`, nonnegative weights imply

\[
W-f=\sum_{\gamma\ne1}\frac{\psi_\gamma}{G}D_\gamma,
\quad\sum_{\gamma\ne1}\frac{\psi_\gamma}{G}\le1.
\]

**Missing mass transfer gate:** verify by directed interval arithmetic

\[
\sup_{P\in F,\,1\le y\le5/3}
\max_{\gamma\ne1:\psi_\gamma(P)>0}|D_\gamma(P)|<7\cdot10^{-6}.
\tag{D}
\]

The volume of the full strip is **exactly `4/25`**, with square root `2/5`. The uniform defect gate would yield

\[
\|W\|\ge\|f\|_\mathrm{strip}-\|W-f\|_\mathrm{strip} > 1.3\cdot10^{-5}-\tfrac25(7\cdot10^{-6})=1.02\cdot10^{-5}.
\]

**Missing residual gate:** enclose `R(P)=(L-λ)W(P)` using the analytic, defect-preserving commutator with **nine active cosets** and interval Taylor models:

\[
\sup_{P\in F,\,y\le5/3}|R(P)|<2\cdot10^{-10}.
\tag{R}
\]

For `u=x²+t²∈[0,1/2]`, the exact inequality `1/(1-u)≤1+2u` gives an **elementary upper bound on the core volume**

\[
\mathrm{vol}(F_{\rm core})=\int_{-.5}^{.5}\int_0^{.5}\left[\frac1{2(1-u)}-\frac9{50}\right]dt\,dx
\le\frac{73}{300}.
\]

Thus (R) implies `||R||₂<√(73/300) × 2e-10 < 1e-10`. Together (D)+(R) yield

\[
\boxed{\|R\|_2/\|W\|_2<10^{-5}},
\]

hence a point of the **Q-odd cuspidal spectrum** lies within `10⁻⁵` of the exact `λ=1+r²`. This does **not** prove it is the first odd eigenvalue, does not prove the global lowest eigenvalue, and does not close any congruence-level or trace-formula exclusion.

## Numerical diagnostics, explicitly non-certified

- New scalar-specialized analytic evaluator agrees with the original 40-column analytic evaluator to `~6.3e-19` at the tested points.
- Root-split Gauss grids `7×6×5`, `11×9×7`, and `15×12×9` evaluate the **original** 40-orbit trial with `||W||_core≈1.90782e-5` and `||R||≈1.14–1.19e-11`. The last two residual estimates are **not monotone**, so do not infer convergence bounds.
- After exact real/t-even projection, the `11×9×7` root-split evaluation changes `||R||` by less than `1e-18`.
- Independent Sobol screens (16,384 core points and 16,384 strip points) for the **projected** trial observe `max |R|≈1.1564e-10`, `max active |Dγ|≈1.695e-12`, and `max |W-f|≈4.38e-13`. They provide **no upper bound on a supremum**.

## Remaining implementation task and stop condition

**Still OPEN:** implement a verified Bessel/Chebyshev radial table, including derivative error bounds, then an adaptive interval branch-and-bound for each `Dγ` on the strip, followed by an interval Taylor-model enclosure of the analytic commutator `R` on the nine-coset Ford core. Use the exact quadratic cutoff surfaces to keep each interval cell analytic; verify no uncovered gaps or ambiguous endpoint cases. Record the maximum outward-rounded interval bound across **every** cell and prove termination. For (R), raw interval evaluation of separately large terms will likely overestimate because of cancellation: use the defect-only identity and Taylor cancellation before applying interval hulls.

**Do not promote spectral flags or claim a `10⁻⁵` enclosure unless both (D) and (R) have passed.** The new seed-strip mass bound is genuine progress; it is not, by itself, a completed spectral certificate.

### Minimal replay

```sh
OPENBLAS_NUM_THREADS=1 python reflection_even_trial.py
python ford_active_cosets.py
python proof_gate_audit.py
python audit_seed_slab_simpson_bounds.py
gcc -O2 -std=c11 -DNT=1024 -DNY=32 -o certify_seed_mass certify_picard_seed_slab_mpfr.c -l:libmpfr.so.6 -lgmp -lm
./certify_seed_mass
# Independent refinement:
gcc -O2 -std=c11 -DNT=1536 -DNY=40 -o certify_seed_mass_fine certify_picard_seed_slab_mpfr.c -l:libmpfr.so.6 -lgmp -lm
./certify_seed_mass_fine
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_residual_support_reduction
OPENBLAS_NUM_THREADS=1 python root_split_diagnostics.py
OPENBLAS_NUM_THREADS=1 python sample_uniform_gates.py
```
