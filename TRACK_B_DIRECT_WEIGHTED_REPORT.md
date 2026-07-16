# Track B direct weighted-residual report

Date: 2026-07-13

## Result

The separated global bound `b0*delta0+b1*delta1` was replaced by an Arb
implementation that preserves the local correlation in

\[
\|B_0d_0+B_1d_1\|_2
\]

and, where possible, evaluates the still sharper cancellation-preserving
Laplacian commutator itself.

This exposed a certified obstruction in the original cusp/core construction
and a substantially better replacement construction.

## Original construction: certified obstruction

The original partition uses `W_B` only in the cusp and blends it with the raw
six-copy field `F` on `1.01 <= y <= 1.20`.  Exact rational-frequency Fourier
Gram integration collapses 128 frequencies to 20 radial Bessel shells.

On the central region

```text
x1 in [-0.39,0.39]
x2 in [ 0.11,0.39]
y  in [ 1.06,1.20]
```

all five core gates are identically one.  There are no face, edge, or floor
partition terms there.  The residual is exactly the cusp commutator

\[
(\Delta\chi_B)(W_B-F)-2\langle\nabla\chi_B,\nabla(W_B-F)\rangle.
\]

A 192-bit, 512-segment second-order Taylor/Whittaker-ODE calculation proves

\[
0.012400710256553562
 < \|R_{\rm cusp}\|_{L^2(\text{central region})}
 <0.01923144927660803.
\]

The complete width budget is

\[
R<0.0101903405004245.
\]

Thus the fixed trial with the original `1.01 -> 1.20` cusp/core blend cannot
close width `<0.1`, even if the missing Humbert atlas were supplied.  This is
a lower bound on a subdomain, so residual contributions elsewhere cannot
cancel it in the global L2 norm.

The non-load-bearing full half-torus midpoint diagnostic is consistent:

```text
requested pointwise weighted norm  ~ 0.03792416
actual cusp commutator norm         ~ 0.02575275
```

Widening the cusp cutoff is not attractive.  At plateau `1.30` the diagnostic
commutator falls to `0.01367`, but the certified witness mass on
`[1.30,1.55]` falls to `0.02344`; the corresponding `0.05*mu` budget is only
`0.001172`.  The estimated residual/budget ratio worsens.

## Replacement: use the projected field in every chart

Define every local core candidate to be the same exact target-projected old
field used in the cusp:

\[
W_B=P_{Q,-}P_{-z,+}P_{\rm old}F.
\]

This finite Fourier field still solves the eigen-equation termwise.  It is
integer-periodic, `-z` even, quarter-turn odd, and repeated in all six old
fibers.  Therefore:

- the cusp/core defect is exactly zero;
- `T1`, `R`, and `TiR` defects vanish algebraically;
- the only nonzero fundamental Picard transition is the inversion `S`.

Independent Arb point diagnostics give:

| relation | six-vector value defect | six-vector first-gradient defect |
|---|---:|---:|
| `T1` | `< 4.0e-51` | `< 8.1e-50` |
| `R` | `< 4.0e-51` | `< 8.1e-50` |
| `TiR` | `< 4.0e-51` | `< 8.1e-50` |
| `S` | `< 8.625e-7` sampled | `< 1.704e-3` sampled |

For the candidate floor gate

\[
q=10u^3-15u^4+6u^5,
\qquad u=\log(x_1^2+x_2^2+y^2)/0.1,
\]

the exact identities

\[
\Delta\log(x_1^2+x_2^2+y^2)=0,
\qquad
|\nabla\log(x_1^2+x_2^2+y^2)|_{\mathbb H}
=\frac{2y}{\sqrt{x_1^2+x_2^2+y^2}}
\]

give a direct floor-shell integrand.  Arb midpoint quadrature in exact
`(x1,x2,log(r^2))` coordinates currently estimates

```text
|| |Delta q| d0 + 2 |grad q| d1 ||_2  ~ 0.01239640
actual cancellation-preserving commutator ~ 0.00810142
target                                      0.01019034
```

The requested triangle bound is now only about 22 percent above target, and
the actual commutator is already below target.  These floor numbers are
diagnostics, not yet continuum certificates.

## Next load-bearing calculation

The floor-width/refinement diagnostics now give:

| width | grid | requested norm | actual commutator |
|---:|---:|---:|---:|
| `0.15` | `8x4x8` | `0.00916188` | `0.00539611` |
| `0.15` | `12x6x12` | `0.01276456` | `0.00808610` |
| `0.20` | `12x6x12` | `0.01016525` | `0.00589625` |
| `0.25` | `12x6x12` | `0.00856394` | `0.00452729` |
| `0.25` | `16x8x16` | `0.00944768` | `0.00516263` |
| `0.25` | `20x10x20` | `0.00985829` | `0.00546536` |
| `0.30` | `16x8x16` | `0.00845702` | `0.00429240` |

Widths `0.15` and `0.20` do not provide credible interval margin.  Width
`0.25` appears to converge too close to the threshold.  Width `0.30` is the
current load-bearing candidate, with a diagnostic requested-bound ratio
`0.8299` and commutator ratio `0.4212`.  Unlike widening the cusp cutoff,
widening the floor gate does not reduce projected mass because the same
`W_B` field is used on both cusp and core charts.

Two rigorous continuum attempts were also run at width `0.30`:

```text
direct independent Arb boxes, 8x4x8     upper ~ 2.03e10
legacy second-order jet / interval Hessian upper ~ 3.29e12
correlated (x1,x2,log(r^2)) jet, 8x4x8     upper ~ 2.55e8
```

The corrected jet uses

\[
y^2=e^s-x_1^2-x_2^2,
\qquad
S(x_1,x_2,y)=(-e^{-s}x_1,e^{-s}x_2,e^{-s}y),
\]

and evaluates the commutator as

\[
-4y^2e^{-s}\{q_{ss}D+q_s(x_1D_{x_1}+x_2D_{x_2}+2D_s)\}.
\]

It also folds the exact quarter-turn-odd four-mode orbits, pairs those
orbits under the floor reflection, and sums equal radial shells before
applying interval absolute values.  Its value at the first cell midpoint
agrees with both the direct action/Jacobian implementation and the legacy
jet, with difference balls below `1e-40`.  The new certified upper is about
`6.5e3` times smaller than the legacy jet bound, but remains useless for the
theorem.

The remaining inflation is now localized.  The largest projected Fourier
coefficient has modulus about `1.7447e4`, while representative midpoint
commutator derivatives are only `1e-2` to `5e-2`.  On the same cell, the
interval-gradient enclosure is `1e8`.  All 56 Bessel evaluations in this
diagnostic were direct Arb evaluations, with no mean-value fallback.  Thus
the present obstruction is cancellation among the large Fourier terms in
the interval Hessian, not Bessel evaluation, shell geometry, or the
`(x_1,x_2,s)` coordinate formula.

Center-affine diagnostics quantify what a genuine Taylor model must prove:

| grid | affine cell-sup envelope | affine-polynomial L2 part |
|---:|---:|---:|
| `8x4x8` | `0.01187715` | `0.00528463` |
| `16x8x16` | `0.00980522` | not run in this audit |

The affine quantities are not certificates because the second-order
remainder is omitted.  The `8x4x8` cell-sup route cannot close even with a
zero remainder.  In contrast, exact cell integration of the affine
polynomial leaves about `0.00491` of L2 norm budget for a rigorous remainder.
This makes a cancellation-preserving Taylor model with exact polynomial L2
integration the smallest credible provable improvement.  Plain interval
Hessians or geometric box refinement are not credible closure mechanisms.

The remaining proof tasks are:

1. represent the already-subtracted scalar `S` residual by a validated
   Taylor/Chebyshev model on each floor cell, summing its polynomial
   coefficients across Fourier modes before bounding the remainder;
2. integrate the squared polynomial part cellwise and prove an L2 bound for
   the Taylor remainder, rather than replacing each cell by a supremum;
3. feed the certified direct residual into `track_b_rung4_integrator.py`.

The finite floor-collar incidence is now certified in
`track_b_partition_geometry.py`: `exp(0.30) < 27/20` is checked by a rational
Taylor-tail bound; `|c|^2>2` is excluded from the closed Picard cell by
`y^2>=1/2`; the `|c|^2=2` cases reduce to the two lower vertices; and all nine
unit-`c` translated `S` spheres are enumerated exactly.  Only `S` enters the
open collar.  Its four closed-collar elliptic edges use cyclic averages of
orders `3,3,2,3`, while vertex averages use the existing exact stabilizers.

The projected-core construction is materially closer to dual certification
than the original five-gate/supremum construction.  It removes the certified
cusp obstruction and reduces the numerical residual problem to one scalar
transition.

## Reproduction

```powershell
python track_b_direct_weighted_arb.py `
  --bits 192 `
  --lower-bound-segments 512 `
  --floor-widths 0.30 `
  --floor-subdivision 16,8,16 `
  --taylor-floor-width 0.30 `
  --taylor-floor-subdivision 8,4,8 `
  --affine-floor-width 0.30 `
  --affine-floor-subdivision 8,4,8
```

`rung4_certified=false` remains mandatory until the floor-shell midpoint
diagnostic is replaced by a complete interval integral.  The finite collar
incidence ledger is GREEN, but the broader normalized partition remains
fail-closed pending its global collar and weight proofs.

## July 2026 global-partition update

The two limitations in the preceding paragraph have since been separated.
The degree-10 floor Taylor integral is certified and stable, and the complete
global normalized partition is now independently verified.  The latter has
exact denominator `Phi=1`, exact partition deviation zero, finite certified
first-gradient/Laplacian bounds on all 2,048 primary cells, a complete global
elliptic-stratum ledger, zero fallback, and zero change under the
`32x16x32` stability refinement.  See
`TRACK_B_GLOBAL_PARTITION_CERTIFICATE.md` and
`track_b_global_partition_result.json`.

`rung4_certified=false` still remains mandatory.  The current overlap
artifact is RED because the complete two-cusp non-floor value/first-gradient
defects have not yet been regenerated and certified against the new
partition identifiers.  The global partition result does not by itself
make an eigenvalue-existence claim.
