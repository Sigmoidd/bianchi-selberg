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

The script now supports a floor-width sweep.  Width `0.15` is the first
candidate to test rigorously: simple derivative scaling from the width `0.10`
diagnostic predicts a requested norm near `0.0097`, but this prediction is
not a certificate.  Unlike widening the cusp cutoff, widening the floor gate
does not reduce the projected mass because the same `W_B` field is used on
both cusp and core charts.

The remaining proof tasks are:

1. run the `0.10,0.15,0.20` floor sweep and select a width with diagnostic
   margin;
2. interval-subdivide the floor shell and certify its direct commutator norm;
3. prove that the selected floor collar meets only the required incident
   chambers, treating elliptic edge averages explicitly;
4. feed the certified direct residual into `track_b_rung4_integrator.py`.

The projected-core construction is materially closer to dual certification
than the original five-gate/supremum construction.  It removes the certified
cusp obstruction and reduces the numerical residual problem to one scalar
transition.

## Reproduction

```powershell
python track_b_direct_weighted_arb.py `
  --bits 192 `
  --lower-bound-segments 512 `
  --floor-widths 0.10,0.15,0.20 `
  --floor-subdivision 8,4,8
```

`rung4_certified=false` remains mandatory until the floor-shell midpoint
diagnostic is replaced by a complete interval integral and the collar
incidence proof is GREEN.
