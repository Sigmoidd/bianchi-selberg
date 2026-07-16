# Track B floor Taylor certificate

## Status

The validated continuum model is implemented.  Global
`rung4_certified=false` remains fail-closed because the normalized global
partition and global weight bounds are not certified.  The first complete
degree-8, `8x4x8` local floor run is also RED because its rigorous Taylor
remainder is too large.

## Validated cell model

Each cell uses centered normalized variables

\[
x_1=c_1+h_1\xi_1,\qquad x_2=c_2+h_2\xi_2,\qquad
s=c_s+h_s\xi_s,\qquad \xi\in[-1,1]^3.
\]

The implementation constructs

\[
G=\sqrt{6J}\,R_S=P+E,\qquad |E|\le\varepsilon_C,
\]

where

\[
J=\frac{e^s}{2(e^s-x_1^2-x_2^2)^2}.
\]

The factor `6` is the exact squared norm of the scalar field repeated in all
six fibers.  Taylor polynomials use outward-rounded complex Arb
coefficients.  Products first accumulate coefficients by multi-index; only
discarded coefficients and final remainder terms are replaced by absolute
upper bounds.

Before this model is bounded, the code:

1. uses the existing old-field projection;
2. folds the exact quarter-turn-odd four-mode orbits;
3. pairs those orbits under floor reflection;
4. sums equal radial Bessel shells;
5. sums every retained multivariate polynomial coefficient.

The radial Taylor coefficients are generated from

\[
y^2r''-yr'+(1+r_0^2)r-\omega^2y^2r=0
\]

using direct Arb values of `K_(ir)` and its first derivative at the cell
center.  Taylor tails use

\[
|K_{ir}^{(n)}(x)|
\le 2^{-n}\sum_{j=0}^n {n\choose j}K_{|-n+2j|}(x),
\]

with real-order Arb Bessel values at the lower endpoint.  These are reported
as analytic majorants, not hidden fallbacks.  A certificate requires the
fallback count to be zero.

## Polynomial integration and remainder

The total-degree monomial polynomial is converted exactly to a tensor
Legendre basis.  Its cell integral is bounded by

\[
A_C=h_1h_2h_s\sum_{a,b,c}|p_{abc}|^2
\frac{8}{(2a+1)(2b+1)(2c+1)}.
\]

For `V_C=8h_1h_2h_s`, the cell contribution is

\[
U_C=\sqrt{A_C}+\varepsilon_C\sqrt{V_C},
\qquad
U_{\rm floor}=\left(\sum_C U_C^2\right)^{1/2}.
\]

The allowed residual budget is not hard-coded.  It is reconstructed from
the certified projected-mass lower bound and requested spectral full-width:

\[
R_{\rm allowed}=\frac{\mu_B\,(0.1)}2.
\]

The only passing comparison is

```text
upper(U_floor) < lower(R_allowed).
```

## Complete degree comparison

Both rows are complete 192-bit Arb runs on width `0.30` and grid `8x4x8`
(256 cells).

| total degree | polynomial L2 upper | remainder L2 upper | floor L2 upper | ratio to budget |
|---:|---:|---:|---:|---:|
| 8 | `0.004952511806261643` | `0.030045824906546633` | `0.032783307947250443` | `3.217096421` |
| 10 | `0.004952510697927052` | `0.000862241231152792` | `0.005410378421602105` | `0.530932055` |

The certified budget lower endpoint is
`0.010190340500424547`.  Degree 10 therefore has certified local margin
`0.004779962078822442` and meets the requested robust target.  It remains
provisional until the separate stability test below succeeds.

The fixed-grid degree contraction is

\[
\rho_p=E_{10}/E_8=0.0286975389704.
\]

This is **strong degree convergence** under the prescribed classification,
and it passes the development milestone by a factor of about 26.1.  The
worst cell remains `[0,3,0]`; its remainder L2 upper falls from
`0.0125628534354` to `0.000368634361511`.  Thus the global reduction is not
caused by the identity of the worst cell changing.

Approximate observed costs on four workers were 1,555 seconds for degree 8
and 7,267 seconds for degree 10.  Worker peaks were about 43--45 MiB, with a
47 MiB parent, so aggregate peak resident memory was approximately
220--227 MiB.  A uniform `16x8x16` degree-8 run has eight times as many cells;
the measured cell cost predicts about 3.45 hours at the same worker count and
similar peak memory.  The rigorous targeted replacement below was used
instead.

## Rigorous refinement comparison

The seven cells comprising the smallest prefix exceeding 90 percent of the
degree-8 squared remainder were each replaced by their eight disjoint
children.  The remaining 249 parents were retained exactly once.  The
resulting 305-leaf audit gives

```text
polynomial L2 upper   0.004952508172016106
remainder L2 upper    0.009300297990944841
floor L2 upper        0.012951062301652950
elapsed               346.97 seconds (4 workers)
```

Hence

\[
\rho_h=0.309537116051.
\]

Refinement also passes the `0.75 E8` convergence milestone, but does not by
itself close the local floor inequality.  On this evidence fixed-grid degree
10 is the cheapest credible closer: its contraction is about 10.8 times
stronger than the one-level targeted refinement contraction.

The adaptive driver supports maximum depth, a chosen subset of the axes,
an L2 threshold, top fraction by cell count, top fraction of squared
remainder, and a target global remainder.  Every replacement event records
one parent and all children; duplicate leaf identifiers and non-one-for-one
parent removal are hard errors.

## Degree-10 stability certificate

The same exact replacement was applied to the seven cells contributing over
90 percent of the degree-10 squared remainder.  It completed in 1,346 seconds
on four workers:

| quantity | closing degree 10 | refined degree 10 |
|---|---:|---:|
| polynomial L2 upper | `0.004952510697927052` | `0.004952510590318300` |
| remainder L2 upper | `0.000862241231152792` | `0.000244153678091229` |
| floor L2 upper | `0.005410378421602105` | `0.005097526161762399` |

Both bounds close.  Their difference is
`0.000312852259839707`.  Ten percent of the smaller positive margin is
`0.000477996207882245`, so the prescribed stability inequality holds.

`track_b_floor_stability.py` independently rereads the closing result, base
ledger, adaptive leaf ledger, and parent-child ledger.  It verifies the
strict inequalities, all local analytic and geometry flags, zero Bessel
fallbacks, unique leaves, exact parent removal/child inclusion, the 90-percent
selection condition, and the margin comparison.  Its result is

```text
stability_check_passed     true
floor_residual_certified   true
rung4_certified            false
```

Thus the validated Taylor remainder is no longer a local floor blocker.
The remaining blockers are the explicitly separate global gates.

## Cell concentration and geometry

For both degrees 8 and 10, 3 cells contribute 50 percent and 7 cells
contribute 90 percent of the rigorous global squared remainder.  The exact
top-20 tables, coordinates, cumulative contributions, and closed-collar
incidences are in `track_b_floor_d8_n8_ledger_analysis.json` and
`track_b_floor_d10_n8_ledger_analysis.json`.

The worst cell is

\[
x_1\in[-1/2,-3/8],\quad x_2\in[3/8,1/2],\quad
s\in[0,3/80].
\]

It touches the inversion sphere, elliptic edges `floor_x1m` and
`floor_x2p`, and vertex `v_mhh`.  Its reflected partner `[7,3,0]` touches
`floor_x1p`, `floor_x2p`, and `v_phh`.  The next four dominant cells have
the same horizontal corner geometry at increasing `s`; only the `k=0`
cells actually touch the inversion sphere and elliptic strata.

The audit also found a serialization defect in the old JSONL coordinate
fields: `arb.union(a,b)` is an error-ball union about zero, not a convex
endpoint serialization.  This did not affect any Taylor model or norm,
which used endpoint pairs directly.  New records store `x1_bounds`,
`x2_bounds`, and `s_bounds`; old ledgers are reconstructed exactly from
their complete uniform-grid indices.

## Omitted indices, cancellation, and Bessel majorants

An instrumented recomputation of the worst cell attributes discarded product
coefficients without changing the certified remainder.  At degree 8 the
largest omitted multi-indices are

```text
(4,4,1), (3,5,1), (5,3,1), (3,4,2), (4,3,2),
(2,6,1), (6,2,1), (2,5,2)
```

At degree 10 the leading pattern shifts to

```text
(5,5,1), (4,6,1), (6,4,1), (5,4,2), (4,5,2),
(7,3,1), (3,7,1), (6,3,2)
```

The evidence is horizontal and mixed-horizontal, not `s` dominated.  True
tensor caps `(px,py,ps)` are now supported.  The justified diagnostic order
is therefore `(10,10,8)` followed, only if needed, by `(12,12,8)`; the
`s`-first alternatives are not justified by the omission ledger.

Worst-cell `(10,10,8)` and control `(8,8,10)` computations were launched
independently.  Each exceeded 1,345 CPU seconds (over 22 minutes wall time)
and about 53 MiB without completing, and was stopped.  This is a censored
cost result, not a residual result.  The true tensor composition order is
`px+py+ps=28`, compared with total order 10 in the successful run; projecting
that measured worst-cell cost to all 256 cells makes these tensor settings
strictly less credible as the immediate closure route.  No tensor residual
value is claimed.

The maximum finite coefficient-addition cancellation condition in the worst
cell is `1.1531042e7`, at multi-index `(0,1,6)`, at both degrees.  There are
also additions whose interval result contains zero (1,320 at degree 8 and
2,291 at degree 10); these are reported separately rather than assigned a
spurious finite condition number.  This arithmetic diagnostic is not used
as a proof bound.

Bessel accounting is:

| degree | direct complex-order | real-order majorants | fallback | Fourier tail |
|---:|---:|---:|---:|---:|
| 8 | 7,680 | 28,160 | 0 | 0 |
| 10 | 7,680 | 33,280 | 0 | 0 |

Per cell, degree 8 evaluates each real integer order 0 through 10 ten times;
degree 10 evaluates orders 0 through 12 ten times.  No backend or fallback
change explains the convergence.

## Cell-total formula

The current valid Minkowski bound is retained.  Merely taking
`L_C <= sqrt(A_C V_C)` in

\[
U_C^2\le A_C+2\varepsilon_C L_C+\varepsilon_C^2V_C
\]

reproduces the existing bound exactly and cannot improve it.  A useful
implementation would therefore require a genuinely sharper validated
polynomial `L1` range/subdivision bound.  Since degree 10 already closes with
a 47-percent budget margin, that additional machinery is not load-bearing
and has not been inserted into the theorem path.

## Gate interpretation

The degree-10 inequality is a local floor result only.  The code continues
to force `floor_residual_certified=false` until the stability criterion is
proved.  Independently, all current outputs force `rung4_certified=false`
while any of these remain false:

```text
global_partition_certified
global_weight_bounds_certified
rung4_integrator_comparison_certified
```

## Reproduction

```powershell
& 'C:\Users\Admin\AppData\Local\Programs\Python\Python313\python.exe' `
  track_b_direct_weighted_arb.py `
  --bits 192 `
  --floor-widths 0.30 `
  --certify-floor `
  --floor-subdivision 8,4,8 `
  --floor-model taylor `
  --floor-degree 8 `
  --floor-workers 4 `
  --floor-audit-jsonl track_b_floor_cells_d8_n8.jsonl `
  --json-out track_b_floor_d8_n8_result.json
```

The adaptive comparison is reproduced by

```powershell
& 'C:\Users\Admin\AppData\Local\Programs\Python\Python313\python.exe' `
  track_b_floor_adaptive.py --adaptive-floor `
  --base-audit track_b_floor_cells_d8_n8.jsonl --base-grid 8,4,8 `
  --degree 8 --bits 192 --workers 4 `
  --adaptive-contribution-fraction 0.90 --adaptive-max-depth 1 `
  --adaptive-axes x1,x2,s --adaptive-target-remainder 0.004
```

Even if a local floor configuration becomes GREEN, it remains provisional
until a degree increase or grid refinement satisfies the formal stability
inequality.  It cannot set global `rung4_certified=true` while the global
partition, global weights, or Rung-4 integration gates remain false.
