# Track B global normalized-partition certificate

## Result

The normalized partition and its theorem-facing first-gradient and
hyperbolic-Laplacian bounds are certified at 192 Arb bits on the complete
closed Humbert reference cell.

The primary run uses `16x8x16` cells in `(x1,x2,s=log rho)`.  The stability
run uses `32x16x32`.  Both prove

```text
minimum denominator lower       1
maximum |grad chi_j|_H upper     114.3421052631578947368421052631579...
maximum |Delta chi_j| upper      14726.475530932594644506001846722...
partition deviation upper       0
unresolved fallback count       0
```

The refinement changes the denominator, gradient maximum, and Laplacian
maximum by exactly zero at the stored precision.  The independent verifier
reconstructs all three extrema from 2,048 primary audit leaves and confirms
the exact-once cell set and all deterministic hashes.

The resulting flags are

```text
floor_residual_certified          true
global_partition_certified        true
global_weight_bounds_certified    true
rung4_certified                   false
```

The last flag remains false.  This certificate is a partition/weight input,
not an eigenvalue-existence theorem.

## Partition construction

Let `H` be the clamped quintic

```text
H(t) = 0                              t <= 0
     = 10 t^3 - 15 t^4 + 6 t^5       0 < t < 1
     = 1                              t >= 1.
```

There is one cusp gate and five complementary core gates.  The 32 core raw
weights are the tensor products of either `a_i` or `1-a_i`, multiplied by
the cusp complement.  Hence the raw denominator is not bounded from samples;
it is simplified in the exact coefficient ring:

```text
Phi = cusp + (1-cusp) product_i (a_i + (1-a_i)) = 1.
```

Consequently `Phi_lower=Phi_upper=1`, `grad Phi=0`, `Delta Phi=0`, and the
normalized partition identity has exact deviation zero.

The non-plateau domain is covered because `y<=6/5` implies
`rho<=1/2+(6/5)^2=1.94`, while Arb proves `exp(0.70)>1.94`.  The remaining
`y>=1.20` region is the analytic cusp plateau with its single constant
weight.  Grid interiors use half-open ownership; closed interfaces are kept
only in the incidence ledger.

## Metric and derivatives

Cells are rectangular in `(x1,x2,s)`, but every stored derivative is in the
original `(x1,x2,y)` coordinates.  The conversion for `s=log rho` is
implemented through its exact gradient and Hessian.  Norms use

```text
|grad f|_H^2 = y^2 (f_x1^2 + f_x2^2 + f_y^2)
Delta f = -y^2(f_x1x1+f_x2x2+f_yy) + y f_y.
```

For the project's positive Laplacian convention, automatic jets and the
explicit quotient rule agree:

```text
Delta(phi/Phi)
 = Delta(phi)/Phi - phi Delta(Phi)/Phi^2
   + 2 <grad phi,grad Phi>_H/Phi^2
   - 2 phi |grad Phi|_H^2/Phi^3.
```

The final sign is negative for this convention.  A regression test compares
the explicit gradient and Laplacian formulas with second-order automatic
differentiation at an Arb point.

The certified derivative bounds use exact factorwise product rules with
`max |H'|=15/8` and `max |H''|<=6`.  This avoids interval dependency across
six factors.  Direct interval jets are retained as diagnostics.  If such a
diagnostic becomes unbounded on a coarse cell, it is explicitly recorded;
it is never substituted for the analytic majorant and is not counted as a
fallback.

The separated `b0,b1` constants deliberately sum the global majorant over
all 33 labels, because Reynolds averaging can activate a label that vanishes
in the canonical representative.  This is conservative:

```text
B0 sup <= 471498.39427516158818...
B1 sup <=   7341.57894736842105...
b0     <= 333399.71191053519708...
b1     <=   5191.28025830060600...
```

## Exact singular-stratum averaging

The global geometry dependency now classifies all positive-height singular
strata of the five-face Humbert cell:

```text
3 self-paired-face fixed axes       orders 2,2,2
4 vertical edge cycles              orders 2,2,2,2
4 floor edge cycles                 orders 3,3,2,3
6 stabilizer-jump vertices          orders 4,6,6,6,12,12
```

Generic stabilizers and all vertex groups are enumerated with Gaussian
integer/rational arithmetic.  For every serialized group, exact matrix
closure and right-multiplication reindexing are rechecked.  The common
Reynolds operator

```text
A_G phi = |G|^-1 sum_(g in G) phi o g
```

is therefore invariant by exact reindexing.  Isometry covariance gives the
corresponding first-gradient transformation and Laplacian invariance.  Group
terms do not become integration leaves, so the half-open reference-cell
ledger counts each positive-measure region once.

## Floor compatibility

The aggregate global floor-core weight simplifies to the same
`H(log(rho)/0.30)` used by the local Taylor certificate.  Both artifacts now
carry

```text
floor_weight_formula_id = track-b-floor-quintic-logrho-width-0.30/v1
floor_geometry_incidence_hash = aec0a1a0...f87d540
```

The integrator rejects different widths, formula identifiers, or geometry
hashes.  It also rejects provisional results, missing ledgers, byte-hash
changes, configuration/definition/formula mismatches, midpoint-only bounds,
or missing per-weight intervals.

## Reproduction

```powershell
& 'C:\Users\Admin\AppData\Local\Programs\Python\Python313\python.exe' `
  track_b_partition_geometry.py `
  --depth 6 `
  --certify-global-partition `
  --bits 192 `
  --partition-subdivision 16,8,16 `
  --partition-degree 8 `
  --partition-workers 4 `
  --partition-audit-jsonl track_b_partition_cells.jsonl `
  --check-refinement 32,16,32 `
  --check-audit-jsonl track_b_partition_cells_refined.jsonl `
  --global-json-out track_b_global_partition_result.json

& 'C:\Users\Admin\AppData\Local\Programs\Python\Python313\python.exe' `
  track_b_global_partition_verify.py `
  --result track_b_global_partition_result.json `
  --audit-jsonl track_b_partition_cells.jsonl `
  --json-out track_b_global_partition_verification.json
```

## Remaining residual channels

The current `track_b_overlap_result.json` is still RED.  It has not yet
certified the complete two-cusp overlap defect transport, stabilizer-averaged
values, or first gradients against the new partition artifact.  Therefore
the remaining work is:

1. regenerate the complete exact two-cusp/face overlap ledger against the
   new partition identifiers;
2. certify the non-floor face, cusp-transition, and tail defect terms and
   their first gradients;
3. combine those channels with the closed floor term in the interval
   integrator;
4. perform the final width/counting comparison required by Rung 4.

No defect-to-spectrum `kappa` theorem and no Hejhal eigenvalue-existence
claim is made here.
