# Track B smooth-gluing audit

Date: 2026-07-13

## Outcome

`rung4_certified=false` remains the only justified state.  The projected-mass
input is GREEN, but the smooth-gluing residual is not presently a theorem
input:

\[
 R\le b_0\delta_0+b_1\delta_1,
 \qquad R<0.0101903405004245
\]

is required for a full lambda width below `0.1`, using the certified
`mu_B > 0.20380681000849094`.

The failure is now localized and fail-closed.  No legacy `eta0` gate is used.

## Dependency status

```text
exact Picard/Humbert geometry
  -> face matrices and two-cusp coordinates                 GREEN
  -> six listed point stabilizers, exact finite enumeration GREEN
  -> epsilon=1/10 curved collar incidence                   RED
  -> edge/vertex transition inventory                       RED
  -> normalized C2 partition and certified b0,b1            RED
  -> six-vector value/gradient overlap boxes                BLOCKED
  -> delta0, delta1 and R                                    BLOCKED
  -> Track-B spectral interval                               BLOCKED
  -> counting/identification                                 independently RED
  -> rung4_certified                                         false
```

## Exact work completed

`track_b_partition_geometry.py` records exact Gaussian-integer face matrices,
the GLUE permutations, the infinity- and zero-cusp coordinate conventions,
and a fail-closed partition ledger.  The five exact codimension-one relations
are `T1`, `T1^-1`, `R`, `Ti R`, and `S`.

The six listed point stabilizers no longer rely on bounded word search.  If
`g=[[a,b],[c,d]]` fixes `(z,y)`, the height denominators for `g` and `g^-1`
give

\[
 |cz+d|^2+|c|^2y^2=1,
 \qquad |-cz+a|^2+|c|^2y^2=1.
\]

At the listed points, `|z|^2<=1/2` and `y^2>=1/2`; hence `|c|^2<=2` and
`|a|,|d|<=2`.  Exhaustive Gaussian enumeration followed by exact rational
fixed-point checks proves stabilizer orders `4,6,6,6,12,12`.  A depth-seven
group-word search independently reproduces every element.

`track_b_overlap_arb.py` implements the downstream calculation without
claiming missing geometry:

- six Acb components are summed before norms are taken;
- transformed values use the stored inverse-right-action permutation;
- gradients use the full PSL2(C) interval Jacobian pullback;
- radial derivatives use `K'_nu=-(K_(nu-1)+K_(nu+1))/2` in Arb;
- finite isotropy averaging and both cusp tags are part of the input schema;
- the weighted residual and threshold comparison use outward Arb arithmetic;
- incomplete geometry produces null defect bounds, never a partial maximum.

`track_b_rung4_integrator.py` verifies matching transition IDs and all theorem
flags before forming the direct `R/mu_B` interval.  It also keeps counting as
a separate hard predicate.

## Remaining mathematical gap

The smallest missing geometric statement is an explicit collar theorem:

> For the epsilon `1/10` thickening of the true curved Humbert cell below
> `y=6/5`, enumerate exactly every incident chamber and singular stratum,
> prove that no other translate meets the collar, and list all edge/vertex
> transition words.

The current 24-simplex FEM mesh cannot prove this statement.  It is an inner
approximation `K_h subset K` and omits the curved sliver `K\K_h`.  Refining
that inner mesh does not by itself become a proof of coverage of `K`.

After the collar theorem, the normalized rational/quintic partition weights
still need interval integration to certify `b0=||B0||_2` and `b1=||B1||_2`.

## Quantitative audit of the proposed partition

The current product-gate majorants are conditional, not theorem inputs:

\[
 b_0<54744.97874481373,
 \qquad b_1<415.3321935811511.
\]

They imply the following necessary targets for closing the *corresponding
supremum-bound proof*:

| allocation of the entire residual budget | required defect |
|---|---:|
| value only | `delta0 < 1.86142012182074e-7` |
| gradient only | `delta1 < 2.45353975875540e-5` |
| equal split, value | `delta0 < 9.3071006091037e-8` |
| equal split, gradient | `delta1 < 1.22676987937770e-5` |

The corrected six-copy collocation diagnostic reports a largest sampled
component defect `7.260191139500237e-5` on `S` (with `R` and `TiR` about
`2.515e-5`).  This is not a certified continuum lower bound, but it is about
390 times the all-value budget and about 780 times the equal-split value
budget.  Therefore merely proving the current conditional constants is very
unlikely to close the interval.  This conclusion is about the present
product-gate/supremum formulation, not about Theorem D(K) itself.

## Ranked bottlenecks

1. **Quantitative stiffness of the candidate partition.**  Its coarse `b0`
   forces a value defect near `1e-7`, inconsistent with the present sampled
   scale near `7e-5`.
2. **Curved-collar incidence and transition completeness.**  Without it no
   maximum over overlaps is a theorem maximum.
3. **Certified pointwise weighted integration.**  The sharper theorem input
   `||B0*d0+B1*d1||_2` has not been implemented; it may avoid the very lossy
   separation into global `b_i delta_i`.
4. **Gradient enclosure and refinement.**  The Arb evaluator exists, but has
   no certified complete box list on which to run.
5. **Counting/identification.**  Even a GREEN Track-B interval does not by
   itself prove the remaining rung-4 counting predicate.

## Smallest credible next changes

1. Import or independently prove the exact Picard cell incidence inventory,
   then certify the epsilon collar by rational/Arb inequalities.  This closes
   a logical gap but does not address the observed quantitative mismatch.
2. Stay within Lemma `Certified automorphization`, but certify the direct
   pointwise norm `||B0*d0+B1*d1||_2` instead of the separated supremum bound.
   This changes neither the theorem nor the arithmetic guarantee.
3. Redesign the partition with wider, geometry-adapted transitions and
   optimize a *proved* derivative norm before recomputing defects.  The
   current epsilon `1/10` product gates should not be promoted merely because
   they are easy to bound.
4. Only after a diagnostic direct weighted residual is below `0.01019` should
   the complete interval subdivision be paid for.  If it stays above the
   threshold, the present trial/partition pair cannot certify Track B.

There is not enough certified data to attach a defensible numerical
probability to any proposed closure.  The evidence does support a firm
negative statement: the atlas proof alone is insufficient, and there is no
evidence that the current `b0*delta0+b1*delta1` formulation can close the
required gap.

## Reproduction

```powershell
python -m py_compile track_b_partition_geometry.py track_b_overlap_arb.py track_b_rung4_integrator.py
python track_b_partition_geometry.py --depth 7
python track_b_overlap_arb.py
python -m unittest test_track_b_integrator.py test_verified_hejhal.py
python track_b_rung4_integrator.py
```

The final two certificate programs intentionally return a blocked/false
result while their inputs are RED.
