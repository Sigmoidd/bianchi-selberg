# Track B overlap-defect certificate

Status: **fail closed pending the certified Humbert partition artifact**.

The implementation in `track_b_overlap_arb.py` computes the two overlap
inputs required by the six-copy/two-cusp theorem:

\[
\delta _0=\operatorname*{ess\,sup}\lVert F_j^{(k)}-F_k\rVert_{\mathbf C^6},
\qquad
\delta _1=\operatorname*{ess\,sup}
 \lVert\nabla F_j^{(k)}-\nabla F_k\rVert.
\]

## What is certified by the evaluator

- All six components are evaluated together.
- The common-fiber convention is the stored inverse right action
  `pi_g`: `defect[c] = F_c(p) - F_pi_g(c)(g p)`.
- First derivatives of the transformed field include the exact interval
  Jacobian pullback `J_g^T`.
- The covector norm is the hyperbolic norm, namely `y` times the Euclidean
  coordinate-covector norm.
- Finite elliptic stabilizer averages are applied before overlaps are
  subtracted.
- Infinity-cusp, zero-cusp, and compact transitions are separately tagged.
- Each Whittaker value and each of its three first derivatives is summed as
  an Acb ball before an absolute value is taken.  The previous absolute
  modal Taylor majorant is not used.
- The radial derivative uses the exact recurrence
  `K'_nu=-(K_(nu-1)+K_(nu+1))/2`, with all three Bessel values evaluated by
  Arb.

## Strict geometry contract

The default input is `track_b_humbert_partition_result.json`, schema
`track-b-humbert-partition-v1`.  Execution is rejected unless all of these
are explicitly true:

- `certified`
- `coverage_certified`
- `local_finiteness_certified`
- `transition_set_complete`
- `stabilizers_complete`
- `two_cusp_coordinates_certified`

Every active transition must provide an exact determinant-one Gaussian
matrix, its six-copy permutation, nonempty overlap boxes, a unique ID, a
cusp-coordinate tag, and any source/target stabilizer averages.  Missing
geometry cannot silently turn into a smaller transition maximum.

## Closure criterion

For the already certified Track-B mass

\[
\mu_B>0.20380681000849094,
\]

a final spectral width below `0.1` requires

\[
R=b_0\delta_0+b_1\delta_1<0.0101903405004245.
\]

The canonical result file is `track_b_overlap_result.json`.  It records raw
and weighted defect bounds, the complete transition-ID list, per-transition
maxima, refinement monotonicity, independent Jacobian checks, theorem
compatibility, and all source hashes.

## Current blocker

No numerical `delta0_upper` or `delta1_upper` is claimed until the partition
producer supplies a certified and complete active-overlap list.  Therefore
`rung4_certified=false` remains mandatory at this stage.  Once that artifact
is present, the evaluator can run without any theorem or guarantee change.
