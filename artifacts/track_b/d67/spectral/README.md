# d=67 threshold route checkpoint

**Spectral exclusion remains OPEN.** The adaptive face-fragment moment map
and tuned exact scalar budget now pass. No interval threshold matrices, verified
positivity certificate, or independent final spectral replay exist.

`docs/D67_TRACK_B_THRESHOLD_THEORY.md` proves a conditional exhaustive
criterion: a single negative-index check at lambda=1 suffices to exclude
all discrete eigenvalues in (0,1). Nonzero cusp modes supply a positive
Robin bound. A rational mapped-core coefficient formulation avoids the
historical erroneous floor-lift lemma and its slivers.

## Closed evidence and scope

- The exact reference mesh has 228 triangles. Its independent verifier
  checks full floor coverage/nonoverlap, all affine Ford triangle images,
  metric and height identities, and periodic boundary subdivisions.
- Normal and optimized reference replay pass. Eight sensitivity controls
  reject altered field/group identities, unearned theorem claims, hashes,
  missing triangles/periodic edges, false pairings and displaced vertices.
- Floating CR probes (diagnostics only) give minimum pencil values
  1.452312794 at 8,322 dofs and 1.515676297 at 66,120 dofs, for mean penalty
  3. The weaker penalty 1/2 gives 1.433774561 on the coarse mesh.
  None of these values is a certified lower bound.
- The rational conservative uniform-mesh error bound fails: gamma squared
  is 72.857106978 at one planar refinement and 12 layers. A favorable
  floating spectrum cannot override this failed input.
- Adaptive local refinement produces 909,276 tetrahedron leaves with
  gamma squared <=0.099999915497 and sigma squared <=0.278166092189.
  With eta=1/5, theta=9/10 and rho=5, the exact scalar error budget has
  c_e>0.24546. This is a **local bound/scalar pass**, conditional on the
  eventual compatible finite space and matrix inequality.
- The compressed refinement plan includes exact rational r layers and
  leaf paths. Its separate verifier checks bindings, prefix freedom and
  exhaustive recursive coverage. It does **not** verify hanging-face
  moment aggregation, coefficient replay, or matrix positivity.

The initial top-mean-only diagnostic failed because that is an unnecessarily
stronger criterion: nearly constant bulk trials can have small top mean.
The threshold form's constant-direction projection uses the core mean plus
one quarter of the top mean. The failure is not evidence of an exceptional
quotient eigenvalue.

## Replay

From repository root with the existing requirements installed:

```sh
python artifacts/track_b/d67/spectral/verify_reference_mesh.py
python -O artifacts/track_b/d67/spectral/verify_reference_mesh.py
python artifacts/track_b/d67/spectral/negative_checks.py
python artifacts/track_b/d67/spectral/verify_adaptive_plan.py
python -O artifacts/track_b/d67/spectral/verify_adaptive_plan.py
python artifacts/track_b/d67/spectral/certify_status.py
python tests/d19_regressions.py replay
python tests/d19_regressions.py protected
```

The status command exits **2**, deliberately reporting the four missing
spectral conditions. It is a fail-closed progress gate, not a final theorem
verifier. No JSON status flag can turn the stored diagnostics into a pass.

Independent computations (the adaptive run takes about eight minutes on
this execution host) can be written to fresh paths to preserve the frozen
reports:

```sh
OPENBLAS_NUM_THREADS=1 python artifacts/track_b/d67/spectral/threshold_probe.py --refinements 1 --layers 12 --rho 3
OPENBLAS_NUM_THREADS=1 python artifacts/track_b/d67/spectral/adaptive_error_bounds.py --target 1/10 --layers 12 --max-leaves 1000000 --output /tmp/d67-adaptive-replay.json
```

The adaptive producer reconstructs geometry rather than taking the stored
leaf coordinates or scalar bounds on trust. Exact gamma/sigma bounds and
leaf paths are reproducible; timing and floating renderings are not proof
inputs. The companion replay plan is saved alongside the requested output.

## Independent adaptive face-moment replay

The later checkpoint now verifies the global prolongation P on the frozen
909,276-leaf plan. `moment_map_verification.json` binds the input, topology,
rows and producer/verifier sources by SHA-256. Replay regenerates the input
from exact geometry, verifies the initial-root ordering and the refinement
coverage, then checks the complete topology and rows. All 2,042,568 master
variables have identity-row witnesses, which proves P has full column rank.
There are 4,084,872 prolongation entries; the largest row has 994 entries.

```sh
python -O artifacts/track_b/d67/spectral/replay_moment_map.py /tmp/d67-moments --rebuild
python artifacts/track_b/d67/spectral/moment_negative_checks.py /tmp/d67-moments
python artifacts/track_b/d67/spectral/centroid_bound.py
python -O artifacts/track_b/d67/spectral/tuned_error_bounds.py --checkpoint /tmp/d67-scalar-checkpoint.json --output /tmp/d67-tuned.json
python -O artifacts/track_b/d67/spectral/index_budget.py
```

The C++ producer propagates barycentric coordinates; its separate verifier
recovers them geometrically from exact integer coordinate minors. It also
reconstructs each leaf from its root/path and checks exact quarter partitions,
physical incidence, canonical partitions and every dyadic row weight.
`replay_moment_map.py` binds the root pairing frames to independently checked
Ford geometry. Binaries are reproducible intermediates, generated outside
the repository. The binary format requires a little-endian GCC-compatible
host with signed 128-bit integer support. Undefined-behavior sanitizer replay
also passed before the final additional canonical-child validation.

`docs/D67_TRACK_B_CENTROID_CR_BOUND.md` proves the sharper self-contained
constant 1661/15000. This gives c_e > 0.21858 with eta=1/10 on the original
coefficient envelope, changing the finite mass factor to 11/10.

The full-mesh vertical polynomial diagnostics identify a remaining coefficient
issue. Retaining half the vertical energy gives a negative restricted form;
retaining three quarters gives a positive restricted form and a positive
floating scalar margin. These are floating diagnostics on a conservative
box/hull envelope. They prove neither positivity of the full matrix nor an
obstruction for the frozen exact envelope. The three-quarter envelope now has
a new exact scalar ledger in `tuned_error_bounds.json`: gamma² <= 0.070917340541
and sigma² <= 0.159181690860. The old-style exact scalar margin is >0.13147.

## Next required evidence

The direct negative-index transfer in `docs/D67_TRACK_B_DISCRETE_INDEX.md`
removes the weighted-mean enclosure requirement and the sigma penalty. Its
exact scalar margin is c0 > 0.21990 (`index_budget.py`). On the verified P
and frozen leaf mesh, assemble and enclose

    C_h = Q_h+(1/2)b_h-(3/4)t_h^2-(11/10)M_h.

Prove C_h has at most one negative direction, for example by certifying
C_h+alpha*z*z^T is positive semidefinite for any explicit rational z and
alpha>0. The candidate z=M_h*1+t_h/4 needs only rational mass and top-area
data, not geometric mean integration. Its polynomial subspace diagnostic
is positive; it is not a full matrix certificate. P already has a verified
trivial representation kernel. Reconstruct all coefficient bounds and matrix
entries in an independent verifier. Only that completed evidence
can connect the scalar pass to the threshold theorem. The earlier
`certify_status.py` and `status.json` remain frozen historical progress gates;
their face-map missing condition is superseded by the later independent replay.

The complete floating sparse matrix is assembled by `adaptive_matrix_probe.py`.
`matrix_assembly.json` records 30,786,620 energy entries and 30,786,664 mass
entries. Top area is 1 and constant energy is 1/2 to floating precision;
the maximum symmetry defect is below 9e-16. All are diagnostics, never
interval enclosures or positivity evidence.

```sh
OPENBLAS_NUM_THREADS=1 python artifacts/track_b/d67/spectral/adaptive_matrix_probe.py /tmp/d67-moments/full
# The optional eigenprobe additionally requires pyamg 5.2.1.
OPENBLAS_NUM_THREADS=1 python -u artifacts/track_b/d67/spectral/adaptive_matrix_probe.py /tmp/d67-moments/full --eigenprobe
```

The two dense rank updates stay factored. The sparse matrices and eigenvectors
are reproducible intermediate files and are not committed. A floating
eigenvalue estimate cannot close the exclusion gate.

The earlier geometry, Gaussian artifacts, d=19 certificates, historical
reports and known failing tests are unchanged. No d=43/d=163 spectral
statement follows from this checkpoint.


## Latest balanced-enclosure and subcell-mass checkpoint

D67_TRACK_B_BALANCED_ENVELOPE.md proves an exact 9/8 inverse-metric comparison
with the original ledger. Combined with the centroid improvement it gives
gamma² <= 0.070073378286 and a direct-index scalar margin > 0.22919.
The comparison uses the original exact coefficient intervals and h0.
The diagnostic's box/hull relaxation is not a certified replacement for
those inputs.

The exact energy-kernel check merges 9,697 identity-connected classes into
one in two passes. Thus the local CR energy has only the constant kernel.
This gives no threshold lower bound.

The full balanced floating matrix still has a negative test direction.
Two levels of integration-only red subdivision sharpen its finite mass
upper bound without changing P or gamma. The total upper mass falls from
about 5.75509 to 4.14379; the previously negative trial becomes positive.
The whole-matrix diagnostic nevertheless reaches a negative estimate
-0.02982259 by iteration 16. It is not a rigorous bound or a completed
spectral replay. The command runner stopped responding before the next
refinement could be completed. latest_progress_status.json records the
remaining gates and keeps the spectral flag false.

Reproduce in a fresh directory after the exact moment replay:

~~~sh
python -O artifacts/track_b/d67/spectral/balanced_envelope.py
g++ -std=c++17 -O2 artifacts/track_b/d67/spectral/moment_kernel_check.cpp -o /tmp/d67-kernel
/tmp/d67-kernel /tmp/d67-moments/full.rows.bin
OPENBLAS_NUM_THREADS=1 python -u artifacts/track_b/d67/spectral/adaptive_matrix_probe.py /tmp/d67-moments/full --balanced --mass-depth 2
OPENBLAS_NUM_THREADS=1 python -u artifacts/track_b/d67/spectral/adaptive_matrix_probe.py /tmp/d67-moments/full --eigenprobe --iterations 20
~~~

The probe additionally needs pyamg 5.2.1. Floating matrix data and eigenvectors
are reproducible intermediates. Neither positive polynomial trials nor a
floating eigenvalue estimate can issue a spectral certificate.
