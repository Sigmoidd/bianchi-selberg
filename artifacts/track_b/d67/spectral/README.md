# d=67 threshold route checkpoint

**Spectral exclusion remains OPEN.** The first remaining gate is an exact
adaptive face-fragment moment map. No interval threshold matrices, verified
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

## Next required evidence

Reconstruct the nested face-fragment partition, every pairing orbit
(including self-pairings), and its exact area-weighted prolongation P.
Then assemble and enclose, on that same P and frozen leaf mesh,

    Q_h+(1/2)b_h-(3/4)t_h^2+(1/2)L_h^2-(6/5)M_h.

Certify this matrix is positive semidefinite after an exact removal of any
representation kernel. Reconstruct all coefficient/mean-vector bounds and
matrix entries in an independent verifier. Only that completed evidence
can connect the scalar pass to the threshold theorem.

The earlier geometry, Gaussian artifacts, d=19 certificates, historical
reports and known failing tests are unchanged. No d=43/d=163 spectral
statement follows from this checkpoint.
