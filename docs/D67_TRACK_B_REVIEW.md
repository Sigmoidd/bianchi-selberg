# d=67 Track B review and replay

## Claim boundary

Geometry is certified for full PSL₂(O₋67), level one. The target of
excluding every discrete eigenvalue in (0,1) remains **OPEN**.
The all-H¹(K) compact-core relaxation is rigorously obstructed by v=b
at s=1/2. This is not a quotient eigenfunction or an impossibility result
for the glued-core method. See `D67_TRACK_B_EXCLUSION_STATUS.md`.

| Package | Evidence | Status |
|---|---|---|
| Fundamental region, pairings, pointwise stabilizers, cusp | Geometry proof and exact independent replay | Closed |
| Cusp reduction to B_s on matching traces | Exclusion status proof, standard analytic dependencies | Conditional exclusion criterion established |
| Guaranteed glued-core discretization theorem | None for this field | Open |
| Validated mesh, quadrature, interpolation and slivers | None for this field | Open |
| Positivity on all parameter windows | None for this field | Open |
| Independent final spectral theorem replay | No final theorem artifact | Open |

## Replay

From repository root, with `requirements-trace.txt` installed:

```sh
python artifacts/track_b/d67/verify_cover.py
python artifacts/track_b/d67/geometry_certify.py
python artifacts/track_b/d67/verify_geometry.py
python -O artifacts/track_b/d67/verify_geometry.py
python artifacts/track_b/d67/geometry_negative_checks.py
python artifacts/track_b/d67/negative_checks.py
python artifacts/track_b/d67/relaxation_obstruction.py
python -O artifacts/track_b/d67/relaxation_obstruction.py
python tests/d19_regressions.py replay
python tests/d19_regressions.py protected
```

Geometry replay reconstructs the candidate set and every polygon before
checking the serialized ledger. The producer is unnecessary for replay.
The verifier independently enumerates vertex groups and their relations,
checks edge incidence and endpoint-group intersections, and validates all
face images by exact affine comparisons against every candidate.
Normal and optimized outputs agree. Eight new geometry mutations and
four existing cover mutations are rejected without changing stored data.

The frozen design and cover retain their original hashes. The probe and
cover reports retain their original narrower scopes. No existing
certificate, known failing test, Gaussian ledger, or historical report
is modified. Only field-independent clipping/group-check algorithms and
cusp analysis are reused from d43. All norm bounds, rows, polygons,
matrices, point groups, area/frequency constants and obstruction are new
for d67. Nothing here certifies d163.

## Baseline results

Full command outputs are in `artifacts/track_b/d67/baseline/`.

| Command | Result |
|---|---|
| `python -m unittest discover -v` | 68 tests, one known failure |
| `python -m unittest discover -s tests -v` | 54 tests, two known failures |
| `python -m unittest test_track_b_two_cusp -v` | 8 tests, one known failure |
| `python tests/mutation/regression_checks.py triage` | exit 0; stored Gaussian verifier remains false |
| `python -O -m groups.systoles` | Pass |
| `python -O -m groups.d19_inventory` | Pass |
| `python tests/d19_regressions.py replay` | Pass; frozen baseline gate satisfied |
| `python tests/d19_regressions.py protected` | Pass; 152 historical files unchanged |
| `python tests/mutation/regression_checks.py links --external` | Initial network denial; retry with network access passes |
| `git diff --check` | Pass |

Preserved blockers, verbatim:

```text
FAIL: test_closing_artifacts_verify_independently
AssertionError: False is not true
collocation_ledger_hash FAIL actual: 8ddb66397840d9ebddca7ac2d4e8810076ec031c0931ea4c21cc0940c66c6486 stored: 383557121919fcddfb8a3a6b8592238f612eca6b1f6bba3bc4583d15611281a4
FAIL: test_status_string_cannot_certify_another_field
AssertionError: ValueError not raised
FAIL: test_new_fields_block_full_assembly
AssertionError: ValueError not raised
UNVERIFIED <urlopen error [Errno 1] Operation not permitted>
```

The last line is an execution-environment issue resolved by the logged
retry. The other failures existed on the base and remain preserved.
No stored hash was repaired to turn a historical failure into a pass.
