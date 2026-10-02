# d67 Track B evidence

**Geometry gate: proved. Spectral exclusion: OPEN.**

[Review guide and baseline record](../../../docs/D67_TRACK_B_REVIEW.md).

The original [frozen design](../../../docs/D67_TRACK_B_DESIGN.md) is
unchanged. The [geometry proof](../../../docs/D67_TRACK_B_GEOMETRY.md)
closes its finite-domain dependency. The
[exclusion status](../../../docs/D67_TRACK_B_EXCLUSION_STATUS.md)
proves a compact-core obstruction to relaxing all face identifications
and identifies the first missing theorem input.

From repository root:

```sh
python artifacts/track_b/d67/verify_cover.py
python artifacts/track_b/d67/geometry_certify.py
python artifacts/track_b/d67/verify_geometry.py
python -O artifacts/track_b/d67/verify_geometry.py
python artifacts/track_b/d67/geometry_negative_checks.py
python artifacts/track_b/d67/negative_checks.py
python artifacts/track_b/d67/relaxation_obstruction.py
python tests/d19_regressions.py replay
python tests/d19_regressions.py protected
```

The replay reconstructs all 1,575 candidates, 37 positive-area floor
patches, 102 edges, 66 vertices, pairing matrices and pointwise stabilizer
groups. The minimum floor height squared is 2/67. Normal and optimized
replays agree. Geometry mutations have their own semantic checks.
The older `face_probe_result.json` remains a diagnostic historical
record; `geometry_result.json` and `verification.json` supersede its
open geometry flags. The cover result certifies only the cover scope.

`relaxation_obstruction.json` certifies B_(1/2)<0 for v=b on all H¹(K).
This trial does not descend to the quotient; it excludes the relaxed
criterion, not the spectral target. No Gaussian or d43 geometry or
numerical constants transfer. Only the exact polygon algorithm, group
relation checks, and analytic cusp reduction are reused.

Historical certificates and known failing tests are preserved. The
frozen d19 replay and protected-file check pass; see their logs.
