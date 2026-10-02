# d43 Track B design freeze and baseline

The d43 fundamental-region geometry gate is now proved by exact replay.
See [D43_TRACK_B_GEOMETRY.md](../../../docs/D43_TRACK_B_GEOMETRY.md)
for the coverage, finite cutoff, face, pairing and stabilizer arguments.
This remains **no spectral exclusion certificate**: d43, d67 and d163
are OPEN. The next field-specific theorem input is compact-core lower-bound
theory and a validated mesh. The PID correspondence and spectral citations
also need source pinning before the final theorem.

Base commit: `86048a49dafb10345993152ef7a69682822ee393`.
Branch: `track-b/d43-design`. Design:
[D43_TRACK_B_DESIGN.md](../../../docs/D43_TRACK_B_DESIGN.md).
`arithmetic.json` binds its exact bytes by SHA-256.

## Recorded baseline

`baseline/status.json` records commands, exit codes and elapsed time;
the corresponding `.log` files preserve exact stdout and stderr.
Dependencies were already satisfied by `requirements-trace.txt`.
The initial dedicated `two-cusp.log` is empty despite recorded exit 1;
the command was rerun and its full output saved in `two-cusp-recheck.log`
(8 tests, 1 failure). Root discovery independently records the same failure.

| Check | Result |
|---|---|
| Root unittest discovery | 68 tests; 1 known failure |
| tests/ unittest discovery | 54 tests; 2 known failures |
| test_track_b_two_cusp | known artifact failure |
| triage | process exit 0; verifier reports `verified: false` |
| Optimized systoles and d19 inventory | PASS |
| Frozen d19 replay | PASS |
| d19 protected artifacts | PASS; 152 protected files, changed [] |
| Links, including external | PASS |
| git diff --check | PASS |

The failing test names and messages are preserved verbatim in the logs:

```text
test_closing_artifacts_verify_independently
AssertionError: False is not true

test_status_string_cannot_certify_another_field
AssertionError: ValueError not raised

test_new_fields_block_full_assembly
AssertionError: ValueError not raised
```

The triage log records `deterministic_hashes_match: false`,
`physical_residual_certified: false`, and `rung4_certified: false`.
Its zero process exit code does not make that artifact pass. The stored
collocation-ledger hash remains unchanged. Mutations remain gated.

## Independent design-input replay

From the repository root:

```sh
python artifacts/track_b/d43/replay_arithmetic.py
python -O artifacts/track_b/d43/replay_arithmetic.py
```

Both outputs were compared byte-for-byte and agree with `arithmetic.json`.
The script checks exact group identity, integral relation, reduced forms,
translation matrices, dual pairings and the rational sufficient cusp
mode bound. Universal norm/shortest-frequency arguments are in the design;
finite samples in the script are consistency checks, not substitutes for
those arguments. The script explicitly reports `exclusion_certified=false`.

The design additionally records an obstruction to directly reusing the
Gaussian floor-inclusion lemma: the stated universal triangle variance
bound d²/4 fails on equilateral triangles. No historical proof file has
been edited. Correct inclusion and lower-bound arguments for the actual
d43 geometry are prerequisites to any matrix certificate.

## Next authorized work package

Adapt and prove the compact-core lower-bound theorem on the certified
19-patch floor with its exact face identifications. Build its field-specific
mesh and error bounds, then verify positivity across the entire exceptional
interval. No conclusion for d67 or d163 transfers from this geometry.
