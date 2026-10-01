# Regression triage and mutation audit

## Scope and environment

Tested d=7 production commit: `d3610366047b3e178942e2c7b9874ca51630b8b5`.
PR [#5](https://github.com/Sigmoidd/bianchi-selberg/pull/5) is a draft targeting
`codex/group-backends-and-proof-docs`, base commit
`70c01bb5f0542141fe0fa61451921273be4fbbe5` (confirmed from PR metadata).
The earlier stop-only report was committed separately as `8e8d40c`.
This revised audit changes only this report and files in `tests/mutation/`.
Production math, frozen certificates, historical reports, and existing
unit tests were not edited.

Runtime: Python **3.12.14**, python-flint **0.9.0**, FLINT **3.6.0**.
Arb is included in this FLINT runtime; a separate Arb version is not exposed
by the Python module. Frozen certificates also record Python 3.12.14 and
python-flint 0.9.0. The optimized arithmetic commands used `python -O`;
other commands used ordinary Python.

## A. Triage of the root-level failure

**Fails on base too; not newly introduced by the d=7 branch.** Both worktrees
ran `python -m unittest test_track_b_two_cusp -v`: 8 tests, 1 failure, 7 passes.
Full captured outputs follow below and are also preserved in
`tests/mutation/head-two-cusp.log` and `tests/mutation/base-two-cusp.log`.

`verify_from_paths` reads:

| Input | Exists | Git tracked | Purpose |
|---|---|---|---|
| `track_b_two_cusp_result.json` | yes | yes | Assembly definition, coefficient enclosures, claimed residuals and hashes |
| `track_b_hejhal_rows.jsonl` | yes | yes | Collocation rows; raw bytes supply the ledger SHA-256 |

These paths resolve relative to the repository working directory, exactly
as the existing test intends. The artifact describes Gaussian integers
`Z[i]`, level `(2+i)`, index 6, with infinity and zero cusps. It does **not**
verify a d=2 or d=7 level-one certificate, and thus does not block Step C.

The failed check is **`deterministic_hashes_match`**, specifically the raw
ledger SHA-256 comparison:

```text
actual: 8ddb66397840d9ebddca7ac2d4e8810076ec031c0931ea4c21cc0940c66c6486
stored: 383557121919fcddfb8a3a6b8592238f612eca6b1f6bba3bc4583d15611281a4
```

Assembly, cusp-data, and coefficient-vector hashes match. All other returned
verification gates pass: exact row count/order, all four blocks, glue,
transformed coordinates, coefficient count, selected-system zero containment,
scaling roundtrip, and all three endpoint comparisons. `verified=False` and
`physical_residual_certified=False` follow from the failed ledger hash;
`rung4_certified=False` is the expected artifact status. The full verifier
output and hash diagnostics are in `tests/mutation/triage.log`. No cause
for the byte-level ledger mismatch beyond this comparison was inferred.

Root discovery (`python -m unittest discover -v`) intentionally matches the
root-level `test_track_b_two_cusp.py` through the default `test*.py` pattern;
it finds **68 tests**. This was the command used in the earlier run, which
failed once. Explicit `tests/` discovery finds **44 different tests**, includes
the d=2/d=7/quotient checks, and does not include the root two-cusp module.
The PR's 44-test claim has that narrower scope. The revised run executed
those 44 tests successfully rather than treating root discovery as complete
coverage of the separate directory.

## B. Remaining regressions

Commands below are from the repository root unless otherwise stated.
`OUT=/workspace/scratch/31144417503a/regression-logs` was the output directory
used in this run; substitute another temporary directory when reproducing.

| Item | Result | Exact operative command(s) |
|---|---|---|
| 2. Optimized arithmetic replay | PASS, all exit 0 | `python -O -m groups.d2_inventory`; `python -O -m groups.d7_inventory`; `python -O -m groups.systoles` |
| 3. Exact frozen replay | PASS | `python tests/mutation/regression_checks.py replay` |
| 3. v1/historical byte identity against PR base | PASS | `python tests/mutation/regression_checks.py protected` |
| 4. Quotient operators | PASS, 4 tests | `python -m unittest discover -s tests -p test_quotients.py -v` |
| 4. Quotient CLI | PASS, 360 vertices, degree 6 | `python examples/finite_quotient.py 7 --level-generator 3 0 --output "$OUT/d7-mod3.json"` |
| 5. d=2 script | PASS | `python examples/d2_certificate.py --output "$OUT/d2-script.json"` |
| 5. d=2 module | PASS | `python -m examples.d2_certificate --output "$OUT/d2-module.json"` |
| 5. d=7 script | PASS | `python examples/d7_certificate.py --output "$OUT/d7-script.json"` |
| 5. d=7 module | PASS | `python -m examples.d7_certificate --output "$OUT/d7-module.json"` |
| 6. Documentation file links and external reachability | PASS in stated scope | `python tests/mutation/regression_checks.py links --external` |
| 6. Whitespace | PASS, no output | `git diff --check` |
| Additional explicit test-directory discovery | PASS, 44 tests | `python -m unittest discover -s tests -v` |

All four CLI JSON objects equal their corresponding frozen certificate
objects in full (checked by loading both JSON files and comparing objects):
`python tests/mutation/regression_checks.py cli --output-dir /workspace/scratch/31144417503a/regression-logs`.
The protected comparison confirms unchanged d=2 v1, d=2 v2, legacy-pre-M0,
M0 regression, mechanical screen, and `old_RIGOR_GAPS.md`. Additionally, **all
134 pre-existing tracked JSON/JSONL artifacts** equal the PR base byte for
byte. `certificates/README.md` was changed by the d=7 PR to document the new
artifacts; it is documentation, not a frozen certificate/historical report.
It was not modified in this audit.

Documentation checking initially covered 40 local inline-link file targets
in `docs/*.md`, README, structure guide, proof ledger and certificate README,
with no broken paths. All 8 unique external URLs returned HTTP 200, including
the DOI redirect (which carries a cookie-error query). This checks file
existence and HTTP reachability, not heading fragments or remote content
correctness. The final check covered 42 local targets and all 9 external URLs (adding
the PR link); all passed. Heading fragments remain outside this check.

Exact bound and endpoint strings are unchanged:

| Group | Bound ball | Lower endpoint ball | Upper endpoint ball |
|---|---|---|---|
| d=2 | `[0.43 +/- 5.45e-3]` | `[0.4245518407392955617061859991770204379294 +/- 4.93e-41]` | `[0.4295747809360909507623917486887391879294 +/- 4.93e-41]` |
| d=7 | `[0.3640 +/- 5.32e-5]` | `[0.3639468761489975954506425851303745125833 +/- 1.95e-41]` | `[0.3640005536935969907832677209097553109232 +/- 3.69e-41]` |

Outward eight-decimal enclosures are d=2 **[0.42455184, 0.42957479]** and
d=7 **[0.36394687, 0.36400056]**. Exact comparisons use the frozen strings,
not these rounded display enclosures.

## C. Gate, predictions, and mutation method

**Gate passed:** optimized arithmetic and both exact frozen-certificate
replays passed before mutation execution.

[Pre-run predictions](../tests/mutation/PREDICTIONS.md) were written before
the mutation campaign. For both groups:

| Mutation | Predicted direction | Reason |
|---|---|---|
| Drop each elliptic class | down | Removes a positive NCE summand. |
| Centralizer denominator ×2 | down | Halves the positive NCE block. |
| Centralizer denominator ×1/2 | up | Doubles the positive NCE block. |
| Remove each/all primitive translation witnesses | zero; arithmetic rejection | Stored norms supply coefficients independently of witness matrices. |
| Volume +1% / -1% | up / down | Positive identity contribution is linear in volume. |
| Shrink exhaustive loxodromic trace box | zero if accepted | Shortest witness survives; no explicit geodesic sum is assembled at this support. |
| Omit cusp/scattering block | up | Baseline aggregate block is negative. |
| Flip cusp/scattering sign | up | Replaces that negative block with its opposite. |

Each mutation starts from a fresh unmodified group and uses only in-memory
monkeypatches or dataclass replacements. `evaluate` and `certificate_payload`
are called with the frozen parameters. Centralizer scaling is run both on
stored class orders and on the coefficient denominator with records intact.
The latter tests trusted numerical code, rather than record binding.

The cusp/scattering block is exactly `CE+Ch0+PARg0+PSI+PHIINT`.
`prime` is a diagnostic already included in PHIINT and is not counted twice.
Truncation replaces the exhaustive systole trace domain `[-6,6]^2` with
`[-1,1]^2`; this is the production loxodromic enumeration that justifies
omitting the geometric sum. There is no separate nonzero geodesic-sum
implementation in this pipeline to truncate. The checked trace count falls
from 18 to 6 for d=2 and from 20 to 6 for d=7, retaining the shortest witnesses.

**D** below means arithmetic rejected the mutation before producing B.
Its table interval comes from a separate diagnostic `evaluate` call with
only `GroupData.require_inventory` bypassed in memory; no diagnostic
certificate was exported. **N** means the native evaluation produced B;
the exporter was then attempted. Calling these diagnostic values a valid
spectral certificate would be incorrect.

Catching checks are native production checks. “Positivity check” in this
table specifically means the exporter's `B < 1` criterion; no independent
spectral positivity theorem is being tested. “Arithmetic replay” means
inventory witness/record verification. Frozen-output equality is a separate
regression check, not counted as a native mutation rejection. Trusted-code
monkeypatches remaining uncaught do not by themselves establish a flaw in
the unmodified mathematical derivation.

Intervals in the table are rounded outward to eight decimals; exact balls,
endpoint balls, exception messages and metadata are retained in
[results.json](../tests/mutation/results.json). Delta uses conservative
interval subtraction of mutated and baseline balls. For byte-identical
midpoints and radii, reported deterministic numerical movement is exactly
zero; this does not narrow the uncertainty of B itself. Radius flags compare
midpoint movement with the **actual Arb radius**, not the radius of the
coarsely formatted `bound_ball` string.

## Mutation results

Flags: **R** = numerical movement below baseline Arb radius; **U** = no native
check caught the mutation and mutated B is provably below 1; **O** = opposite
to prediction. Every row compares in the predicted direction (including
predicted zero movement).

| Group | Mutation | Baseline B | Mutated B | Delta | Check caught | Source | Flags |
|---|---|---|---|---|---|---|---|
| d=2 | drop:0:d2 order 2: multiplier A2 | [0.42455184, 0.42957479] | [0.20124532, 0.20626827] | [-0.22832946, -0.21828357] | arithmetic replay | D | — |
| d=2 | drop:1:d2 order 2: multiplier S2 | [0.42455184, 0.42957479] | [0.31289858, 0.31792153] | [-0.11667620, -0.10663031] | arithmetic replay | D | — |
| d=2 | drop:2:d2 order 3: alpha | [0.42455184, 0.42957479] | [0.16641187, 0.17143482] | [-0.26316291, -0.25311702] | arithmetic replay | D | — |
| d=2 | drop:3:d2 order 3: alpha inverse | [0.42455184, 0.42957479] | [0.16641187, 0.17143482] | [-0.26316291, -0.25311702] | arithmetic replay | D | — |
| d=2 | centralizer denominator ×2 | [0.42455184, 0.42957479] | [-0.00106802, 0.00395493] | [-0.43064280, -0.42059691] | none | N | U |
| d=2 | centralizer denominator ×0.5 | [0.42455184, 0.42957479] | [1.27579154, 1.28081449] | [0.84621676, 0.85626265] | positivity check | N | — |
| d=2 | centralizer records ×2 | [0.42455184, 0.42957479] | [-0.00106802, 0.00395493] | [-0.43064280, -0.42059691] | arithmetic replay | D | — |
| d=2 | centralizer records ×0.5 | [0.42455184, 0.42957479] | [1.27579154, 1.28081449] | [0.84621676, 0.85626265] | arithmetic replay | D | — |
| d=2 | translation:0:d2 order 2: multiplier A2 | [0.42455184, 0.42957479] | [0.42455184, 0.42957479] | [0, 0] | arithmetic replay | D | R |
| d=2 | translation:1:d2 order 2: multiplier S2 | [0.42455184, 0.42957479] | [0.42455184, 0.42957479] | [0, 0] | arithmetic replay | D | R |
| d=2 | translation:2:d2 order 3: alpha | [0.42455184, 0.42957479] | [0.42455184, 0.42957479] | [0, 0] | arithmetic replay | D | R |
| d=2 | translation:3:d2 order 3: alpha inverse | [0.42455184, 0.42957479] | [0.42455184, 0.42957479] | [0, 0] | arithmetic replay | D | R |
| d=2 | remove all translation witnesses | [0.42455184, 0.42957479] | [0.42455184, 0.42957479] | [0, 0] | arithmetic replay | D | R |
| d=2 | volume +1% | [0.42455184, 0.42957479] | [0.43577701, 0.44079996] | [0.00620222, 0.01624811] | none | N | U |
| d=2 | volume -1% | [0.42455184, 0.42957479] | [0.41332667, 0.41834962] | [-0.01624811, -0.00620222] | none | N | U |
| d=2 | trace box [-6,6] → [-1,1] | [0.42455184, 0.42957479] | [0.42455184, 0.42957479] | [0, 0] | none | N | R, U |
| d=2 | omit cusp/scattering | [0.42455184, 0.42957479] | [0.89924853, 0.89924854] | [0.46967375, 0.47469670] | none | N | U |
| d=2 | flip cusp/scattering | [0.42455184, 0.42957479] | [1.36892229, 1.37394524] | [0.93934751, 0.94939340] | positivity check | N | — |
| d=7 | drop:0:d7 order 2: A2 | [0.36394687, 0.36400056] | [-0.00055704, -0.00050335] | [-0.36455759, -0.36445023] | arithmetic replay | D | — |
| d=7 | drop:1:d7 order 3: inverse classes merged | [0.36394687, 0.36400056] | [-0.00276375, -0.00271006] | [-0.36676430, -0.36665694] | arithmetic replay | D | — |
| d=7 | centralizer denominator ×2 | [0.36394687, 0.36400056] | [-0.00166039, -0.00160671] | [-0.36566095, -0.36555358] | none | N | U |
| d=7 | centralizer denominator ×0.5 | [0.36394687, 0.36400056] | [1.09516140, 1.09521509] | [0.73116085, 0.73126821] | positivity check | N | — |
| d=7 | centralizer records ×2 | [0.36394687, 0.36400056] | [-0.00166039, -0.00160671] | [-0.36566095, -0.36555358] | arithmetic replay | D | — |
| d=7 | centralizer records ×0.5 | [0.36394687, 0.36400056] | [1.09516140, 1.09521509] | [0.73116085, 0.73126821] | arithmetic replay | D | — |
| d=7 | translation:0:d7 order 2: A2 | [0.36394687, 0.36400056] | [0.36394687, 0.36400056] | [0, 0] | arithmetic replay | D | R |
| d=7 | translation:1:d7 order 3: inverse classes merged | [0.36394687, 0.36400056] | [0.36394687, 0.36400056] | [0, 0] | arithmetic replay | D | R |
| d=7 | remove all translation witnesses | [0.36394687, 0.36400056] | [0.36394687, 0.36400056] | [0, 0] | arithmetic replay | D | R |
| d=7 | volume +1% | [0.36394687, 0.36400056] | [0.37510399, 0.37515768] | [0.01110344, 0.01121080] | none | N | U |
| d=7 | volume -1% | [0.36394687, 0.36400056] | [0.35278975, 0.35284344] | [-0.01121080, -0.01110344] | none | N | U |
| d=7 | trace box [-6,6] → [-1,1] | [0.36394687, 0.36400056] | [0.36394687, 0.36400056] | [0, 0] | none | N | R, U |
| d=7 | omit cusp/scattering | [0.36394687, 0.36400056] | [0.77810715, 0.77810716] | [0.41410660, 0.41416028] | none | N | U |
| d=7 | flip cusp/scattering | [0.36394687, 0.36400056] | [1.19221375, 1.19226744] | [0.82821320, 0.82832056] | positivity check | N | — |

## Explicit flags and limits

- **Below radius: 10 cases.** Removing each primitive translation or the
  entire witness set (5 d=2 cases, 3 d=7 cases), plus both truncated trace
  domains, causes zero deterministic movement in B. Translation matrices
  are exercised by arithmetic guards, but not by the numerical coefficient
  formula; stored primitive norms remain intact. The omitted loxodromic
  numerical term is not exercised at this support. The truncated domain
  itself is accepted by the monkeypatched verifier.
- **Opposite prediction: none observed** in 32 completed cases.
- **Uncaught with B<1: 10 cases.** For each group: coefficient denominator
  ×2 with exact records intact; volume +1%; volume -1%; truncated trace box;
  omitted cusp/scattering block. These are flagged without fixes.
- Stored centralizer normalization changes, missing classes, and missing
  translation witnesses trigger arithmetic guards: **18 cases**. Scaling
  the trusted coefficient formula leaves the stored records valid, so its
  denominator-×2 case is uncaught; denominator-×1/2 fails B<1 for both groups.
- **Four B<1 rejections:** coefficient denominator ×1/2 and cusp/scattering
  sign flip, for each group. **No completeness rejection** occurred for
  the truncated trace domain.

The d=7 coefficient-denominator ×2 mutation is accepted with an entirely
negative B interval; the exporter checks the upper threshold B<1 and does
not reject that negative interval in this run. This is included in the U flags.

These results measure the chosen parameters and specified perturbations.
They do not prove universal mutation coverage, mathematical independence
of the checks, or anything about other fields or levels. Nothing flagged
here was fixed.

## Verbatim native mutation rejection messages

```text
ArithmeticError: complete d7 class count
ArithmeticError: d7 record differs: finite_centralizer_order
ArithmeticError: inventory record count
ArithmeticError: inventory record differs: finite_centralizer_order
ArithmeticError: translation witness is missing or outside the specified group
ValueError: this test function does not prove B < 1
```

Per-case attribution is retained in results.json.

## Exact reproduction commands

Run from the d=7 checkout; use a separate base worktree for comparison:

```sh
OUT=/workspace/scratch/31144417503a/regression-logs
mkdir -p "$OUT"
python -m unittest test_track_b_two_cusp -v > "$OUT/head-two-cusp.log" 2>&1
python tests/mutation/regression_checks.py triage > "$OUT/triage.log"
git fetch origin codex/group-backends-and-proof-docs
git worktree add --detach /workspace/scratch/31144417503a/bianchi-base 70c01bb5f0542141fe0fa61451921273be4fbbe5
```

From `/workspace/scratch/31144417503a/bianchi-base`:

```sh
python -m unittest test_track_b_two_cusp -v > /workspace/scratch/31144417503a/regression-logs/base-two-cusp.log 2>&1
```

From the d=7 checkout again (commands are independent; continue after a
nonzero exit to collect all regressions):

```sh
python -m unittest discover -s tests -v
python -O -m groups.d2_inventory
python -O -m groups.d7_inventory
python -O -m groups.systoles
python tests/mutation/regression_checks.py replay
python tests/mutation/regression_checks.py protected
python -m unittest discover -s tests -p test_quotients.py -v
python examples/finite_quotient.py 7 --level-generator 3 0 --output "$OUT/d7-mod3.json"
python examples/d2_certificate.py --output "$OUT/d2-script.json"
python -m examples.d2_certificate --output "$OUT/d2-module.json"
python examples/d7_certificate.py --output "$OUT/d7-script.json"
python -m examples.d7_certificate --output "$OUT/d7-module.json"
python tests/mutation/regression_checks.py cli --output-dir "$OUT"
python tests/mutation/regression_checks.py links --external
git diff --check
```

Only after optimized arithmetic and exact frozen replay pass:

```sh
python tests/mutation/run.py "$OUT/mutations.json" > "$OUT/mutations.log" 2>&1
```

`run.py` independently rechecks frozen B, both endpoint strings and parameters
before any mutations for each group; it raises on baseline mismatch. Read
PREDICTIONS.md before running. Do not overwrite `tests/mutation/results.json`
unless deliberately recording a new audit. The first harness pass incorrectly
used Arb `==` for nonzero-radius ball identity; this was corrected to exact
midpoint/radius comparison and rerun. The final retained results contain
32 cases, including both centralizer-record and formula perturbations.

## Full targeted unittest output: d=7 branch

```text
test_all_four_blocks_and_slow_direct_rows_agree (test_track_b_two_cusp.TrackBTwoCuspTests.test_all_four_blocks_and_slow_direct_rows_agree) ... ok
test_closing_artifacts_verify_independently (test_track_b_two_cusp.TrackBTwoCuspTests.test_closing_artifacts_verify_independently) ... FAIL
test_each_certification_gate_fails_closed (test_track_b_two_cusp.TrackBTwoCuspTests.test_each_certification_gate_fails_closed) ... ok
test_forced_bessel_failure_is_fail_closed (test_track_b_two_cusp.TrackBTwoCuspTests.test_forced_bessel_failure_is_fail_closed) ... ok
test_identity_reduction_commutes_with_matrix (test_track_b_two_cusp.TrackBTwoCuspTests.test_identity_reduction_commutes_with_matrix) ... ok
test_scaling_and_physical_residual_round_trip (test_track_b_two_cusp.TrackBTwoCuspTests.test_scaling_and_physical_residual_round_trip) ... ok
test_sigma0_direct_and_specialized_actions_agree (test_track_b_two_cusp.TrackBTwoCuspTests.test_sigma0_direct_and_specialized_actions_agree) ... ok
test_spectral_parameter_dependence (test_track_b_two_cusp.TrackBTwoCuspTests.test_spectral_parameter_dependence) ... ok

======================================================================
FAIL: test_closing_artifacts_verify_independently (test_track_b_two_cusp.TrackBTwoCuspTests.test_closing_artifacts_verify_independently)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/workspace/scratch/5332e74d4d90/bianchi-selberg/test_track_b_two_cusp.py", line 148, in test_closing_artifacts_verify_independently
    self.assertTrue(verify_from_paths(result_path, ledger_path)["verified"])
AssertionError: False is not true

----------------------------------------------------------------------
Ran 8 tests in 0.134s

FAILED (failures=1)
```

## Full targeted unittest output: PR base

```text
test_all_four_blocks_and_slow_direct_rows_agree (test_track_b_two_cusp.TrackBTwoCuspTests.test_all_four_blocks_and_slow_direct_rows_agree) ... ok
test_closing_artifacts_verify_independently (test_track_b_two_cusp.TrackBTwoCuspTests.test_closing_artifacts_verify_independently) ... FAIL
test_each_certification_gate_fails_closed (test_track_b_two_cusp.TrackBTwoCuspTests.test_each_certification_gate_fails_closed) ... ok
test_forced_bessel_failure_is_fail_closed (test_track_b_two_cusp.TrackBTwoCuspTests.test_forced_bessel_failure_is_fail_closed) ... ok
test_identity_reduction_commutes_with_matrix (test_track_b_two_cusp.TrackBTwoCuspTests.test_identity_reduction_commutes_with_matrix) ... ok
test_scaling_and_physical_residual_round_trip (test_track_b_two_cusp.TrackBTwoCuspTests.test_scaling_and_physical_residual_round_trip) ... ok
test_sigma0_direct_and_specialized_actions_agree (test_track_b_two_cusp.TrackBTwoCuspTests.test_sigma0_direct_and_specialized_actions_agree) ... ok
test_spectral_parameter_dependence (test_track_b_two_cusp.TrackBTwoCuspTests.test_spectral_parameter_dependence) ... ok

======================================================================
FAIL: test_closing_artifacts_verify_independently (test_track_b_two_cusp.TrackBTwoCuspTests.test_closing_artifacts_verify_independently)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/workspace/scratch/31144417503a/bianchi-base/test_track_b_two_cusp.py", line 148, in test_closing_artifacts_verify_independently
    self.assertTrue(verify_from_paths(result_path, ledger_path)["verified"])
AssertionError: False is not true

----------------------------------------------------------------------
Ran 8 tests in 0.138s

FAILED (failures=1)
```
