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

## Follow-up: classify each zero-movement mutation

This follow-up instruments the same 10 zero-movement perturbations. New
measurements are in [zero-movement-results.json](../tests/mutation/zero-movement-results.json);
the original 32-case results are retained unchanged. Values and radii below
are from actual Arb objects, not inferred from coarse bound formatting.
Classifications are **1** = vanishes at these support parameters,
**2** = movement below the term's Arb radius, and **3** = numerical target
not reached in the native pipeline. **No case is attributed to class 2:**
none of these measurements establishes a small nonzero effect hidden by
interval uncertainty.

The earlier “below radius” flag was a numerical flag, not a causal diagnosis.
In particular, the eight unchanged translation-mutation B values were
**diagnostic values after arithmetic-gate bypass**, not native certificates.
Removing a primitive translation does not remove its separately stored norm
or its NCE contribution. Those NCE terms are positive and nonzero; they are
not Fourier values at the primitive translation length. Their lengths being
large would not make these elliptic contributions vanish.

### Eight translation-witness mutations: class 3

Each native run reached the inventory verifier once, rejected the missing
witness, and executed **zero coefficient calls**. A diagnostic rerun with
only `require_inventory` bypassed executed the unchanged coefficient formula,
which reads the stored norm and does not read `primitive_translation`.
The requested missing-witness-to-numerical-term dependency is not present.
The missing-witness coefficient counters below show that the diagnostic
really received the mutated records; this does not establish a valid
certificate or a witnessed primitive norm.

| Group | Mutation | Target term | Baseline Arb value | Actual Arb radius | Class | Native inventory / coefficient calls | Diagnostic calls using missing witness |
|---|---|---|---|---|---|---|---|
| d=2 | translation:0:d2 order 2: multiplier A2 | NCE class contribution | `[0.22330651696674861957860859578 +/- 5.72e-30]` | `2.628107700942912017008170743997461889385e-30` | 3 | 1 / 0 | 2 |
| d=2 | translation:1:d2 order 2: multiplier S2 | NCE class contribution | `[0.11165325848337430978930429789 +/- 2.88e-30]` | `1.331799042661808211387130940640445732110e-30` | 3 | 1 / 0 | 2 |
| d=2 | translation:2:d2 order 3: alpha | NCE class contribution | `[0.25813996549286154668341755604 +/- 4.72e-30]` | `3.202093125951373209623528650298175110640e-30` | 3 | 1 / 0 | 2 |
| d=2 | translation:3:d2 order 3: alpha inverse | NCE class contribution | `[0.25813996549286154668341755604 +/- 4.72e-30]` | `3.202093125951373209623528650298175110640e-30` | 3 | 1 / 0 | 2 |
| d=2 | translation_all | NCE | `[0.8512397064358460227347480057 +/- 5.95e-29]` | `1.086773945256828479684271354262546653064e-29` | 3 | 1 / 0 | 8 |
| d=7 | translation:0:d7 order 2: A2 | NCE class contribution | `[0.36450391157636372869396769320 +/- 6.12e-30]` | `2.265924208700321664376410260785994729110e-30` | 3 | 1 / 0 | 2 |
| d=7 | translation:1:d7 order 3: inverse classes merged | NCE class contribution | `[0.36671062069376735078119507415 +/- 6.11e-30]` | `2.593358177754285112055414442760299321000e-30` | 3 | 1 / 0 | 2 |
| d=7 | translation_all | NCE | `[0.73121453227013107947516276735 +/- 4.92e-30]` | `4.859282386454606776431824703546294050109e-30` | 3 | 1 / 0 | 4 |

Verbatim rejection in every one of these eight native runs:

```text
ArithmeticError: translation witness is missing or outside the specified group
```

### Two truncated trace domains: class 1 at the frozen support

The target is the loxodromic/geodesic sum omitted because the Fourier support
lies at/below the proven shortest length. It is **zero by the support proof**,
with exact radius zero; it is not a stored `Evaluation.terms` entry and no
geodesic-sum numerical routine is being claimed. The range monkeypatch did
execute: 12 range replacements and 3 backend systole verification calls per
full evaluate/export run. The surviving shortest trace is unchanged.

| Group | Mutation | Target term | Baseline value / Arb radius | Class | Support radius S=2kδ | Shortest removed trace | Its length | Patched range / verification calls |
|---|---|---|---|---|---|---|---|---|
| d=2 | trace box [-6,6] → [-1,1] | omitted loxodromic sum | 0 / 0 (support-derived) | 1 | `[1.315640939027891809232073683233465999365 +/- 1.48e-40]` | `[-2, -1]` | `[1.76274717403908605046521864996 +/- 2.38e-30]` | 12 / 3 |
| d=7 | trace box [-6,6] → [-1,1] | omitted loxodromic sum | 0 / 0 (support-derived) | 1 | `[1.265948638401894754679233301430940628052 +/- 2.43e-40]` | `[-2, 1]` | `[1.48602212487692709247916872095 +/- 6.59e-30]` | 12 / 3 |

For both groups, the shortest removed length is **provably strictly above**
the support upper endpoint. Thus every removed trace is outside support.
This explains zero numerical movement without asserting that the shortened
enumeration still proves completeness. The full, unmodified systole theorem
is what justifies the frozen certificate.

### Class-1 reruns with an inside-support target

Pre-run prediction: shorten one excluded trace's synthetic length to S/2.
The complete verifier should reject it as shorter than the witness; the
truncated verifier should fail to visit it. These predictions were written
in the harness before execution. This mutates `cosh_length` in memory for
one exact `(A, radicand)` pair, with S and the certificate parameters fixed.
The chosen pair has no representative in the truncated box. For d=7 the
shortest removed pair also occurs in the retained box, so a different
excluded pair was used to isolate non-visitation.

This is a **synthetic inconsistent trace-length fixture**, not a claim that
an actual geodesic in either group has that shorter length. Its Fourier
kernel g(S/2) is provably positive. No orbital weight was invented and no
valid geodesic certificate was exported from a rejected run.

| Group | Target trace | Original length | Mutated length S/2 | g(S/2) |
|---|---|---|---|---|
| d=2 | `[-2, -1]` | `[1.76274717403908605046521864996 +/- 2.38e-30]` | `[0.6578204695139459046160368416167329996824 +/- 2.65e-41]` | `[0.2533619344345034938073001098 +/- 6.73e-29]` |
| d=7 | `[-3, 0]` | `[1.92484730023841378999103565370 +/- 4.68e-30]` | `[0.6329743192009473773396166507154703140259 +/- 2.11e-41]` | `[0.2633071541939694153636436388 +/- 3.91e-29]` |

| Group | Search | Shortened-length branch calls | Range replacements | Result / B |
|---|---|---|---|---|
| d=2 | complete | 1 | 0 | Rejected before B: `ArithmeticError: trace (-2,-1) not proved longer than witness` |
| d=2 | truncated | 0 | 12 | Accepted; exact baseline B unchanged: `[0.43 +/- 5.45e-3]` |
| d=7 | complete | 1 | 0 | Rejected before B: `ArithmeticError: trace (-3,0) not proved longer than witness` |
| d=7 | truncated | 0 | 12 | Accepted; exact baseline B unchanged: `[0.3640 +/- 5.32e-5]` |

The inside-support truncated reruns are **class 3 for the injected short-length
branch**: counters confirm it was never visited. The positive interior
kernel does not become a geometric contribution because the production
pipeline has no geodesic-sum assembly route at these parameters. Consequently
these reruns cannot provide an independently computed nonzero geodesic B
term. They expose the exact reachability limit, without fixes.

## Independent volume check

Computed independently of `fields.quadratic.volume`, `zetaK2`, and production
`character_data` using

\[
V=\frac{|D|^{3/2}\zeta_K(2)}{4\pi^2},\qquad
\zeta_K(2)=\frac{\pi^2}{6}\sum_{n\ge1}\frac{\chi_D(n)}{n^2}.
\]


The checker uses explicit independent residue tables for D=-8 and D=-7.
It sums 100,000 terms for D=-8 and 99,995 terms for D=-7 (whole character
periods), in Arb at 100-bit precision. For both characters, the partial
sums over a period lie between 0 and 2 and the full-period sum is zero.
Abel summation therefore encloses the remaining L-series tail in
`[0, 2/(N+1)^2]`. The volume computation uses integer `|D|*sqrt(|D|)` and
zeta(2)=π²/6, rather than the production Hurwitz-zeta implementation or its
fractional-power call. This independently checks the numerical L-value and
normalization against the requested formula, not the formula's underlying
covolume theorem.

| Group | D | Production volume ball | Independent volume ball | Actual independent Arb radius | Check |
|---|---|---|---|---|---|
| d=2 | -8 | `[1.0038410033411981372723648858 +/- 4.78e-29]` | `[1.003841003 +/- 4.36e-10]` | `9.428045887396874213948194665135815739632e-11` | PASS: contains production ball |
| d=7 | -7 | `[0.8889149278163532635989041542 +/- 2.02e-29]` | `[0.8889149278 +/- 9.36e-11]` | `7.717509953580198311939852828800212591887e-11` | PASS: contains production ball |

Both comparisons overlap and, more strongly, the independent enclosure
contains the entire production enclosure. The +1% and -1% scaled production
volumes are disjoint from the independent enclosure in both groups: this
independent check would reject all four earlier volume perturbations. It
was not part of the native production gate in the original campaign and
does not retroactively change those rows' “none” classification.

Reproduce this follow-up from the repository root:

```sh
python tests/mutation/zero_movement.py /workspace/scratch/31144417503a/regression-logs/zero-movement.json > /workspace/scratch/31144417503a/regression-logs/zero-movement.log 2>&1
```

The run completed successfully for all 10 original zero-movement cases,
four inside-support probes, and two independent volume checks. Original
production code, existing tests, certificates and historical reports remain
unchanged. No flagged behavior was fixed.

## Follow-up: primitive translation consumers and executing mutations

The primitive translations are **per-class matrix witnesses**, not a separate
numerical translation-set collection. Each class separately stores its norm,
finite-centralizer order, and elliptic order. Numerical assembly consumes
those scalars. It does not derive a fresh norm from `primitive_translation`.
The arithmetic proof establishes the relationship between those records.

### Consumers in order

“Arithmetic replay” below means the standalone `verify_inventory()` routine
used by `python -O -m groups.d2_inventory` / `groups.d7_inventory`. The profiled
follow-up called the same routine under ordinary Python; it did not claim
new optimized runs. `GroupData.verify_inventory()` is a **group-bound wrapper**
used by the certificate pipeline; it is different from the standalone
field-inventory routine. Full call-count dictionaries for d=2 and d=7 are in
[translation-consumer-results.json](../tests/mutation/translation-evidence/translation-consumer-results.json).

The order includes two branches: arithmetic binding first, then scalar
numerical assembly. Arithmetic source construction is repeated when replayed.

| Order | Consumer (function, file) | What flows through it | Runs in |
|---|---|---|---|
| 1 | `expected_classes`, `groups/d2_inventory.py` or `groups/d7_inventory.py` | Creates each translation matrix and the separate exact norm record | both |
| 1a | `QuadraticOrder.matrix`, `power`, `QuadraticOrder.mul`, `groups/relative_orders.py` | Builds the matrix from an explicit unit/power; these are source helpers | both |
| 2 | `d2_group` / `d7_group`, `groups/builtins.py` | Wraps the dictionaries into `EllipticClass` records inside `GroupData.elliptic_classes` | certificate pipeline |
| 3 | `GroupRegistry.get`, `groups/registry.py`; `get_group`, `groups/data.py` | Resolves the factory or passes an already constructed group | certificate pipeline |
| 4 | `evaluate`, `core/assemble.py` → `GroupData.require_inventory` → `GroupData.verify_inventory`, `groups/data.py` → `GroupRegistry.verify`, `groups/registry.py` | Sends group records to the inventory gate before B assembly | certificate pipeline |
| 5 | `InventoryBackend.verify` → `_verify_element_witnesses`, `groups/registry.py` | Reads each T directly; checks existence, membership, commutation | certificate pipeline |
| 5a | `GroupKey.contains`, `groups/identity.py`; `MatrixOps.validate/det/mul`, `groups/matrix.py` | Exact matrix shape/determinant/membership and commutator operations used on T | certificate pipeline for membership; shared matrix arithmetic also runs in arithmetic replay |
| 6 | `verify_group_records`, `groups/d2_inventory.py` / `groups/d7_inventory.py` | Invokes the field proof, then compares actual T and all recorded fields with the regenerated expected classes | certificate pipeline |
| 7 | `verify_inventory`, `groups/d2_inventory.py` / `groups/d7_inventory.py` | Runs number-theory/unit completeness checks, then witness verification | both |
| 8 | `verify_class_witnesses`, `groups/d2_inventory.py` | Regenerates T; checks unit matrix matches T, primitive norm, determinant and commutation | both |
| 8 | `verify_witnesses_and_splitting`, `groups/d7_inventory.py` | Regenerates T; checks determinant/commutation and the generator-derived individual norm record | both |
| 8a | `RelativeOrder.abs_square`, field inventory module; `MatrixOps.det/mul`, `groups/matrix.py` (d=2 aliases `mat_det/mat_mul`) | Exact primitive norm and matrix operations used inside those witness checks | both |
| 9 | `verify_group_records`, field inventory module (return phase) | Compares proof class list back to the supplied records, including `primitive_translation` | certificate pipeline |
| 10 | `LevelOneBackend.geometry`, `core/backends/level_one.py` → `GroupData.analytic_data`, `groups/data.py` | Creates geometry and cached `C_ell`; calls scalar coefficient routine, without consuming T | certificate pipeline |
| 11 | `EllipticClass.coefficient`, `groups/data.py` | Reads `norm`, `norm_denominator`, `m`, `finite_centralizer_order`; computes log(N)/(4·centralizer·sin²). Does not read T | certificate pipeline |
| 12 | `evaluate`, `core/assemble.py` | Calls coefficient again directly, multiplies the sum by `g0` to produce NCE, adds NCE to final B | certificate pipeline |
| 13 | `certificate_payload`, `core/certificate.py` | Repeats `InventoryBackend.verify` on the records; repeats geometry/support checks and checks B<1 | certificate pipeline |
| 14 | `report_payload`, `core/certificate.py` | Repeats geometry; serializes B/terms. `certificate_payload` attaches proof class records including T and recomputes the reported elliptic coefficient | certificate pipeline |
| 15 | `main`, `examples/group_certificate.py` | Serializes/writes the certificate payload; both compatibility certificate commands delegate here | certificate pipeline |
| 15a | Inventory module `__main__`, `groups/d2_inventory.py` / `groups/d7_inventory.py` | Serializes/prints the standalone proof class records, including T; no B | arithmetic replay |
| — | Witness-to-norm recomputation inside `EllipticClass.coefficient` | No such code path exists | neither |
| — | Numerical use of cached `geometry.constants['C_ell']` as the assembly NCE coefficient | Assembly recomputes coefficients directly; this cached value is not used by `_terms` | neither |

The standalone proof is not a replay of arbitrary supplied `GroupData`.
Missing records injected into `GroupData` are caught by the group-bound gate;
the standalone routine independently regenerates its own class records.

Observed counts for a fresh group followed by full evaluate/export:

| Function | d=2 certificate | d=7 certificate | Standalone replay, each field |
|---|---|---|---|
| Group factory | 1 | 1 | 0 |
| `expected_classes` | 3 | 3 | 1 |
| Field `verify_inventory` | 2 | 2 | 1 |
| `verify_group_records` | 2 | 2 | 0 |
| `_verify_element_witnesses` | 2 | 2 | 0 |
| `GroupData.analytic_data` | 3 | 3 | 0 |
| `EllipticClass.coefficient` | 20 | 10 | 0 |

### Reruns at consumers that execute

Two placements keep the distinction between record deletion and numerical
sensitivity explicit. All use the frozen k/δ/R and verify baseline B plus
both endpoint strings against the frozen certificate before mutating.

**A. Original witness-removal semantics at the executing guard.** A wrapper
on `_verify_element_witnesses` injects the missing T records immediately
before the original consumer runs. All eight are rejected before B exists;
mutated B and delta are therefore **unavailable, not zero**. No guard bypass
is used in these reruns.

| Group | Mutation | Baseline B | Mutated B / delta | Check | Wrapper / original guard / missing records / downstream coefficient calls |
|---|---|---|---|---|---|
| d=2 | remove d2 order 2: multiplier A2 | [0.42455184, 0.42957479] | unavailable / unavailable | arithmetic replay | 1 / 1 / 1 / 0 |
| d=2 | remove d2 order 2: multiplier S2 | [0.42455184, 0.42957479] | unavailable / unavailable | arithmetic replay | 1 / 1 / 1 / 0 |
| d=2 | remove d2 order 3: alpha | [0.42455184, 0.42957479] | unavailable / unavailable | arithmetic replay | 1 / 1 / 1 / 0 |
| d=2 | remove d2 order 3: alpha inverse | [0.42455184, 0.42957479] | unavailable / unavailable | arithmetic replay | 1 / 1 / 1 / 0 |
| d=2 | remove all translations | [0.42455184, 0.42957479] | unavailable / unavailable | arithmetic replay | 1 / 1 / 4 / 0 |
| d=7 | remove d7 order 2: A2 | [0.36394687, 0.36400056] | unavailable / unavailable | arithmetic replay | 1 / 1 / 1 / 0 |
| d=7 | remove d7 order 3: inverse classes merged | [0.36394687, 0.36400056] | unavailable / unavailable | arithmetic replay | 1 / 1 / 1 / 0 |
| d=7 | remove all translations | [0.36394687, 0.36400056] | unavailable / unavailable | arithmetic replay | 1 / 1 / 2 / 0 |

Every original guard call raises, verbatim:

```text
ArithmeticError: translation witness is missing or outside the specified group
```

**B. Numerical contribution deletion at the executing scalar consumer.**
The eight selected translation targets are identified by their matrix values
at `EllipticClass.coefficient`; the monkeypatch returns zero for their
associated coefficients. This deliberately tests **removal of the associated
NCE contribution**, not removal of a matrix witness with its scalar norm left
intact. The arithmetic records remain valid, so this downstream trusted-code
mutation is not caught by witness/record verification. No production code
was edited, and this is not a claim about valid mutated spectral theorems.

Predictions were recorded in the harness before these runs: original guard
placement rejects; scalar contribution placement decreases B by the selected
positive NCE contribution. All eight numerical runs move downward, as predicted.
The five invocations per class are three cached-geometry evaluations, direct
assembly, and the exporter's reported coefficient calculation. Counters prove
that the suppressed numerical target was actually reached.

| Group | Mutation | Baseline B | Mutated B | Delta | Check rejected | Coefficient calls / suppression hits | Zero movement? |
|---|---|---|---|---|---|---|---|
| d=2 | remove d2 order 2: multiplier A2 | [0.42455184, 0.42957479] | [0.20124532, 0.20626827] | [-0.22832946, -0.21828357] | none | 20 / 5 | no |
| d=2 | remove d2 order 2: multiplier S2 | [0.42455184, 0.42957479] | [0.31289858, 0.31792153] | [-0.11667620, -0.10663031] | none | 20 / 5 | no |
| d=2 | remove d2 order 3: alpha | [0.42455184, 0.42957479] | [0.16641187, 0.17143482] | [-0.26316291, -0.25311702] | none | 20 / 5 | no |
| d=2 | remove d2 order 3: alpha inverse | [0.42455184, 0.42957479] | [0.16641187, 0.17143482] | [-0.26316291, -0.25311702] | none | 20 / 5 | no |
| d=2 | remove all translations | [0.42455184, 0.42957479] | [-0.42668787, -0.42166492] | [-0.85626265, -0.84621676] | none | 20 / 20 | no |
| d=7 | remove d7 order 2: A2 | [0.36394687, 0.36400056] | [-0.00055704, -0.00050335] | [-0.36455759, -0.36445023] | none | 10 / 5 | no |
| d=7 | remove d7 order 3: inverse classes merged | [0.36394687, 0.36400056] | [-0.00276375, -0.00271006] | [-0.36676430, -0.36665694] | none | 10 / 5 | no |
| d=7 | remove all translations | [0.36394687, 0.36400056] | [-0.36726766, -0.36721397] | [-0.73126821, -0.73116085] | none | 10 / 10 | no |

Intervals are rounded outward to eight decimals; exact balls and counter
outputs are retained in the linked JSON and
[consumer-counters.log](../tests/mutation/translation-evidence/consumer-counters.log).
**Zero-movement numerical-consumer mutations: 0/8.** All eight pass the
native B<1 threshold and export, with valid original records but altered
numerical coefficients; no additional independent check is implied.

The earlier eight zero-movement **diagnostic** runs remain in the preceding
section: native coefficient calls 0; diagnostic missing-witness coefficient
calls 2 per single removal, 8 for all d=2 witnesses, and 4 for all d=7 witnesses.
Those diagnostics retain the scalar norms. They are not a native zero-B-change
result and are not silently replaced by these contribution-deletion runs.

### Complete versus truncated length lists up to support

Lists were dumped from the production rational trace-pair and Arb length
functions with the exact original filters (`A<=9`; real |trace|<=2 excluded),
using `[-6,6]^2` and `[-1,1]^2`. They retain trace representatives/multiplicities;
they are not a primitive conjugacy-class geodesic inventory. Inclusion requires
`length.upper() <= support.lower()`; exclusion requires the opposite strict
separation. **Boundary-ambiguous candidates: 0** for both groups and searches.

| Group | Support radius | Complete candidate count | Truncated candidate count | Complete lengths up to support | Truncated lengths up to support | Diff |
|---|---|---|---|---|---|---|
| d=2 | `[1.315640939027891809232073683233465999365 +/- 1.48e-40]` | 18 | 6 | `[]` | `[]` | empty; exit 0 |
| d=7 | `[1.265948638401894754679233301430940628052 +/- 2.43e-40]` | 20 | 6 | `[]` | `[]` | empty; exit 0 |

The support lies strictly below the shortest trace length, including the d=7
fraction=1 case because delta is rounded downward. The empty lists therefore
match. This says nothing about whether the truncated domain proves the
systole; the complete search still has 12 additional candidates for d=2 and
14 for d=7, all outside support.

The four dumps are
[d2 complete](../tests/mutation/translation-evidence/d2-complete-up-to-support.json),
[d2 truncated](../tests/mutation/translation-evidence/d2-truncated-up-to-support.json),
[d7 complete](../tests/mutation/translation-evidence/d7-complete-up-to-support.json),
and [d7 truncated](../tests/mutation/translation-evidence/d7-truncated-up-to-support.json).
Each file contains exactly `[]` plus a newline. Full candidate lengths,
including those excluded by support, are also in translation-consumer-results.json.

Exact commands used from the repository root:

```sh
python tests/mutation/translation_consumers.py /workspace/scratch/31144417503a/translation-consumers > /workspace/scratch/31144417503a/translation-consumers.log 2>&1
diff -u /workspace/scratch/31144417503a/translation-consumers/d2-complete-up-to-support.json /workspace/scratch/31144417503a/translation-consumers/d2-truncated-up-to-support.json
diff -u /workspace/scratch/31144417503a/translation-consumers/d7-complete-up-to-support.json /workspace/scratch/31144417503a/translation-consumers/d7-truncated-up-to-support.json
```

All completed. Production math, certificates, historical reports and
existing tests remain unchanged; no flagged behavior was fixed.
