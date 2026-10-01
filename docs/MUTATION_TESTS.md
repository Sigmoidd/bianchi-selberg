# Regression gate and mutation tests

## Run status

Stopped at the first regression failure, as instructed. No mutation was
applied or run. No production math, certificate, or historical report was
edited. This report does not claim that any certificate was replayed or that
protected files were compared against an earlier revision.

Checkout tested: `codex/d7-certificate`, commit
`d3610366047b3e178942e2c7b9874ca51630b8b5`.

## Regression results

| Item | Result | Exact command used / scope |
|---|---|---|
| 1. Full unittest suite | FAIL; stopped during root discovery | `python -m unittest discover -v` — 68 tests, one failure. Separate discovery under `tests/` was not run, so full-suite coverage was not completed. |
| 2. Exact arithmetic replay under `python -O` | NOT RUN | Stopped at item 1. |
| 3. Frozen d=7 and d=2 v2 replays; endpoint preservation; v1 and historical byte comparisons | NOT RUN | Stopped at item 1. |
| 4. Quotient-operator checks | NOT RUN | Stopped at item 1. |
| 5. d=2 and d=7 CLI as script and module | NOT RUN | Stopped at item 1. |
| 6. Documentation links and `git diff --check` | NOT RUN | Stopped at item 1. |

The first command exited with status 1. Its failure block and summary are
reproduced verbatim:

```text
======================================================================
FAIL: test_closing_artifacts_verify_independently (test_track_b_two_cusp.TrackBTwoCuspTests.test_closing_artifacts_verify_independently)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/workspace/scratch/5332e74d4d90/bianchi-selberg/test_track_b_two_cusp.py", line 148, in test_closing_artifacts_verify_independently
    self.assertTrue(verify_from_paths(result_path, ledger_path)["verified"])
AssertionError: False is not true

----------------------------------------------------------------------
Ran 68 tests in 0.970s

FAILED (failures=1)
```

## Mutation results and predictions

Mutation work was not started because the regression gate failed. No
pre-run predictions were authored, no baseline or mutated B intervals were
measured, and no deltas or catching checks were assessed.

| Requested mutation | Groups | Status |
|---|---|---|
| Drop each elliptic class in turn | d=2, d=7 | NOT RUN |
| Scale centralizer normalization by 2 and by 1/2 | d=2, d=7 | NOT RUN |
| Remove one primitive translation / shrink the translation set | d=2, d=7 | NOT RUN |
| Perturb volume by +1% and -1% | d=2, d=7 | NOT RUN |
| Truncate geodesic/loxodromic enumeration below the proven completeness bound | d=2, d=7 | NOT RUN |
| Omit cusp/scattering block; separately flip its sign | d=2, d=7 | NOT RUN |

All requested flags remain **unassessed**: movement below the Arb radius,
movement opposite to prediction, and uncaught mutations still yielding
B < 1. Absence of measured flags is not a passing result.

## Reproduce the stopping failure

From the repository root at the tested commit:

```sh
python -m unittest discover -v
```

No mutation reproduction command exists because no mutation harness was
created or executed. The failure was not diagnosed or fixed in this run.
