# d=19 closure and regression record

For Γ=PSL₂(O₋₁₉), the frozen full trace bound is
**B ∈ [0.89665300, 0.89668143] <1**, excluding discrete Laplace eigenvalues
in (0,1), subject to the documented mathematical/analytic dependencies.
[D19_INVENTORY_PROOF.md](D19_INVENTORY_PROOF.md) and
`groups/d19_inventory.py` prove and replay the complete two-class inventory.
[The certificate](../certificates/d19-k2.json) retains both endpoint balls:

```text
lower [0.8966530054607105611801803788891943976403 +/- 2.36e-41]
upper [0.8966814228100802083070580593100699469567 +/- 2.99e-41]
```

This branch is stacked on the d=11 draft, commit
`56c65c339b16f32584e705160d1dc5a62f285cb8`. The new backend adds d=19
without changing existing analytic formulas or prior inventories. All 152
existing JSON/JSONL artifacts and the historical `old_RIGOR_GAPS.md` are
byte-identical to that base. Existing failing tests are unchanged.

## Proof and independent checks

Both maximal relative orders have class number one; their complete unit
norm images are {±1}. There is one involution class with full finite
centralizer 4 and one order-3 class with centralizer 3. Its inverse merges
by an explicit SL conjugator. Primitive norms are 57799+13260√19 and
45601+6040√57. The base PID proof uses Minkowski, not norm-Euclideanity.
[Published counts](D19_PUBLISHED_COUNTS.md) agree with the inventory.

The no-geodesic claim is load-bearing. The global trace cutoff, exhaustive
proved box and strict Arb comparisons establish systole
acosh((5+√77)/4)=1.907925539233777287698… . The downward-rounded support
is 1.907925539233777056935…; absence inside it follows from that lower
bound, not only an empty search list.

The [search manifest](../certificates/d19-search.json) has 49 screen samples
(k=2,…,8; fractions .8,.9,.95,.975,.99,.999,1; R=40), five refinements
(k=2; fractions .98,.99,.995,.999,1; R=128), then k=2, frac=1, R=256.
The selected upper endpoint is the smallest among these **55 evaluations**,
not a global optimum. δ is binary `0x1.e86dcee236553p-2`, precision 100 bits.

The independent primitive-matrix trace path at 200 bits computes coefficients
0.72861712087236581934… and 1.26898134505847833834…; both balls are
contained in the production 100-bit balls, with radii 9.2104e−31 and
3.0951e−30 respectively. This checks lengths through a different code path,
while relying on the proved class list, orders and centralizer sizes.
The independent volume series sums χ₋₁₉ through N=99997 with an Abel tail,
and applies |D|^(3/2) ζ_K(2)/(4π²). Its enclosure contains production
volume 2.65314813111069769324…; neither ±1% perturbation overlaps it.
Full balls and commands are in `tests/d19-results/coefficients.log` and
`volume.log`.

## Regressions actually run

Environment: Python 3.12.14, python-flint 0.9.0, FLINT/Arb 3.6.0.
Root discovery intentionally covers root-level tests; explicit `tests/`
discovery is a separate suite. The two suites must both be invoked.

| Item | Result | Exact command |
|---|---|---|
| Root unittest suite | FAIL: 68 tests, one existing ledger failure | `python -m unittest discover -v` |
| tests/ suite | FAIL: 54 tests, two existing obsolete incomplete-field expectations | `python -m unittest discover -s tests -v` |
| New d=19 tests | PASS: five tests | `python -m unittest tests.test_d19_inventory -v` |
| Optimized arithmetic | PASS for each d=2,7,11,19 | `for d in 2 7 11 19; do python -O -m groups.d${d}_inventory; done` |
| Global trace lower bounds | PASS | `python -O -m groups.systoles` |
| Frozen d=2 and d=7 B/endpoints | PASS, exact equality | `python tests/mutation/regression_checks.py replay` |
| Frozen d=11 B/endpoints | PASS, exact equality | `python tests/d11_regressions.py replay` |
| Frozen d=19 B/endpoints | PASS, exact equality | `python tests/d19_regressions.py replay` |
| Independent matrix coefficients | PASS | `python tests/d19_regressions.py coefficients` |
| Independent volume | PASS | `python tests/d19_regressions.py volume` |
| Protected artifacts | PASS, 152 unchanged | `python tests/d19_regressions.py protected` |
| Quotient operators | PASS: four tests | `python -m unittest tests.test_quotients -v` |
| Script and module CLI | PASS, full JSON equality to each frozen certificate | Commands below |
| Document links | PASS: local inline paths and external HTTP responses; fragments not checked | `python tests/mutation/regression_checks.py links --external` |
| Whitespace | PASS | `git diff --check` |

An initial quotient command used a nonexistent module,
`python -m unittest tests.test_finite_quotients -v`, and failed with
`ModuleNotFoundError: No module named 'tests.test_finite_quotients'`.
The corrected existing module command above passes; the original output
is retained in `tests/d19-results/quotients-command-error.log`.
The first concurrent tests/ log ended mid-output; the complete sequential
rerun is `tests/d19-results/tests-suite-final.log`. No pass is inferred
from an incomplete log. Its trailing space was stripped for `git diff --check`;
the complete failure output is unmodified.

For CLI reproduction, from root:

```sh
mkdir -p /tmp/d19-regression-cli
for d in 2 7 11 19; do
  python examples/d${d}_certificate.py --output /tmp/d19-regression-cli/d${d}-script.json
  python -m examples.d${d}_certificate --output /tmp/d19-regression-cli/d${d}-module.json
done
python tests/mutation/regression_checks.py cli --output-dir /tmp/d19-regression-cli
python tests/d11_regressions.py cli --output-dir /tmp/d19-regression-cli
python tests/d19_regressions.py cli --output-dir /tmp/d19-regression-cli
```

The exact expanded commands/exits are also in `tests/d19-results/cli-status.json`.
To reproduce the parameter search without modifying frozen artifacts:

```sh
python examples/d19_search.py --output-dir /tmp/d19-search-replay
```

## Failures retained verbatim

Root test `test_closing_artifacts_verify_independently`, line 148:

```text
self.assertTrue(verify_from_paths(result_path, ledger_path)["verified"])
AssertionError: False is not true
```

Triage reads existing, tracked `track_b_two_cusp_result.json` and
`track_b_hejhal_rows.jsonl`. Only the deterministic ledger hash check fails:
actual `8ddb66397840d9ebddca7ac2d4e8810076ec031c0931ea4c21cc0940c66c6486`,
stored `383557121919fcddfb8a3a6b8592238f612eca6b1f6bba3bc4583d15611281a4`.
It is Gaussian level (2+i), not a d=2/7/11/19 artifact.
`python tests/mutation/regression_checks.py triage` prints full verification,
paths, tracked status and discovery counts. The same eight root Track-B
tests were run on the d=11 base in a detached worktree and produce the
same one failure.

The tests/ failures are `test_status_string_cannot_certify_another_field`
(line 61, d=11) and `test_new_fields_block_full_assembly` (line 108,
loop [11,19]); each says:

```text
AssertionError: ValueError not raised
```

Both also fail on the d=11 base. The loop stops at d=11 there and here;
d=19 now also has a valid backend, so its old expectation is obsolete too.
No failing test was edited. Full output is preserved in `tests/d19-results/`.
Base reproduction:

```sh
git worktree add --detach /tmp/bianchi-d19-base 56c65c339b16f32584e705160d1dc5a62f285cb8
(cd /tmp/bianchi-d19-base && python -m unittest test_track_b_two_cusp -v)
(cd /tmp/bianchi-d19-base && python -m unittest tests.test_d2_inventory.D2InventoryTests.test_status_string_cannot_certify_another_field tests.test_trace_core.TraceCoreTests.test_new_fields_block_full_assembly -v)
```

The existing mutation harness results are appended to
[MUTATION_TESTS.md](MUTATION_TESTS.md), with predictions, classifications,
counters and soundness limits. Coefficient values are protected by the
proof note and arithmetic replay, not B<1, which cannot catch undercounts.
The remaining three fields are scoped in
[TRACK_B_LAST_THREE_HANDOFF.md](TRACK_B_LAST_THREE_HANDOFF.md).
