# d=11 proof and regression record

The level-one certificate is **B in [0.52924365, 0.52928221] <1** for
PSL₂(O₋₁₁). The arithmetic proof is
[D11_INVENTORY_PROOF.md](D11_INVENTORY_PROOF.md); its executable backend is
`groups/d11_inventory.py`. The [frozen certificate](../certificates/d11-k2.json)
retains both endpoint balls. This is an ordinary mathematical proof with
executable checks, dependent on the documented trace formula and Arb;
it is not a proof-assistant formalization.

This extension adds the d=11 inventory and registry entry. Existing shared
analytic formulas, d=2/d=7 inventories, all existing certificates and historical
JSON/JSONL reports are unchanged. All existing test files are also unchanged.
The branch is stacked on the d=7 draft, commit
`bba4a2cd07600c0cb460d9df02ad076ccc5fc5bb`.

## Proof-to-replay map

| Argument | Exact replay | Conclusion |
|---|---|---|
| D11-L1 | QuadraticRing relation and base norm | Euclidean PID, one cusp, no cuspidal elliptics |
| D11-L2 | `verify_orders_and_lattices` | Both relative orders maximal, every small prime principal, class numbers one |
| D11-L3 | `verify_units` | Proved finite boxes; complete μ₄ε₂^Z and μ₆ε₃^Z; determinant images {±1} and {1} |
| D11-L4 | `verify_witnesses_and_splitting` plus L3 | One involution and two distinct inverse order-3 element classes |
| D11-L5 | Exact matrices, subgroup closure, L3 minimality | Primitive norms 199+60√11 and 23+4√33; finite centralizers 4,3,3 |
| Orbital factor | D2-L7 cylinder derivation, coefficient test | log N/(4q sin²(π/m)), including the involution flip |
| No geodesic in support | SYSTOLES global trace cutoff and exact finite minimum | Support 1.53439443650263873664… below proved systole 1.53439443650263888666… |
| Analytic bound | Shared certified quadrature and tails | Frozen B upper endpoint below 1 |

The geodesic claim uses the rigorous trace lower bound and proved finite box
in [SYSTOLES.md](SYSTOLES.md), not merely an empty length list.

The [search manifest](../certificates/d11-search.json) records 49 screen
samples (k=2,…,8; fractions .8,.9,.95,.975,.99,.999,1; R=40), five refinements
(k=2; fractions .98,.99,.995,.999,1; R=128), and the final (k=2, frac=1,
R=256). The final certified upper endpoint is the smallest among these
**55 evaluations**, not a global optimum over continuous parameters or test
functions. δ is binary `0x1.88ce12e3f1746p-2`, rounded down by the existing engine.

## Independent numerical checks and soundness

`tests/d11_regressions.py coefficients` computes lengths from each primitive
translation matrix's trace using
cosh ℓ=(|tr T|²+|(tr T)²−4|)/4. It uses neither the stored norm radical nor
`EllipticClass.coefficient`. At 200-bit precision, each independent coefficient
ball is contained in the production 100-bit ball:

| Class | Coefficient midpoint (rounded) | Production Arb radius |
|---|---:|---:|
| Involution with flip | 0.37415285576579761224 | 4.5895627333e−31 |
| Order 3 | 0.42535205237034455722 | 9.3352430623e−31 |
| Order 3 inverse | 0.42535205237034455722 | 9.3352430623e−31 |

[Full coefficient output](../tests/d11-results/coefficients.log) records the
balls and differences. This check independently exercises the translation-to-
length calculation; it still depends on the proved class list, m and q. It
is not an independent recomputation of the parabolic/scattering coefficients.

The independent volume path sums an explicit χ₋₁₁ character table through
N=99990 and bounds its tail by Abel summation, using prefix bounds [0,3].
With ζ(2)=π²/6 it evaluates |D|^(3/2)ζ_K(2)/(4π²). Its interval
[1.382608308 +/- 2.50e−10] contains the production ball
[1.3826083079026458736716533445 +/- 4.58e−29]; neither ±1% perturbation
of the production volume overlaps it. See the
[volume output](../tests/d11-results/volume.log).

Dropping a positive contribution makes B smaller and can leave B<1 true.
Consequently B<1 does not protect against undercounts. Completeness guards
reject removed witnesses and records; the proof note, record-bound arithmetic
replay and independent coefficient tests protect coefficient values. Changing
trusted coefficient code requires reviewing those derivations and tests.
The new test exercises 16 record perturbations (four alterations and one
class removal per class, plus removing the involution flip); all are rejected
before a bound is evaluated. This is not a new full mutation campaign.
The earlier d=2/d=7 mutation classifications remain in
[MUTATION_TESTS.md](MUTATION_TESTS.md), unchanged.

## Commands actually run

Python 3.12.14; python-flint 0.9.0; FLINT 3.6.0 (including Arb).
Commands below run from the repository root. `$out` denotes the actual directory
`/workspace/scratch/31144417503a/d11-exploration/regressions`; a new directory
may be used for reproduction. Commands producing reports write new output
paths, leaving frozen files untouched. The initial search was run with
`--output-dir certificates` to create the new d=11 artifacts only.

| Check | Result | Exact command |
|---|---|---|
| Complete parameter comparison | PASS, 55 evaluations | `python examples/d11_search.py --output-dir certificates` |
| New d=11 tests | PASS, 5 tests | `python -m unittest tests.test_d11_inventory -v` |
| tests/ discovery | FAIL, 49 tests, 2 failures below | `python -m unittest discover -s tests -v` |
| Root discovery | FAIL, 68 tests, existing Track-B failure | `python -m unittest discover -v` |
| d=11 optimized arithmetic | PASS | `python -O -m groups.d11_inventory` |
| d=2 optimized arithmetic | PASS | `python -O -m groups.d2_inventory` |
| d=7 optimized arithmetic | PASS | `python -O -m groups.d7_inventory` |
| Global systoles | PASS, all supported fields | `python -O -m groups.systoles` |
| Frozen d=11 B and both endpoint balls | PASS, exact strings | `python -O tests/d11_regressions.py replay` |
| Frozen d=2/d=7 B and both endpoint balls | PASS, exact strings | `python -O tests/mutation/regression_checks.py replay` |
| Matrix-trace coefficients | PASS, production balls contain independent balls | `python tests/d11_regressions.py coefficients` |
| Independent volume | PASS | `python tests/d11_regressions.py volume` |
| d=11 base artifacts | PASS, 144 protected files unchanged | `python tests/d11_regressions.py protected` |
| Earlier historical artifacts | PASS, all 134 original JSON/JSONL unchanged | `python tests/mutation/regression_checks.py protected` |
| Quotient operators | PASS, 4 tests | `python -m unittest tests.test_quotients -v` |
| d=11 script | PASS | `python examples/d11_certificate.py --output $out/d11-script.json` |
| d=11 module | PASS | `python -m examples.d11_certificate --output $out/d11-module.json` |
| d=2 script | PASS | `python examples/d2_certificate.py --output $out/d2-script.json` |
| d=2 module | PASS | `python -m examples.d2_certificate --output $out/d2-module.json` |
| d=7 script | PASS | `python examples/d7_certificate.py --output $out/d7-script.json` |
| d=7 module | PASS | `python -m examples.d7_certificate --output $out/d7-module.json` |
| d=11 CLI complete JSON equality | PASS | `python tests/d11_regressions.py cli --output-dir $out` |
| d=2/d=7 CLI complete JSON equality | PASS | `python tests/mutation/regression_checks.py cli --output-dir $out` |
| Markdown file links | PASS, local file targets and all external URLs | `python tests/mutation/regression_checks.py links --external` |
| Whitespace | PASS | `git diff --check` |

Saved outputs are in [tests/d11-results](../tests/d11-results).
For safe search reproduction use
`python examples/d11_search.py --output-dir /tmp/d11-search`.

## Failures preserved verbatim

The tests/ suite has two **new compatibility failures**. Both tests use d=11
as an unproved group. Registering the completed inventory invalidates that
fixture assumption: setting its status to `self-contained` no longer forges
a proof, and complete d=11 assembly is now admitted. Both tests pass on the
exact stacked base; neither test was edited. This branch is **not all green**.
The full suite output is [tests-final.log](../tests/d11-results/tests-final.log).

```text
FAIL: test_status_string_cannot_certify_another_field (test_d2_inventory.D2InventoryTests.test_status_string_cannot_certify_another_field)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/workspace/scratch/5332e74d4d90/bianchi-selberg/tests/test_d2_inventory.py", line 61, in test_status_string_cannot_certify_another_field
    with self.assertRaisesRegex(ValueError, "no self-contained inventory verifier"):
         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AssertionError: ValueError not raised

======================================================================
FAIL: test_new_fields_block_full_assembly (test_trace_core.TraceCoreTests.test_new_fields_block_full_assembly)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/workspace/scratch/5332e74d4d90/bianchi-selberg/tests/test_trace_core.py", line 108, in test_new_fields_block_full_assembly
    with self.assertRaisesRegex(ValueError, "inventory is incomplete"):
         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AssertionError: ValueError not raised
```

The root-level Track-B failure remains pre-existing, unrelated to the d=11
certificate. `verify_from_paths` fails only the ledger-hash comparison for
`track_b_two_cusp_result.json` and `track_b_hejhal_rows.jsonl`, both existing
and git-tracked. Actual ledger SHA-256:
`8ddb66397840d9ebddca7ac2d4e8810076ec031c0931ea4c21cc0940c66c6486`;
stored SHA-256:
`383557121919fcddfb8a3a6b8592238f612eca6b1f6bba3bc4583d15611281a4`.

```text
FAIL: test_closing_artifacts_verify_independently (test_track_b_two_cusp.TrackBTwoCuspTests.test_closing_artifacts_verify_independently)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/workspace/scratch/5332e74d4d90/bianchi-selberg/test_track_b_two_cusp.py", line 148, in test_closing_artifacts_verify_independently
    self.assertTrue(verify_from_paths(result_path, ledger_path)["verified"])
AssertionError: False is not true
```

The root suite intentionally discovers root-level test files; tests/ discovery
runs the trace-engine tests separately. There is no claim that root discovery
includes all tests/ files. Triage command:
`python tests/mutation/regression_checks.py triage`.

The exact d=11 base was checked in a separate worktree:

```sh
git worktree add --detach /workspace/scratch/31144417503a/bianchi-d11-base bba4a2cd07600c0cb460d9df02ad076ccc5fc5bb
cd /workspace/scratch/31144417503a/bianchi-d11-base
python -m unittest test_track_b_two_cusp -v
python -m unittest tests.test_d2_inventory.D2InventoryTests.test_status_string_cannot_certify_another_field tests.test_trace_core.TraceCoreTests.test_new_fields_block_full_assembly -v
```

Track-B fails on base too (8 tests, 1 failure); the two compatibility tests
pass there (2 tests). Full outputs:
[base Track-B](../tests/d11-results/base-track-b.log),
[base compatibility](../tests/d11-results/base-compatibility.log).
No existing failure, certificate, ledger or historical report was repaired.
