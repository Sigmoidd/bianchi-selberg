# Handoff: Track B for d=43,67,163

Target: for each full level-one Γ=PSL₂(O₋d), prove that the discrete
Laplace spectrum contains no eigenvalue in (0,1), with independently
replayable evidence. None of these three fields is closed by this handoff.
Start with d=43, then reuse only field-independent arguments for 67 and 163.

The completed trace-formula fields are d=2,7,11,19; their certificates
and regression commands are in the corresponding review guides. This
handoff does not launch a computation for the remaining three fields.
The existing mechanical trace screens at k=2, frac=.999 fail to put B<1;
that observation is not an impossibility theorem for all trace tests.

## What Track B currently provides

Read `TRACK_B_GLOBAL_HEJHAL_DEFECT_AUDIT.md`,
`TRACK_B_TWO_CUSP_HEJHAL_CERTIFICATE.md`, `track_b_theorem_defect.py`, and
`theorem_DK_sixcopy.tex` first. The current verified finite physical solve
is for a Gaussian congruence group, level (2+i), index six and two cusps.
Its geometry and numerical artifacts are not for any of the three targets.
The theorem is a working draft. The audit records a rank failure at cutoff
M=4, order-one drift of low coefficients, incompatible floor/mass witnesses,
and missing continuum residual ledgers; `rung4_certified` is false.

The finite Hejhal solve checks finitely many rows. Even a proved small
continuum defect with positive mass localizes spectrum near the trial
parameter; it does not alone exclude every eigenvalue in (0,1). Closing
the stated target needs a proved exhaustive lower-bound/counting criterion.
The separate `independent_exclusion/` route contains a compact-core cusp
reduction, Crouzeix–Raviart lower-bound theory and verified matrix positivity.
Read `DESIGN.md`, `lower_bound_theory.md`, `PROOF.md`, and `m3_certify.py`
in that directory to assess this route. Its existing Gaussian geometry and
constants also require new field proofs; no Gaussian certificate transfers.

## Field inputs and first deliverable

| d | Integral relation for τ=(1+√−d)/2 | D | Cusp covolume of O | Current registered inventory |
|---:|---|---:|---|---|
| 43 | τ²=τ−11 | −43 | √43/2 | incomplete |
| 67 | τ²=τ−17 | −67 | √67/2 | incomplete |
| 163 | τ²=τ−41 | −163 | √163/2 | incomplete |

`fields/arithmetic.py`, `groups/arithmetic.py`, `groups/systoles.py` and
`groups/identity.py` already supply field arithmetic, class-number checks,
trace minima and exact group identity. The fields have class number one
and units ±1; there is one level-one cusp with no rotational quotient.
They are not norm-Euclidean. Obtain arithmetic/cusp data from proved PID
and lattice arguments, never from Gaussian generators or nearest-point
division assumptions.

For d=43, first write `docs/D43_TRACK_B_DESIGN.md`: specify exact group,
metric and Laplacian convention, compact core and cusp truncation, target
spectral interval, theorem dependencies and the final exclusion criterion.
Freeze this before expensive numerical work. Choose the exclusion method
explicitly; an approximate eigenvalue search is a diagnostic only.

## Work packages and acceptance evidence

| Order | Work | Acceptance evidence | Existing starting points |
|---:|---|---|---|
| 1 | Prove fundamental domain, all face pairings, stabilizers, cusp lattice and dual frequencies | Exact matrices and relation checks; coverage/nonoverlap proof; precise cusp invariance and correct metric weights | `groups/matrix.py`, `groups/arithmetic.py`; Gaussian `track_b_two_cusp_data.py` as interface example only |
| 2 | Establish field-specific compact-core exclusion theorem | Cusp zero-mode reduction, nonzero-mode lower bounds, constant-mode handling and essential-spectrum threshold 1; explicit coercivity criterion | `independent_exclusion/DESIGN.md`, `lower_bound_theory.md` |
| 3 | Assemble a field-specific mesh and matrices | Validated geometry/quadrature, boundary identification, interpolation/error constants and finite-dimensional lower bounds | `independent_exclusion/cr_prototype.py`, `m3_certify.py` |
| 4 | Use spectral/Hejhal probes to guide refinement, if useful | Correct O-dual Whittaker modes, normalization, cutoff convergence and resolved Arb Bessel values | `track_b_two_cusp_hejhal.py`, `track_b_cutoff_ladder.py`; no reuse of six-copy glue |
| 5 | Certify exclusion over the entire required parameter domain | Verified positivity/lower bounds on every window including near 0 and 1; justified cusp/infinite tails and discretization error | `independent_exclusion/m3_certify.py`; `track_b_global_partition_arb.py` for interval bookkeeping patterns |
| 6 | Independently replay final theorem evidence | Reconstruct matrices/geometry rather than trusting serialized entries; all hashes and field/level/trial identities match; report every failed condition | `track_b_two_cusp_verify.py`, `track_b_global_partition_verify.py` as verifier examples |
| 7 | Audit, freeze, then carry method to 67 and 163 | Per-field artifact, proof note, commands, honest pass/fail and reused-versus-new argument map | Existing `tests/mutation/` for applicable shared checks; field-specific negative tests for new Track-B gates |

If using the continuum-defect approach, additionally prove gluing regularity,
all-face residual bounds, cusp decay/tails and positive projected mass for
exactly the same frozen trial, parameter and cutoff. A pointwise residual
or a stable midpoint vector is insufficient. Connect that result to an
exhaustive exclusion/counting theorem before claiming the target.

Do not force the trace mutation harness onto a different theorem: its B is
a trace bound. Retain it for shared trace code; new Track-B completeness,
geometry, mesh and positivity guards need their own sensitivity checks.

## Baseline and preservation commands

From repository root, install `requirements-trace.txt`, then record:

```sh
python -m unittest discover -v
python -m unittest discover -s tests -v
python -m unittest test_track_b_two_cusp -v
python tests/mutation/regression_checks.py triage
python -O -m groups.systoles
python -O -m groups.d19_inventory
python tests/d19_regressions.py replay
python tests/d19_regressions.py protected
python tests/mutation/regression_checks.py links --external
git diff --check
```

The existing root-level Track-B artifact test has a ledger SHA mismatch,
also failing on the base branch. It reads tracked `track_b_two_cusp_result.json`
and `track_b_hejhal_rows.jsonl`, both Gaussian level (2+i) artifacts. Preserve
these files and the failing test; do not recompute a stored hash merely to
make it pass. Two tests under `tests/` still assume d=11/d=19 cannot assemble
and now fail because these fields have complete inventories. Their exact
failure logs remain part of the d=19 regression record.

New work should use separate branches and per-field artifact directories.
Keep all existing certificates, historical reports and known failing tests
unchanged. Report blockers verbatim; continue independent checks. A failure
of a frozen baseline replay gates mutations. Stop any theorem claim at the
first unproved input even if a numerical approximation looks favorable.
