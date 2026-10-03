# Swarm handoff: complete the class-number-one spectral certifications

## Launch instruction

Read this handoff and the cited sources. Launch one coordinator and six workers
with the ownership below. Preserve published checkpoints, reconstruct pruned
inputs, and complete each proof dependency rather than stopping at numerical
probes. Start with d=67 implementation, parallel historical audits and field
inputs. Publish compact evidence at every completed phase. Keep all unresolved
spectral flags false. Report and resolve actual mathematical blockers; do not
weaken the theorem, alter the target group or mask failed replay evidence.

## Mandate and starting point

For each full level-one group Gamma=PSL2(O_-d), establish that the discrete
Laplace spectrum has no eigenvalue in (0,1), with independently replayable
evidence. The class-number-one field list is d=1,2,3,7,11,19,43,67,163.
Use the hyperbolic metric and Laplacian convention stated in each proof.
Do not silently substitute a congruence subgroup, character, or finite cover.

Repository: https://github.com/Sigmoidd/bianchi-selberg

Resume branch: `track-b/d67-spectral-exclusion`.
Published code checkpoint before this handoff:
`299375b908b7a90e74b2f50da39dd8bc70b50826`.
Merged spectral checkpoint PR: https://github.com/Sigmoidd/bianchi-selberg/pull/10
Draft swarm handoff continuation: https://github.com/Sigmoidd/bianchi-selberg/pull/11
Geometry PR: https://github.com/Sigmoidd/bianchi-selberg/pull/9
Use the latest commit containing this handoff; inspect changes since the
pinned checkpoint before using an existing ledger.

The immediate proof implementation priority is d=67. Begin the d=43/d=163
field-input work in parallel, but transfer numerical certificates only after
each field's dependencies are separately proved. Audit the historical d=1/d=3
claims in parallel. Existing d=2,7,11,19 certificates must remain unchanged.
The deliverable is a proof and replay command for each field, not a collection
of approximate eigenvalues or manually edited status flags.

## Current evidence inventory

| d | Evidence at the pinned checkpoint | Required next action |
|---:|---|---|
| 1 | Historical Picard trace and independent FEM exclusions are claimed; registered inventory is legacy | Audit each theorem separately. The FEM route invokes the defective general triangle-variance lemma; repair or prove its particular inclusion argument. Modern trace exporter rejects the legacy inventory. |
| 2 | Self-contained arithmetic inventory and frozen full certificate | Replay and preserve. |
| 3 | Historical numerical enclosure; legacy normalization is `unresolved-flips-frozen` | Resolve full integral centralizers, maximal finite subgroup, primitive norm and conjugacy multiplicity, then recertify. A Klein-four witness alone is insufficient. |
| 7 | Self-contained arithmetic inventory and frozen full certificate | Replay and preserve. |
| 11 | Self-contained arithmetic inventory and frozen full certificate | Replay and preserve. |
| 19 | Self-contained arithmetic inventory and frozen full certificate | Replay and preserve. |
| 43 | Registered trace inventory incomplete; no final field-specific exclusion identified | Prove field inputs and a complete exclusion route. |
| 67 | Exact Ford geometry, adaptive moment map, rank/kernel checks and conditional scalar budgets | Close coefficient binding, full finite enclosures, global index and independent replay. |
| 163 | Registered trace inventory incomplete; no final field-specific exclusion identified | Prove field inputs and a complete exclusion route. |

Read `RIGOR_GAPS.md`, the current registry/source gates, and the field review
guides. README contains stale d=19-open prose and inconsistent d=3 theorem
claims. Resolve these against evidence; do not use a headline as a proof input.
For d=3 read `docs/NORMALIZATION_ISSUE.md` and `bianchi_omega_arb.py`.

## d=67 theorem and remaining blockers

Read, in this order:

1. `docs/D67_TRACK_B_GEOMETRY.md` and `docs/D67_TRACK_B_THRESHOLD_THEORY.md`.
2. `docs/D67_TRACK_B_DISCRETE_INDEX.md`.
3. `docs/D67_TRACK_B_CENTROID_CR_BOUND.md` and `docs/D67_TRACK_B_VERTICAL_ENVELOPE.md`.
4. `docs/D67_TRACK_B_BALANCED_ENVELOPE.md` and `docs/D67_TRACK_B_SUBCELL_MASS.md`.
5. `artifacts/track_b/d67/spectral/README.md`, `moment_map_verification.json`,
   `moment_kernel_check.json`, and `balanced_eta13_budget.json`.

For the normalized glued core with top Y=2, the threshold form is

\[
C(v)=Q_K(v)/A-M_K(v)/A+\tfrac12b(v)-\tfrac34t(v)^2,
\qquad A=\sqrt{67}/2.
\]

The continuum theorem excludes all residual and cuspidal discrete eigenvalues
in (0,1) if the negative index of C is at most one. On the verified CR space,

\[
C_h=Q_h+\tfrac12b_h-\tfrac34t_h^2-(1+\eta)M_h^{sharp},
\qquad C(v)\ge C_h(Iv)+c_0 Q_h(v-Iv),
\]

where c0=1-(1+1/eta)gamma². The exact inherited balanced bound
gamma²<0.070073378286 permits eta=1/13 and mass inflation 14/13:
c0>0.0189727 (approximately 0.01897270399955813). This is conditional on using the same exact original
coefficient intervals and h0 as its scalar ledger. The direct index criterion
requires neither a weighted bulk-mean enclosure nor a sigma penalty.

It is sufficient to certify index_-(C_h)<=1. One possible witness is
C_h+alpha*z*z^T PSD, for explicit rational alpha>0 and rational z.
A candidate is z=(1+eta)M_h^sharp*1+t_h/4, alpha=1/2. Other rational witnesses
or a directly verified inertia bound are permitted; prove the implication.

Four gates remain:

| Gate | Acceptance evidence |
|---|---|
| Coefficient/scalar identity | Reconstruct every leaf's coefficient interval, h0, lower energy matrix and original parent error-mass bound. Use the inherited exact hull bounds, or independently certify the scalar budget for the new envelopes. |
| Finite form enclosure | Independently enclose all energy, sharpened mass and top-fragment contributions; propagate through verified P; bind the same finite space, eta and rational stabilization vector. |
| Global negative index | Rigorous PSD/inertia evidence for the actual enclosed finite form, including permutation, arithmetic error and every pivot/Schur-complement bound. |
| Independent full replay | Reconstruct geometry, topology, coefficients, assembly and index from fresh inputs; exercise sensitivity controls; issue a field/group-specific theorem result only after every gate passes. |

`adaptive_matrix_probe.py` is a floating diagnostic. Its height maximum uses
coordinate-box containment for the stationary center; the inherited exact
scalar ledger uses projected-hull containment. The resulting h0 can differ.
Do not certify the diagnostic matrix merely by attaching the inherited gamma
report. Exact moment rank, a constant-only energy kernel, and exact subcell
Gram conservation do not prove finite threshold positivity.

## Swarm organization: coordinator plus six workers

Use isolated git worktrees and disjoint file ownership. Publish an interface
contract before concurrent implementation. Producer and independent verifier
must not share coefficient/assembly or inertia code that they purport to check.

| Worker | Ownership and task | Required deliverable |
|---|---|---|
| A: theorem auditor | New proof notes; threshold/index transfer; class-number-one scope | Audit exhaustive interval coverage, constants, quotient traces and candidate stabilization. Generalize cusp/truncation coefficients before transferring fields. Review all worker claims. |
| B: finite-form producer | New d=67 coefficient/enclosure producer and ledgers | Exact rational or rigorous interval lower energy/upper mass forms on all 909,276 leaves. Bind scalar identity. Preserve original interpolation-error mass bounds while sharpening finite mass. |
| C: global index certificate | New solver/certificate format and resource plan | Prove index<=1 with explicit arithmetic-error treatment. Evaluate ordering, block/streamed Schur complements, domain decomposition or a rigorously bounded iterative method. No expensive factorization before resource admission. |
| D: independent replay | Separate verifier and new negative controls | Reconstruct B/C evidence without importing their mathematical producer routines; reject malformed, mismatched, incomplete and overoptimistic evidence. Final theorem command must fail closed. |
| E: historical theorem audit | New d=1/d=3 audit branches; normalization/inventory proofs | Reconcile d=1 trace versus FEM routes; resolve d=3 normalization with exact group witnesses and a self-contained inventory. Preserve historical artifacts and record corrections separately. |
| F: field inputs | Separate d=43 and d=163 branches/directories | Freeze each design; exact fundamental domain, pairings, stabilizers, cusp lattice/dual minimum, truncation and scalar bounds. Complete d=43 input chain, then d=163, and use audited shared theorem interfaces. |

Coordinator owns integration, the all-nine manifest, baseline protection,
dependency admission and PR descriptions. Run only one large d=67 matrix or
factorization job per 8-GiB execution host. Lightweight proof/arithmetic work
may proceed concurrently. More workers must not multiply large arrays in RAM.

Worker B can use the following exact construction for local finite mass:
subdivide only for integration; for parent CR basis phi and subcell vertex
value vectors f_i,

\[
\int_S\phi\phi^T=\frac{|S|}{20}
[(\sum_i f_i)(\sum_i f_i)^T+\sum_i f_i f_i^T].
\]

Use concavity to bound H below by its subcell vertex minimum, bound r below
by its vertex minimum, and take m_S=(4-H_min)/(2g_min²),
g_min=(1-r_min)H_min+4r_min. Upward dyadic rounding of each density is valid;
cap it by the independently proved parent upper bound if useful. Basis values
are exact dyadics for red subdivision. Prove matrix-order mass domination,
not entrywise domination: individual basis products may be negative.
This sharpens M_h without changing P or the original interpolation-error
gamma. Independently test exact constant-density Gram conservation and mass
PSD/domination. A sampled ledger is not a complete finite enclosure ledger.

## Durable evidence and interrupted diagnostics

Workspace maintenance removed the previous transient checkout, raw matrices,
dependencies and a running single-mode diagnostic. The GitHub checkpoint
survived and has been re-cloned. Do not assume those raw files remain anywhere.
Restore geometry/topology by replay from the published sources and frozen plan.

Completed process outputs observed before maintenance, preserved separately
in `artifacts/track_b/d67/spectral/pre_handoff_observations.json`:

- Depth-three floating finite mass total: 3.9009611474015085; depth two:
  4.143787464796889; original parent bound: 5.75509150845232.
- Depth-three eta=1/13 seven-mode run: minimum estimate
  +0.0020864741982754556, residual 0.002068555513554218, 43 history entries.
  This is not a lower bound and does not establish positivity.
- A later 100-iteration single-mode run was interrupted. Its last visible
  iteration-21 estimate was +0.00148711, residual 0.00143669. No completed
  report or certificate is available. Do not use it as a pass.
- A floating complete boxed scalar diagnostic found gamma² approximately
  0.06911955763127157. It does not resolve the exact binding gate.
- An AMD/QDLDL symbolic probe counted 409,801,053 strict lower factor entries
  for the sparse union graph. Int64 indices plus double values require
  6,556,816,848 bytes; int64 plus 128-bit point values require 9,835,225,272
  bytes, before matrices, workspace or interval radii. Treat these as resource
  diagnostics, not inertia evidence; regenerate and bind the graph/ordering.

The interrupted host had 8 GiB RAM and 32 GiB disk. A direct in-core rigorous
factorization with that ordering has no demonstrated memory admission.
Investigate reduced fill, compact indices and rigorous block/streamed methods
before requesting a bigger run. If another resource is necessary, present
the concrete measured requirement; do not silently launch paid compute.

All observed values above are recovered from conversation/process output;
raw binary hashes were not retained. They guide reproduction, not proof.

## Frozen moment bindings

| Artifact | SHA-256 |
|---|---|
| Exact moment input | `71276662d562e35a6508bfce9f00e40ffe8da9ddd0029ae9a0fb45f1ff6a84b6` |
| Compressed adaptive plan | `8d17dd752d3c58ba2b18c850690f3bbf8f61840497a5483967fb3781d53cee18` |
| Topology binary | `5778bae2f5c4358c85e0efb37b933736c5c26d62c46b925e9b9b39353b5e59da` |
| Moment rows binary | `aff22736c1d8d356927ece7e9fae8d0d2b82583bd1a5c389e86cf5f5db34a5f1` |
| Moment producer source | `c89d5fc8d1c9ece13a566821551e47799557a4c4e2cd1d7d5e903861637b397f` |
| Moment verifier source | `a3b7b6555ff58a1defe8ded34ee9451185bc2dfeb95537d86c7b90202bfbdc5a` |
| Original scalar report | `5486f73bc02fba040393cdd14fd9f757bf6a06dfd9f1869d7dc9dc03ae0cd50a` |

Independently reproduce 909,276 leaves, 2,042,568 master variables,
4,084,872 prolongation entries, maximum row support 994, and one final
constant-energy component. Frozen r levels come from the plan; never generate
new floating geometric layers and call them the same proof input.

## Restart and baseline commands

Install `requirements-trace.txt` for its pinned python-flint/mpmath packages.
It does not install NumPy, SciPy or PyAMG. Record versions for the diagnostic
environment; PyAMG 5.2.1 was used previously. The binary moment pipeline
requires a little-endian GCC-compatible platform with signed 128-bit integers.

Before moment regeneration, assembly or factorization, record RAM/disk admission
including sparse-product temporaries, AMG hierarchy and simultaneous arrays.
Do not overlap large jobs. Run from repository root; the following Bash block
stops on failure and records a session log. Also record individual exit codes
in the phase manifest:

```bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
S=artifacts/track_b/d67/spectral
D67_RUN=$(mktemp -d /tmp/d67-swarm-replay.XXXXXX)
exec > >(tee "$D67_RUN/session.log") 2>&1
df -h "$D67_RUN"
python -B artifacts/track_b/d67/verify_geometry.py
python -B -O artifacts/track_b/d67/verify_geometry.py
python -B "$S/verify_reference_mesh.py"
python -B -O "$S/verify_reference_mesh.py"
python -B "$S/verify_adaptive_plan.py"
python -B -O "$S/verify_adaptive_plan.py"
python -B -O "$S/replay_moment_map.py" "$D67_RUN" --rebuild --output "$D67_RUN/moment-replay.json"
python -B "$S/moment_negative_checks.py" "$D67_RUN" --output "$D67_RUN/moment-negative-checks.json"
g++ -std=c++17 -O2 "$S/moment_kernel_check.cpp" -o "$D67_RUN/kernel-check"
"$D67_RUN/kernel-check" "$D67_RUN/full.rows.bin" > "$D67_RUN/kernel.json"
python -B - "$D67_RUN/kernel.json" <<'PY_CHECK'
import json, sys
if json.load(open(sys.argv[1]))["constant_energy_kernel_proved"] is not True:
    raise SystemExit("FAIL: energy kernel is not proved constant-only")
PY_CHECK
python -B -O "$S/balanced_envelope.py" > "$D67_RUN/eta10-budget.json"
```

Before the next expensive phase, persist a manifest binding the commit,
source/input/output hashes, exact commands, dependency versions, exit codes,
peak RSS/disk, logs and proof flags. Upload compact verified phase evidence
immediately. Never replace a failed condition with an edited status flag.
The kernel executable currently returns zero even for multiple components;
the explicit JSON condition above is necessary.

Read-only scalar/integration-rule replay:

```bash
python -B - "$D67_RUN" <<'PY'
import json, sys
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, "artifacts/track_b/d67/spectral")
from balanced_envelope import budget
from verify_subcell_mass_rule import verify
out=Path(sys.argv[1])
(out/"eta13-budget.json").write_text(json.dumps(budget(Fraction(1,13)),indent=2)+"\n")
(out/"subcell-rule.json").write_text(json.dumps({"checks":[verify(d) for d in range(4)],
    "full_finite_mass_enclosure_verified":False,
    "matrix_positivity_verified":False,"spectral_exclusion_certified":False},indent=2)+"\n")
PY
```

Some old CLIs, including `verify_subcell_mass_rule.py`, `certify_status.py`
and `negative_checks.py`, overwrite frozen reports in the checkout. Use their
read-only function interfaces or disposable worktrees; do not alter historical
reports just to obtain a fresh log. Historical `certify_status.py` exits 2 and
lists superseded gates; preserve it but do not make it the new final verifier.

Optional diagnostic after resource admission:

```bash
python -B -u "$S/adaptive_matrix_probe.py" "$D67_RUN/full" --balanced --mass-depth 3 > "$D67_RUN/assembly.log" 2>&1
python -B -u "$S/adaptive_matrix_probe.py" "$D67_RUN/full" --eigenprobe --eta 1/13 --iterations 100 > "$D67_RUN/eigenprobe.log" 2>&1
```

The assembler stores a legacy eta=1/10 z; eigenprobe correctly recomputes z
from the loaded M for its requested eta. Bind the actual invocation and the
matrix hashes. Assembly flags do not reassemble matrices in eigenprobe mode.
The example starts seven trial modes. A single-mode continuation requires a
bound (2042568,1) initial-vector .npy file and --initial-vectors PATH; there is
no single-mode CLI flag. These commands are diagnostics and can never set a
spectral certificate flag.

Preservation and existing-field replay:

```bash
python -B tests/mutation/regression_checks.py replay
python -B tests/d11_regressions.py replay
python -B tests/d19_regressions.py replay
python -B tests/d19_regressions.py protected
python -B -O -m groups.d2_inventory
python -B -O -m groups.d7_inventory
python -B -O -m groups.d11_inventory
python -B -O -m groups.d19_inventory
git diff --check
git status --short
```

Preserve known baseline failures: the Gaussian Track-B artifact ledger SHA
mismatch and stale d=11/d=19 incomplete-inventory tests. Report their source
and exact logs; never rewrite historical hashes or weaken tests to mask them.
Use each field's review guide for its full required replay, beyond this list.

## Transfer and final acceptance

d=43 and d=163 need independent geometry, face pairings/stabilizers, cusp
lattice and dual minima, truncation, mesh/P, scalar ledger and global index.
These fields are PID/class number one, not norm-Euclidean; do not assume
Gaussian nearest-point division. Transfer only audited analytic lemmas and
software interfaces. In particular, d=67's Y=2 and Robin coefficient 1/2 do
not transfer automatically to d=163: 8*pi/sqrt(163)<3, unlike d=67. Derive the
appropriate truncation-dependent threshold form and constants first.

Each completed field must supply:

- Frozen exact group identity, analytic dependencies and source/artifact hashes.
- Complete route-specific evidence: inventory/orbital normalization,
  admissible-test positivity/support, quadrature and tails for trace proofs;
  geometry, interpolation/coefficient/error and cusp/infinite-tail bounds
  for CR/index proofs. Existing trace certificates need no FEM ledger.
- A genuine exhaustive exclusion/index certificate with rigorous arithmetic.
- A clean independent reconstruction command, normal and optimized replay,
  and mutations for missing contributions, identities, bindings and positivity.
- An honest report of passed/failed conditions; a failed dependency prevents
  a theorem claim even if another computation is favorable.

The all-nine manifest may mark the project complete only when every field has
such evidence. Publish source, compact certificates and verified logs at each
completed phase. Save restartable long-run checkpoints with strong bindings;
regenerate large deterministic binaries rather than trusting transient storage.
PR #10 was already merged when this handoff was published. Keep PR #11
and subsequent spectral-proof follow-ups draft until their actual theorem
gates pass; merging a checkpoint does not certify a theorem. Publishing
checkpoints is authorized; merging unrelated PRs or paid compute is outside
this handoff. Do not broaden to congruence levels or other fields.
