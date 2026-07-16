# Track B observable collocation milestone

Date: 2026-07-15

## Result

The exact-row observability problem is repaired:

```text
M=4 cutoff_rank_certified = true
M=5 cutoff_rank_certified = true
all four cusp blocks present
480 full physical rows independently reconstructed
retained-mode Bessel fallbacks = 0
```

The coefficient-stability milestone fails:

```text
delta_2,3 <= 2.343595953986195
delta_3,4 <= 2.933840152773255   (required < 0.25)
delta_4,5 <= 11.303563492682994  (preferred < 0.1)
```

Consequently no `track-b-frozen-trial/v1` manifest was emitted and
`trial_frozen=false`, `global_hejhal_defect_certified=false`,
`rung4_certified=false`, and `dual_certification=false` remain correct.

## Why the old six-point family failed

The old M=4 normalized matrix had numerical rank 73 rather than 80.  A
256-bit Arb-midpoint Gram eigensystem gives a seven-dimensional diagnostic
nullspace.  Between 73.9% and 87.5% of every reported null basis vector lies
in modes first introduced at M=4.  Dominant shells include

```text
infinity: 4
zero: 17/5, 18/5, 20/5
```

The least-sensitive row families are consistently the parabolic `T1` rows
and the infinity-to-infinity block.  This agrees with the earlier observation
that the null vectors were concentrated in newly added high-shell modes.
The complete basis decomposition by cusp, shell, mode, real/imaginary part,
unit orbit, and row family is recorded in
`track_b_M4_nullspace_audit.json`.  The basis is diagnostic: nearly degenerate
null vectors can rotate within their seven-dimensional subspace.

The failure was therefore insufficient phase/height/fiber observability in
the small structured point family, not insufficient Arb precision.

## Exact replacement family

The replacement uses 20 hand-auditable rational points.  Every height exceeds
one, every point is in the closed Humbert core, and the coordinates avoid
`x1=0`, `x2=0`, and `x1=+/-x2`.  The denominators and height layers deliberately
do not form a symmetry-closed grid.

Each point expands to all four exact pairing relations and all six source
fibers, giving 480 physical rows.  No row is synthetic.

```text
exact phase aliases through M=5: 0
phase coherence diagnostic:      0.6047035377
distinct exact height signatures: 20
deduplicated source/target heights: 40
```

The radial profiles are still highly coherent: the worst M=4/M=5 diagnostic
is approximately 0.999868 between cusp-zero shells `17/5` and `18/5`.
Nevertheless, the combined phase, height, cusp, and fiber information gives
full candidate rank.  This distinction is why radial coherence alone was not
used as a rejection gate.

The canonical collocation-family hash is

```text
604a4118548db52277e9c37cbed151646e6c3de38035950e4456b0091c924446
```

## Adaptive enrichment and certified solves

The enrichment loop starts with the exact normalization row and one
independent row from each cusp block.  It repeatedly computes the
least-observed right-singular direction, scores unused exact physical rows by
their action on that direction, and appends the deterministically best row.
The final square subsystem is reconstructed in Arb and passed to the existing
power-of-two scaling plus Neumann/Rump contraction solver.  All unselected
physical rows remain in the independent verification set.

| M | unknowns | candidate rows | selected physical | omitted rows | certified sigma-min lower | contraction upper | full residual L2 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 44 | 480 | 43 | 437 | 1.0671e-2 | 8.92e-14 | 6.8221e-5 |
| 3 | 52 | 480 | 51 | 429 | 9.8815e-3 | 1.28e-13 | 4.5664e-5 |
| 4 | 80 | 480 | 79 | 401 | 1.1640e-3 | 2.01e-12 | 1.2280e-5 |
| 5 | 100 | 480 | 99 | 381 | 1.4706e-4 | 2.16e-11 | 1.7710e-6 |

The selected-row L2 residual is about `6.88e-16` at M=4 and `1.20e-14` at
M=5.  The omitted-row residuals are respectively `1.2280e-5` and `1.7710e-6`.
Thus the full residual is not being replaced by the square-system residual;
all omitted rows are present and independently reconstructed.

Condition estimates worsen from about `3.6e5` at M=2/M=3 to `5.2e6` at M=4
and `1.15e8` at M=5.  The certified rank remains positive, but coefficient
sensitivity increases rapidly.

## Where coefficient instability lives

The exact normalization `a_infinity,(1,0)=1` is identical at every cutoff and
fixes both global scale and global phase.  Restriction maps are exact
label-preserving sparse maps; no fitted phase rotation was used.

For M=3 to M=4:

```text
cusp infinity delta <= 0.6740
cusp zero     delta <= 2.8554
```

The instability is primarily a cusp-zero redistribution.  The newest shared
cusp-zero shell contributes about 1.3498, so the drift is not confined to
modes newly added at M=4.

For M=4 to M=5:

```text
cusp infinity delta <= 8.5666
cusp zero     delta <= 7.3746
```

The previously retained infinity norm-4 shell contributes about 8.3253 and
the retained cusp-zero norm-20 shell about 7.1519.  Adding the M=5 tail causes
order-one changes in old modes at both cusps.  Fiber-permutation inconsistency
is not applicable because the coefficients are shared cusp-channel
coefficients rather than six independent fiber coefficients.

## Independent least-squares crosscheck

To test whether adaptive square-row choice caused the drift, a separate
floating-point diagnostic solved the full 480-row constrained least-squares
problem at every cutoff.  It also gives nonconvergent low modes:

```text
M=3 -> M=4 delta = 1.7234
M=4 -> M=5 delta = 6.9164
```

The corresponding constrained condition estimate rises from about 156 at M=3
to 1,051 at M=4 and 27,709 at M=5.  This diagnostic is not a certificate, but
it falsifies the hypothesis that changing only the square-row selector will
restore convergence.

## Remaining physical bottleneck

The collocation family now observes the finite spaces, but the normalized
fixed-`r` coefficient family is not convergent.  The evidence does not support
freezing any current cutoff.  The next falsifiable diagnostic should track the
spectral parameter with cutoff—locate the minimum singular direction of the
same full oversampled operator as a function of `r` at M=3,4,5—before another
continuum calculation.  If the minimizing `r_M` does not stabilize, the
current truncated formulation is not producing a convergent trial family. If
it does stabilize, each selected `r_M` still requires a fresh interval solve
and exact cutoff comparison.

No floor, mass, partition, or continuum certificate was recomputed.  The old
floor and projected-mass artifacts remain valid for their original trial and
are explicitly recorded as incompatible with this current trial in
`track_b_legacy_trial_compatibility.json`.

## Outputs

- `track_b_collocation_family.py`
- `track_b_M4_nullspace_audit.py`
- `track_b_observable_hejhal.py`
- `track_b_observable_verify.py`
- `track_b_observable_stability.py`
- `track_b_oversampled_ls_diagnostic.py`
- `track_b_legacy_trial_compatibility.py`
- `track_b_collocation_family_result.json`
- `track_b_M4_nullspace_audit.json`
- `track_b_observable_M2_result.json` through `track_b_observable_M5_result.json`
- independent M=4 and M=5 verification JSON artifacts
- per-cutoff full-row and enrichment JSONL ledgers
- `track_b_observable_stability_result.json`
- `track_b_trial_freeze_status.json`
- `track_b_oversampled_ls_diagnostic.json`
- `track_b_legacy_trial_compatibility.json`
- `test_track_b_observable_collocation.py`
