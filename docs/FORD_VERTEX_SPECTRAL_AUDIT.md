# Ford-vertex spectral audit

Prepared 2026-10-11. The exact Ford geometry is proved and independently reconstructed. A bounded numerical odd candidate near r=6.62205583 survives M=8. Hard vertex jets reduce dimensions and condition estimates, but do not produce a smaller held-out residual. No global eigenvalue certificate or newspace claim is made.

## Evidence tiers

- **PROVED:** exact vertex images/reduction words, stabilizers and six-fiber actions, value/first-jet projectors, scalar nodal geodesics, full corner-vertex jet vanishing, and 14 genuinely added geometric equations. The doubled-partition derivative transfer is a theorem conditional on a valid input global atlas/section.
- **CERTIFIED:** independently reconstructed 192-bit Arb nonzero minors establish finite vertex-constraint ranks at both recorded fixed parameters. The original 2,048-leaf partition is reverified unchanged. A new conditional infinity-cusp witness has positive mass lower bound, independently rechecked on a finer 256-bit grid.
- **EXPERIMENTAL:** all finite physical matrix ranks, conditioning, least squares, singular vectors, parameter scans and coefficient drift. The experimental candidate has exact rational frozen coefficients, but has no verified global gluing certificate.
- **BLOCKED:** complete original/transformed/mixed continuum overlap coverage, global Dom(Delta) compatibility, same-trial non-floor and cusp/core residual channels, full projected mass, tail/channel compatibility and final residual-to-mass enclosure. Existing existence and Rung-4 flags remain false.

## Source and preserved foundations

The attached bundle was recovered at b5f29bd1b314e4d385e5af01805729ba27e3f3e2, whose required parent f233ef4ba7cde6e9811cd375e5135e10afaafdf3 was present. GitHub had no global-symmetry branch; a tree-identical mirror was published at b337589612006e9e8752ab074e73914899ecd645. Both source trees are a2469d0616caebc61d35c6a3615185069fff6db4. The PR stacks on that mirror. The bundle source remains the provenance reference in the exact ledger.

The existing global defect, partition, partition leaves, two-cusp result and six-copy algebra JSON remain byte-identical to the bundle source. `global_hejhal_defect_certified`, `rung4_certified` and `hejhal_existence_certified` remain false. No old floor or mass artifact is aggregated into this trial.

## Comparable mathematical inputs

All three families use the repository mode sets N(nu)<=M and N(beta)<=5M, identical physical rows, exact coset permutations, the same r, and normalization a_infinity,(1,0)=1. The baseline is the full six-copy coefficient formulation; the odd control uses its established C4 reduction; the adapted family additionally satisfies the independent vertex first-jet equations.

Training uses 32 distinct deterministic rational points in the true curved half-cell, with seed 173; the independent hold-out uses 16 disjoint points, seed 947. There are 768 training and 384 held-out generator/fiber rows, including all four cusp blocks. Points and hashes are serialized. Every family retains equal physical row weights. RMS is ||A c||/sqrt(number of physical rows), not hyperbolic L2 norm.

The 144 raw vertex-jet equations are first reduced in exact arithmetic modulo parabolic, Q and translation identities to 14 equations. Zero first-jet columns are proved in Q(zeta10) and retained as kernel directions before equilibration. Without these two steps, scaling roundoff-sized implied rows/columns fabricated ranks 18 and 10/12; those preliminary results were rejected. No such results enter the committed evidence.

## Exact constraint ranks and finite spaces

| M | Full unknowns | Odd unknowns | Certified added rank at both fixed r | Permitted odd+jets unknowns | Infinity normalization available |
|---:|---:|---:|---:|---:|---|
| 2 | 44 | 11 | 9 | 2 | no |
| 3 | 52 | 13 | 11 | 2 | no |
| 4 | 80 | 20 | 14 | 6 | yes |
| 5 | 100 | 25 | 14 | 11 | yes |
| 6 | 116 | 29 | 14 | 15 | yes |
| 7 | 128 | 32 | 14 | 18 | yes |
| 8 | 152 | 38 | 14 | 24 | yes |

At M2/M3 the remaining two nonzero fields have zero (1,0) coefficient, so the prescribed normalization is unavailable; the finite space itself is not zero. At the fixed parameters, Arb minor lower bounds equal exact algebraic rank upper bounds. No interval rank claim is made over entire scan regions, and the numerical nullspace basis is not an interval-certified coefficient solve.

The local allowed odd dimensions (value, full first gradient) are (1,1) at v_00, (2,2) at each of the three linked order-six representatives, and (0,0) at both corners. Their simultaneous component jet dimensions are 2,4,0. Thus the three order-six representatives supply one coupled four-dimensional space. The two corner nodal lines are translation-equivalent. Exact details are in FORD_VERTEX_GLOBAL_EQUIVARIANCE.md.

## Fixed-parameter held-out comparison

| r | M | Baseline RMS | Odd RMS | Odd+vertex jets RMS |
|---|---:|---:|---:|---:|
| 6.62208 | 2 | 1.42738e-06 | 1.34409e-06 | — |
| 6.62208 | 3 | 1.36361e-06 | 1.32421e-06 | — |
| 6.62208 | 4 | 2.63845e-07 | 1.34646e-07 | 6.04634e-06 |
| 6.62208 | 5 | 4.26895e-09 | 1.25095e-09 | 1.9674e-09 |
| 6.62208 | 6 | 4.26905e-09 | 1.25208e-09 | 1.96774e-09 |
| 6.62208 | 7 | 4.26952e-09 | 1.25213e-09 | 1.96788e-09 |
| 6.62208 | 8 | 5.31045e-09 | 1.54148e-09 | 2.03003e-09 |
| 6.7439020359331625 | 2 | 1.88574e-06 | 2.2806e-06 | — |
| 6.7439020359331625 | 3 | 1.79057e-06 | 2.15234e-06 | — |
| 6.7439020359331625 | 4 | 5.72199e-07 | 7.16245e-07 | 8.05877e-05 |
| 6.7439020359331625 | 5 | 6.36245e-07 | 7.18842e-07 | 8.33469e-07 |
| 6.7439020359331625 | 6 | 6.33216e-07 | 7.17602e-07 | 8.30921e-07 |
| 6.7439020359331625 | 7 | 6.37828e-07 | 7.18081e-07 | 8.2537e-07 |
| 6.7439020359331625 | 8 | 8.25103e-07 | 7.52087e-07 | 8.24092e-07 |

At r=6.62208, M8, hard jets raise the held-out RMS by 1.31693 times. They also worsen M4->M5 drift substantially. This comparison does not support a residual improvement from hard vertex jets.

| M, r=6.62208 | Baseline physical condition | Odd physical condition | Odd+jets physical condition |
|---:|---:|---:|---:|
| 4 | 813.844 | 272.781 | 23.4685 |
| 5 | 106021 | 24955.9 | 8676.52 |
| 6 | 106023 | 24956.3 | 8677.39 |
| 8 | 230175 | 25779.2 | 8833.45 |

Condition estimates are midpoint diagnostics in the serialized bases; column-equilibrated estimates are also stored. Better conditioning or fewer columns is not evidence of a better global residual.

## Spectral parameter scan with a fixed coefficient metric

The final objective is the smallest singular value of the physical operator on an orthonormal full-coefficient basis. It is invariant under unitary rotation of the constrained kernel. No column equilibration is used to choose r_M. Equilibrated matrices remain conditioning diagnostics only. An independent review rotated the same M5 kernel 12 ways: the discarded equilibrated objective changed by a factor 3.51, while its unscaled singular value stayed fixed. The final scans were rerun after correcting this.

Each original region has 17 deterministic grid samples; bounded local minimization refines each detected interior minimum and endpoints are retained. Region bounds are [6.55,6.69] and [6.69,6.81]. An endpoint optimum is not a tracked interior spectral candidate. The M7/M8 extension scans only the promising Picard region, retaining both fixed-r comparisons.

| M | Baseline r_M | Odd r_M | Odd+jets r_M |
|---:|---:|---:|---:|
| 2 | 6.55000000000 (edge) | 6.55000000000 (edge) | — |
| 3 | 6.55000000000 (edge) | 6.55000000000 (edge) | — |
| 4 | 6.62368958572 | 6.62957142752 | 6.62177989222 |
| 5 | 6.62256629221 | 6.62204908061 | 6.62183276359 |
| 6 | 6.62256629186 | 6.62204908017 | 6.62183276355 |
| 7 | 6.62256629129 | 6.62204908014 | 6.62183276344 |
| 8 | 6.62225277068 | 6.62205582791 | 6.62183251068 |

| M | Prior-region baseline r_M | Prior-region odd r_M | Prior-region odd+jets r_M |
|---:|---:|---:|---:|
| 2 | 6.690000000 (edge) | 6.690000000 (edge) | — |
| 3 | 6.690000000 (edge) | 6.690000000 (edge) | — |
| 4 | 6.690000000 (edge) | 6.690000000 (edge) | 6.810000000 (edge) |
| 5 | 6.810000000 (edge) | 6.690000000 (edge) | 6.690000000 (edge) |
| 6 | 6.810000000 (edge) | 6.690000000 (edge) | 6.690000000 (edge) |

The former r~6.7439 region has no robust interior minimum through M6. Its unprojected high-cutoff singular directions have very large normalized coefficient norms and worse held-out errors, consistent with weakly observable tail directions rather than a stable spectral family. No spectral exclusion is inferred from the failed scan.

## Coefficient convergence without fitted phase or rescaling

All drift uses the exact cusp/mode-labelled restriction embedding and the fixed (1,0)=1 normalization. Newly added mode norms, common-mode drift, cusp contributions and symmetry defects are retained in the JSON. No optimal phase rotation or norm adjustment is applied.

| Step | Fixed r=6.62208 odd delta | Fixed r=6.62208 jets delta | Scanned Picard odd delta | Scanned Picard jets delta |
|---|---:|---:|---:|---:|
| 2->3 | 0.390437 | — | 0.0400944 | — |
| 3->4 | 0.600243 | — | 0.0599352 | — |
| 4->5 | 1.02029 | 14.5175 | 0.892919 | 20.6136 |
| 5->6 | 0.00162095 | 0.00111173 | 0.000722886 | 0.000373197 |
| 6->7 | 4.72622e-05 | 0.000395786 | 4.07705e-05 | 0.000131125 |
| 7->8 | 0.00374917 | 0.00181909 | 0.00267789 | 0.00179389 |

M5,M6,M7 have identical infinity mode sets; M8 adds the norm8 shell. The bounded numerical gate requires two final common-mode deltas<0.1, two r changes<1e-4, interior minima, full numerical rank and held-out RMS<1e-6. Both odd families pass through M8; the unprojected family fails the final r-change threshold. This is numerical stabilization under a finite test, not certified cutoff convergence or evidence of a new level eigenvalue.

At M8 the scanned odd held-out RMS is 1.43777e-09; the adapted RMS is 1.53841e-09. Hard jets slightly improve the final parameter/drift diagnostic, but worsen held-out physical residual. They are useful exact necessary conditions, not an observed global spectral stabilization mechanism.

## Conditional analytic progress and remaining global obligations

ford_experimental_candidate.json freezes the M8 odd candidate at r=6.6220558279118809, complete mode order and exact rational coefficients, normalization, code/partition/atlas hashes and experiment identities. It is expressly not a track-b-frozen-trial/v1 or a certified global field. Existing verified-solve and global freeze gates are not promoted.

For any global odd extension agreeing with this infinity Fourier field on 1.3<=y<=1.4, Parseval on the torus modulo R gives

```math
||W_-||_2^2 >= 2 integral_{1.3}^{1.4} K_{ir}(2 pi y)^2 dy/y.
```

The coefficient orbit (±1,0),(0,±1) has exact moduli 1. The cusp slice embeds because a nonparabolic H element has |c|²>=5 and yprime<=1/(5y). Direct 192-bit Arb boxes on 128 height segments certify mu>2.1836190630e-6 under this explicit extension hypothesis. An independent 256-bit/512-segment replay proves a stronger lower bound>2.5304156803e-6 and validates the conservative witness. The projected-mass flag remains false because the global extension itself has not been constructed and certified. Once that matching global extension is proved, this cusp lower bound suffices for the spectral theorem denominator; an integral of its full mass is not additionally required.

Global odd projection does not require pointwise Q-equivariant weights: double the 33 labels to 66 half-weights. Isometry covariance preserves the original B0 and B1 conservative bounds. The finite incidence ledger includes all vertex, edge and face group words and Q conjugates, but it is not a continuum overlap cover. Original/Q mixed overlaps, boundary slivers, gluing/Domain compatibility and complete same-trial value/gradient defects remain open and fail closed.

A physical row RMS cannot be substituted for R_-=||(Delta-lambda)W_-||_2. No total residual, enclosure radius or eigenvalue interval is formed. The old floor/mass artifacts concern a different parameter/coefficient trial and remain incompatible.

The smallest next mathematical task is to turn the canonical M8 odd field into a complete globally glued Dom(Delta) field using the exact stabilizer averages and a fully accounted overlap cover, then certify its signed S-transition commutator with cancellation-preserving Taylor L2 integration. Start with the canonical floor-collar witness before expanding the same-trial ledger to all face/cusp channels. The conditional cusp witness gives a real mass budget once that extension is valid; loose global supremum bounds will squander it. Cutoff diagnostics help select a trial, while the spectral theorem ultimately requires this one field, its mass and its global residual.

## Hash portability, replay and independent results

New two-cusp results use schema v2 with sha256-utf8-lf/v1. Historical v1 digests are never rewritten: one audited CRLF ledger digest is explicitly mapped to its canonical LF digest under track-b-two-cusp-crlf-compatibility/v1. Unknown policies, unregistered mismatches, content changes and bare CR fail closed. The independent verifier reports the compatibility evidence rather than silently changing historical identity.

| Validation | Result |
|---|---|
| Existing baseline before edits | 107 passed, 1 historical CRLF hash failure, 8 subtests |
| Final regression and new tests | 130 passed, 14 subtests |
| Independent SymPy geometry reconstruction | 156 checks passed; no producer imports |
| Original-row numerical trial replay | 115 checks passed; parameter/coefficient/source/gate corruptions rejected |
| Independent bounded-extension replay | 55 checks passed |
| Arb vertex-constraint minors | exact rank lower=upper at both fixed r, M2 through M8 |
| Original partition reverification | 2,048 leaves pass; bounds and source bytes unchanged |
| Conditional cusp witness replay | 256 bits, 512 independent leaves; positive lower bound |

Install the tested stack with `python -m pip install -r requirements-ford-research.txt`, then reproduce from this branch:

```bash
python ford_vertex_equivariance.py
python ford_vertex_verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python ford_equivariant_collocation.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python ford_spectral_extension.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python ford_constraint_rank_arb.py
python ford_stability_verify.py
python ford_stability_verify.py --input ford_spectral_extension.json --output ford_extension_independent_verification.json
python ford_candidate_audit.py
python ford_cusp_mass_verify.py
python -m pytest -q tests test_track_b_*.py test_verified_hejhal.py
```

The floating scans are deterministic within a fixed numerical stack, while small BLAS/version differences may change last digits. Exact algebra, immutable historical hashes, interval determinant witnesses, conditional mass bounds and fail-closed gates are the replayable proof-facing outputs. No PR merge is authorized or performed.
