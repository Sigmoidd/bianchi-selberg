# Ford-vertex Astra emergency recovery record (2026-10-11)

**Status: archive / recovery, NOT a recovered implementation or spectral certificate.**

- Repository: \`Sigmoidd/bianchi-selberg\`.
- Intended missing feature branch: \`research/ford-vertex-equivariance\`.
- Published prerequisite branch: \`research/global-cusp-symmetry\` at commit \`b337589612006e9e8752ab074e73914899ecd645\` when checked on 2026-10-11.
- Earlier user-uploaded prerequisite bundle: \`global-cusp-symmetry.bundle\` with head \`b5f29bd1b314e4d385e5af01805729ba27e3f3e2\` and prerequisite ancestor \`f233ef4ba7cde6e9811cd375e5135e10afaafdf3\`.
- Preservation purpose: the user exhausted Astra tokens and wants to avoid losing unpublished research.
- Search result: no pushed \`research/ford-vertex-equivariance\` branch or matching PR, no Ford-vertex source bundle in currently accessible uploaded files. The original Astra ephemeral execution environment is not accessible from this conversation.

## A. Mathematical foundation

For \`H=Gamma_0(2+i) <= PSL_2(Z[i])\`, define \`Q(z,y)=(iz,y)\`. Q normalizes H and Q^2 belongs to H; the induced \`U_Q\` is a global unitary involution commuting with the Laplacian. Both cusp constant terms are Q-even, yielding a global cuspidal projector \`P_-=(I-U_Q)/2\`. This avoids mistaking a continuous-spectrum wave packet for an odd cuspidal quasimode. The prerequisite global-cusp-symmetry work established exact fifth-root character algebra, bundle-level six-copy Q action, and Gaussian conjugate-prime correspondence.

The corresponding Fourier symmetry reduces the number of coefficient unknowns fourfold; local fifth-root averaging alone is *not* a global projector.

## B. Ford geometry to preserve

Published source \`track_b_partition_geometry.py\` contains the exact Picard reference cell and special points. The cell is \`x1 in [-1/2,1/2]\`, \`x2 in [0,1/2]\`, \`rho=x1^2+x2^2+y^2 >= 1\`.

| Point | z | y² | Full Picard stabilizer order | Six-fiber stabilizer orbit lengths |
|---|---|---:|---:|---|
| v_00 | 0 | 1 | 4 | 2+2+2 |
| v_0h | i/2 | 3/4 | 6 | 3+3 |
| v_mh0 | -1/2 | 3/4 | 6 | 3+3 |
| v_ph0 | 1/2 | 3/4 | 6 | 3+3 |
| v_mhh | (-1+i)/2 | 1/2 | 12 | 6 |
| v_phh | (1+i)/2 | 1/2 | 12 | 6 |

The Q-odd scalar function has exact zero value and zero first derivatives along the vertical lines over \`z=0\` and \`z=(±1+i)/2\`: they satisfy \`iz=z+lambda\` with \`lambda in Z[i]\`. The two corner lines are quotient-equivalent. On the six-copy induced field, use **full fiber transport** before interpreting these as zero-value equations; it is false that every fiber component must independently vanish at the same reference-cell coordinate.

The existing \`TRACK_B_GLOBAL_PARTITION_CERTIFICATE.md\` reports an independently verified, 192-bit Arb global normalized partition and complete singular-stratum incidence; it does **not** certify a global eigenvalue.

## C. Astra Ford-vertex results reported in prior research conversation

**These numbers are carried over from the previous session's final report; underlying unpublished implementation/run artifacts have not been recovered or independently replayed in this recovery branch.**

- Exact Ford vertex/stabilizer six-fiber actions and Q-odd value/first-derivative projectors were reportedly computed and verified.
- A symmetry-adapted collocation scan ran to Fourier cutoff **M=8**.
- Primary tracked spectral parameter near **r = 6.62205583**, close to the known Picard level-one mode.
- M=7→8 coefficient drift: **0.00268** baseline vs **0.00179** with vertex-jet constraints (roughly **33% reduction**).
- **Held-out physical residual increased for the constrained trial.** Smaller coefficient drift is not a proof of better approximation, convergence, or eigenvalue existence.
- Secondary target near \`r=6.7439\` lacked a convincing tracked singular-value minimum.
- Reported finite constraint ranks and permitted dimensions: \`M2: 9 / 1\`; \`M3: 11 / 1\`; \`M4: 14 / 6\`; \`M5: 14 / 11\`; \`M6: 14 / 15\`. Independently replay these numbers before claiming certified status.
- No compatible certified same-trial global residual, both-cusp overlap/gradient bounds, projected-mass witness, spectral enclosure, or rung-4 closure was obtained.

## D. What is *not* included

This archive does **not** contain the Ford-vertex source tree, the actual collocation matrices, M2–M8 coefficient vectors, raw held-out residuals, full scan tables, original test logs, original Ford-vertex JSON results, or an Astra commit SHA. Those were not present in GitHub or the accessible user uploads. This note is not a substitute for the original source bytes.

Expected outputs named in the original handoff, **not confirmed present**: \`ford_vertex_equivariance.py\`, \`ford_vertex_equivariance.json\`, \`ford_equivariant_collocation.py\`, \`ford_spectral_stability.json\`, \`docs/FORD_VERTEX_GLOBAL_EQUIVARIANCE.md\`, \`docs/FORD_VERTEX_SPECTRAL_AUDIT.md\`.

## E. Rescue and resumption checklist

If the suspended original Astra workspace can still be opened: **save it before doing more computation**. Immediately run \`git status --short\`, \`git log -5 --oneline --decorate\`, save all unstaged/untracked files, commit research files on the actual feature branch, push it or export a \`git bundle\` with complete prerequisites, and save scripts plus raw M2–M8 scan JSON. Preserve original exact inputs and hashes.

If the original workspace cannot be recovered, begin a **new implementation**, not a falsely labeled retrieval:

1. Restore the published global-cusp-symmetry branch and verify the exact Q induced action.
2. Recompute each exact Ford vertex Q-image/reduction word and induced six-fiber stabilizer; verify allowed value/jet spaces with rational/cyclotomic arithmetic.
3. Determine if vertex jets supply independent constraints beyond the already Q-odd Fourier basis.
4. Compare baseline versus constrained scans **at exactly identical** mode sets, spectral parameters, normalization, sampling families and held-out rows for M2–M8. Track \`r_M\`, coefficient drift, and held-out residuals; capture actual numeric data in file artifacts.
5. Promote only a stable candidate with improved generalization to full Arb same-trial partition/overlap/gradient/tail/projected-mass/residual certification.

Key repository references: \`six_copy_DK_algebra.py\`, \`theorem_DK_sixcopy.tex\`, \`track_b_partition_geometry.py\`, \`TRACK_B_GLOBAL_PARTITION_CERTIFICATE.md\`, \`TRACK_B_OBSERVABLE_COLLOCATION_REPORT.md\`, \`track_b_observable_hejhal.py\`, \`track_b_overlap_arb.py\`, \`TRACK_B_DIRECT_WEIGHTED_REPORT.md\`, \`track_b_global_hejhal_defect.py\`.

## F. Strict research-status gates

Keep \`global_hejhal_defect_certified=false\` and \`rung4_certified=false\` until verified analytic closure. The reported constrained trial's higher held-out residual is an unresolved concern. Collocation diagnostics do not bound \`||(Delta-lambda)W||_2\`. Never mix certificates belonging to different coefficient vectors, cutoff parameters, or partition hashes. Do not infer newspace status by deleting only the invariant fiber line: prime-level oldspace has two degeneracy lifts.

**This branch is an emergency archive. Do not merge it into the main theorem branch as a completed Ford-vertex PR.**