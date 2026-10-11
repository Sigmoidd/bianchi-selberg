# Global cusp symmetry audit (Γ₀(2+i))

Branch `research/global-cusp-symmetry` off `main` @ `f233ef4ba7cde6e9811cd375e5135e10afaafdf3`
(merge of PR #2).  Companion: `docs/GLOBAL_CUSP_SYMMETRY_THEOREM.md`.

**Outcome.** No global eigenvalue enclosure is produced, and no certification flag is changed.
`global_hejhal_defect_certified`, `rung4_certified`, `hejhal_existence_certified` and every
`*_certified=false` channel in `track_b_global_hejhal_defect_result.json` are untouched (asserted
by `FailClosedTests`).  What is established is a proved global symmetry theorem and exact algebra;
what is measured is experimental.

## 1. Repository state recorded

* HEAD `f233ef4`; remote branches: `codex/d2-exact-inventory`, `codex/d7-certificate`,
  `codex/d11-certificate`, `codex/d19-certificate`, `codex/field-interface-and-rigor-safety`,
  `codex/group-backends-and-proof-docs`, `ladder-theorems-34`, `track-b/d43-geometry`,
  `track-b/d67-design`, `track-b/d67-spectral-exclusion`, `main`.
* No other branch changes `track_b*`, `TRACK_B*`, `six_copy*` or `theorem_DK*` relative to `main`
  (the d11/d19/d43/d67 branches are other fields; `ladder-theorems-34` is older, 2026-07-12).
  Nothing supersedes `TRACK_B_GLOBAL_HEJHAL_DEFECT_AUDIT.md` (2026-07-15) for Γ₀(2+i).
* Existing proof statuses preserved: `six_copy_DK_algebra_result.json` regenerates byte-identically;
  Track B ledger as in the 2026-07-15 audit (M=1 physical assembly certified; defect, continuum,
  floor/mass compatibility, tails, reprojection not certified).

## 2. Test outcomes (independent of the new work)

| Run | Result |
|---|---|
| Existing suite before any change (`tests/`, `test_track_b_*.py`, `test_verified_hejhal.py`) | 82 passed, 1 failed, 8 subtests |
| Same suite + `tests/test_global_cusp_symmetry.py` | 107 passed, 1 failed, 8 subtests (25 new tests) |
| `six_copy_DK_algebra.py` regression | output byte-identical to committed JSON |

The single failure, `test_track_b_two_cusp.py::test_closing_artifacts_verify_independently`,
**predates this work** and is not caused by mathematics: the verifier reports
`deterministic_hashes_match = false` only because the committed artifacts' recorded hashes were
produced from CRLF (Windows) bytes while this checkout is LF.  On CRLF copies of
`track_b_two_cusp_result.json` and `track_b_hejhal_rows.jsonl` (temporary files, not committed)
`deterministic_hashes_match = true` and `verified = true`.  Likewise the recorded
`trial_sha256` of the projected-mass certificate (`c77b711b…`) equals the SHA-256 of the CRLF form
of `six_copy_hejhal_balanced_coeffs.json` (LF form: `3f572f36…`).  Existing artifacts were not
modified; a `.gitattributes` or LF-normalized hashing would remove this platform dependence.

## 3. Evidence by tier

**PROVED** (exact; replay with `python3 global_cusp_character_algebra.py`, `global_quarter_turn.py`,
`global_conjugate_prime.py`, `global_oldspace.py`): everything listed in Theorems A1–A3, B1–B3,
C1–C2, D1 (algebraic part), E (algebraic part) of the theorem document.  Analytic steps B4, E2,
F1 are proved there from named standard theorems; the one citation (discreteness of the cuspidal
spectrum, EGM Ch. 6) is given at chapter level and its theorem number has not been pinned.

**CONDITIONAL**: F3 (a glued trial W is exactly odd only if the partition of unity and gluing are
Q-equivariant — unchecked); E Gram statements (for Hecke eigenfunctions); closedness of
Ran B₁+Ran B_π (needs a uniform bound toward Ramanujan, cited not verified); any newspace
statement (needs certified level-one odd spectral data, which the repository lacks).

**EXPERIMENTAL** (floating SVD / least squares on Arb midpoints; `global_cusp_symmetry_experiment.py`,
`global_cusp_symmetry_scan.py`, `global_trial_compatibility.py`; none feeds a certificate):

*Observability (144 rows = the repository's 6 collocation points, r = 6.7439020359331625).*

| M | unknowns full | rank full | unknowns per parity | rank per parity | min σ (col-equilibrated) full / even / odd |
|--:|--:|--:|--:|--:|--|
| 1 | 24 | 24 | 6 | 6 / 6 | 6.1e-2 / 4.1e-1 / 3.5e-1 |
| 2 | 44 | 44 | 11 | 11 / 11 | 1.3e-2 / 8.8e-2 / 5.4e-2 |
| 3 | 52 | 52 | 13 | 13 / 13 | 7.8e-3 / 5.7e-2 / 5.1e-2 |
| 4 | 80 | **73** | 20 | 20 / 20 | 1.2e-16 / 1.0e-2 / 5.5e-3 |
| 5 | 100 | **73** | 25 | 25 / 25 | 1.9e-17 / 9.8e-3 / 3.9e-3 |
| 6 | 116 | **73** | 29 | 29 / 29 | 7.4e-19 / 6.8e-3 / 2.3e-3 |

The documented M=4 deficiency is reproduced and has a plain explanation: the six points supply
a numerical rank of 73 at each of M = 4, 5, 6 (80 − 73 = 7, matching the seven near-null directions of
the audit), so a system with more unknowns than that is underdetermined.  Column equilibration does
not change the rank, so it is not a scaling artifact.  The unprojected least-squares "residual"
at M ≥ 4 (~1e-19) is an exact fit of an underdetermined system and carries no information.
Parity sectors (≤ 29 unknowns) stay overdetermined and well conditioned.  **But the same fix is
available without symmetry**: with 24 or 60 random exact points the unprojected system has full
rank 116/116 at M = 6 (scan file).  Symmetry is a 4× cheaper alternative, not the only cure.

*Cutoff stability (normalized least squares, c_(1,0) = 1; L² difference of solutions on common
modes; smaller-cutoff norm in brackets).*

| system | 1→2 | 2→3 | 3→4 | 4→5 | 5→6 |
|--|--|--|--|--|--|
| odd (6 pts) | 0.63 (2.7) | 0.87 (2.7) | 0.78 (2.9) | 7.6 (3.8) | 9.5 (21.7) |
| even (6 pts) | 2.2 (2.4) | 1.5 (6.1) | 3.8 (4.8) | 18.0 (16.5) | 36.9 (40.5) |
| unprojected (6 pts, unreliable at M≥4) | 0.57 | 1.8 | 2.2 | 4.8 | 0.17 |
| odd, 60 random pts | 0.24 | 0.54 | 1.01 | 1.56 | 0.23 |
| unprojected, 60 random pts | 0.18 | 0.46 | 1.07 | 1.78 | 0.27 |

Neither parity restriction nor more collocation points stabilizes the low coefficients: drift is
O(1) for M = 3→5 in every system.  The original audit's finding (no numerical evidence of cutoff
convergence) is **not** repaired by this branch.

*Parity of the trials.* With 60 random points the odd and unprojected least-squares residuals
coincide to within ~20% (e.g. M=1: 1.2e-4 vs 1.0e-4) while the even-sector residual is about 3× larger
and nearly flat in M (3.3e-4 → 2.1e-4 over M = 1…6): the approximate solution at r = 6.7439 is essentially Q-odd.  The recorded
mass/floor trial (r = 6.62208, M_∞ = 8, M_0 = 40) has coefficient-norm fractions odd = 0.9999997,
even = 7e-4 (R-even fraction 1.0).  These are coefficient-norm proxies, **not** L² mass fractions.

*Work accounting.* Unknowns shrink exactly 4× (24→6, 80→20, 100→25, 116→29); the dense interval
solve shrinks by roughly 4³.  Distinct Arb K_{ir} evaluations do **not** shrink (204 for 80
unknowns at M=4, merged already by the cache because rotation-equivalent modes share |μ|).
Whether the pairing-face and continuum overlap calculations halve under Q depends on a
Q-equivariant partition; that has not been checked, so no reduction is claimed there.

*Field-level check.* On the repository's own `component_row`, F_{q(c)}(QP) = ε F_c(P) holds in
Arb (every difference enclosure contains 0) for both parities, and the same test rejects a
wrong sign and the identity in place of q (`field_level_intertwining_check`).

## 4. Answers to the six questions

1. **Does the fifth-root conjugation give a global reduction beyond the cusp-local decomposition?**
   Not from fiber algebra: the commutant of the full image is span{I,J}, [A_T, ρ(S)] ≠ 0, and no
   nonzero equivariant constant fiber operator kills both singular directions (Lemma A2,
   Corollary A3).  The map k ↦ 3k on characters comes from q, which is only half of the global
   symmetry.  A global reduction exists through U_Q = (q, base rotation): factor 4 on
   Fourier–Whittaker unknowns (factor 2 beyond the point reflection R).
2. **Is the odd projector rigorously cuspidal on the entire Γ₀(2+i) quotient?** Yes (Theorem B4),
   by two independent arguments (vanishing constant terms at both cusps; Eisenstein/residual
   spectrum is Q-even), with discreteness of the cuspidal spectrum cited from the standard theory
   (theorem number unpinned).  Non-emptiness of the odd sector is not claimed.
3. **Does it reduce certification work or improve certified residual bounds?** Work: yes in
   linear-algebra dimension (4×), no in distinct Bessel evaluations, unknown for face/continuum
   integrals.  Bounds: exactly R₋ ≤ R (Corollary F2) but no certified improvement exists, and the
   projected mass μ₋ has no certificate.  Cutoff drift is unchanged (§3).
4. **Is the conjugate-prime correspondence independently verified?** The algebra is replayed
   exactly (levels, action on ℍ³, cusp lattices, fiber transport, C Q C⁻¹ = Q⁻¹) and Theorem D1
   proves the unitary intertwiner and isospectrality.  I did **not** run a second numerical
   Hejhal/spectral computation for π̄ as a cross-check.
5. **Can genuine newspace be separated from both degeneracy lifts?** Algebraically yes: both
   lifts and the Fricke element are verified, P_old,full exists, commutes with U_Q and Δ, and
   P_new,₋ is an orthogonal projector with cuspidal range (Theorem E2).  Computationally no:
   applying P_old,full needs certified level-one odd spectral data, so newness stays uncertified.
6. **First remaining obstruction to a global eigenvalue certificate?** In pipeline order
   (the 2026-07-15 flow chart) it is **cutoff stability of the normalized physical solve**: no
   cutoff at which the low coefficients stop moving has been found, with or without the parity
   restriction, and no Arb-verified solve exists at M ≥ 4 for the odd-sector system.  Independent of
   that, the trial identity split remains a hard blocker: the physical solve (r = 6.7439020359331625,
   M = 1, coefficient hash `549511ef…`) and the floor/mass certificates (r = 6.62208,
   M_∞ = 8, M_0 = 40, trial `c77b711b…`) are different fields; the floor artifact carries no trial
   hash or spectral parameter (`python3 global_trial_compatibility.py` → `compatible: false`).
   Behind those: continuum d₀/d₁ ledgers, non-floor faces, tails, reprojection, and the projected
   mass μ₋ for the same W.

## 5. Suggested next steps (not done here)

1. Pick one canonical trial and regenerate floor and mass for it (or the physical solve at the
   existing r), with `trial_sha256`, r, partition and theorem hashes shared; hash over LF-normalized bytes.
2. Add independent non-symmetric collocation points (24+ already give full rank) and test whether
   any cutoff gives a stable odd-sector solve; only then run an Arb-verified solve on the 20–29
   odd-sector unknowns.
3. Check whether the Track B partition and gluing maps are Q-equivariant; if so, integrate over a
   Q-fundamental domain and certify μ₋² = ½(‖W‖² − Re⟨W, U_Q W⟩) for the same W.
4. Obtain certified Q-odd level-one spectral data before saying anything about newness.

## 6. Files added

`global_cusp_gaussian.py`, `global_cusp_character_algebra.py`, `global_quarter_turn.py`,
`global_conjugate_prime.py`, `global_oldspace.py`, `global_trial_compatibility.py`,
`global_cusp_symmetry_experiment.py`, `global_cusp_symmetry_scan.py`,
`tests/test_global_cusp_symmetry.py`, the two docs, and the `global_*_result.json` outputs of the
runs above.  No existing file was modified.
