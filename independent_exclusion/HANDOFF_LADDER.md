# HANDOFF: certify Γ₀(𝔭) for N𝔭 = 9 and N𝔭 = 13 (CONGRUENCE.md §8)

Self-contained brief for a fresh agent. Read this, then CONGRUENCE.md
(especially §§1–2, 7, 8), then skim m3p_certify.py. Everything else is
reference.

## 1. What exists (all in this directory, all verified runs reproducible)

Project: rigorously exclude Laplace eigenvalues in (0,1) on Bianchi
quotients Γ\ℍ³ by verified finite elements — NO trace formula. Method
chain (per quotient): Lax–Phillips cusp reduction to a positivity
criterion on a compact core (DESIGN.md for level 1, CONGRUENCE.md §2 for
congruence level) → guaranteed Crouzeix–Raviart lower bound reducing it
to finitely many checks (Theorem G1 in lower_bound_theory.md; Theorem
G1𝔭 in CONGRUENCE.md §7) → verified execution (arb ball arithmetic for
all constants/enclosures; Rump BIT 46 (2006) Thm 2.3 + Cor 2.4 +
Lemma 2.5 for matrix positive-definiteness, implemented literally).

Certified theorems so far:
- λ₁ ≥ 1 for PSL(2,ℤ[i]) (level 1): m3_certify.py, master doc PROOF.md.
- **No eigenvalues in (0,1) for Γ₀((2+i)), N𝔭 = 5** (index 6, 2 cusps):
  m3p_certify.py; frozen certificate in PROOF.md "Theorem 2" and
  CONGRUENCE.md §6. Cusp-0 subtleties audited in 0cuspchecks.md.

**YOUR TASK: the next two rungs — 𝔭 = (3) (inert, N𝔭 = 9, index 10) and
𝔭 = (3+2i) (split, N𝔭 = 13, index 14).** Expected outcome: two new
certified theorems "Γ₀(𝔭)\ℍ³ has no eigenvalue in (0,1)". Plan of
record: CONGRUENCE.md §8 (this file operationalizes it).

## 2. File map

- `DESIGN.md` — level-1 criterion with proofs (Lemmas 0–3, exclusion
  theorem). The mathematical root; criterion generalizes per cusp.
- `lower_bound_theory.md` — Theorem G1 (level 1), Lemma I1
  (self-contained CR constant κ_sc = √(1/π²+1/15), arbitrary tets),
  Lemma R, Rump appendix (§5). NOTE the remark: τ_h = 0 (see §5 below).
- `CONGRUENCE.md` — congruence architecture: Reference-Cell Principle,
  Facts A–B (copy tiling, collar assembly), multi-cusp criterion (§2),
  gluing combinatorics (§3), 24-split mesh (§4), **Theorem G1𝔭 (§7)**,
  **ladder plan (§8)**.
- `0cuspchecks.md` — audit C1–C9 of the cusp-0 chain. C2 (fold-tiling)
  and C6 (precise invariance) are the proofs you'll reuse verbatim for
  the new primes.
- `congruence_prototype.py` — float pipeline for Γ₀(2+i): reference
  24-split mesh builder (`build_reference`), 𝔽₅ combinatorics
  (`build_gluing`), 6-copy assembly, dense + sparse-LOBPCG margin
  sweeps. Margin measured: μ(λ) ≥ 4.39 on (0,1).
- `m3p_certify.py` — the level-𝔭 certificate: arb reference-cell data
  (element enclosures via `weighted_moments`, constants), exact scatter,
  float parameter search (scalar conditions + constant-direction
  heuristic + LOBPCG pencil filter), Rump PSD per window. THIS is the
  file you generalize.
- `m3_certify.py` — level-1 certificate; imports reused by m3p
  (`tet_arb_data`, `mid_rad`, `upper`, `amax`, `amin`, `a_yf`, EPS, ETA).
- `PROOF.md` — master document (Theorems 1, 2, dependency list). Add
  Theorems 3, 4 here when certified.
- `refs/` — Rump 2006, Carstensen papers (local PDFs + txt).

Environment: Windows, PowerShell; python 3.13 with numpy 2.4, scipy
1.17, python-flint (arb). Machine RAM ≈ 16 GB — this budget drives the
mesh decisions below. Long runs: use `python -u`, background them, and
print with flush=True (already done in m3p_certify).

## 3. Mathematics you may take as proved (do not re-derive)

For ANY prime 𝔭 of ℤ[i] with N𝔭 < 4π²Y² ≈ 61.7 (Y = 1.25 throughout):

- Facts A–B hold with right-coset reps {identity} ∪ {γ_k = (0,−1;1,k)}
  where k runs over Gaussian-integer lifts of ℤ[i]/𝔭:
  σ₀⁻¹γ_k = (1,k;0,1) is height-preserving, so the level-𝔭 domain is
  exactly (N𝔭+1) isometric copies of the level-1 reference cell K
  (truncated at Y) plus two product collars T×(Y,∞) and T₀×(Y,∞) at the
  SAME height Y. (0cuspchecks.md C2–C3; the fold-tiling argument only
  needs residues closed under negation ✓.)
- Cusp data: Λ₀ = 𝔭, fold present (C1); |T_∞| = 1/2, |T₀| = N𝔭/2;
  β_∞ = 2(1−s)/Y², β₀ = 2(1−s)/(N𝔭·Y²); t_∞(1) = 1/2, t₀(1) = N𝔭/2;
  shared κ_c(s) = 1/((1+s)Y²) ⇒ constraint rank-one in t_∞ + t₀;
  shortest dual vector of Λ₀ = 1/√N𝔭 ⇒ nonzero-mode bound
  4π²Y²/N𝔭 > 1 ✓ for N𝔭 ∈ {9, 13}.
- Precise invariance of both collars at height Y: inherited from the
  level-1 Shimizu lemma because σ₀ = S ∈ Γ (C6). Prime-independent.
- Theorem G1𝔭 (CONGRUENCE.md §7) with **Lemma D0**: the CR interpolant
  reproduces t_α exactly (top plane tiled by mesh faces + face-mean
  property) ⇒ NO trace error terms, boundary blocks in D_h uninflated,
  σ_h = √6′·α_h^ref where √6′ generalizes to √(N𝔭+1) (disjoint copies —
  check this exponent when editing: the (D3) constant is
  √(#copies)·α_h^ref).
- Exact-main-term refinement: M_ex, Q_rem via Taylor moments
  (`weighted_moments`), pc sandwich only on error terms; c_e carries the
  (1−(1/θ₄−1)ρ_w) factor.
- Rump certificate applies to any standard-model Cholesky execution
  ("any library routine", op. cit.).

Convention warning (0cuspchecks.md C8): σ₀ is NOT unit-covolume
normalized; all formulas carry |T_α| explicitly. Keep it that way.

## 4. Code changes required (small, localized)

1. **Residue ring** (`congruence_prototype.build_gluing` and its use in
   m3p_certify): currently hardcoded 𝔽₅ = ℤ/5 with i ↦ 3 and inverses
   via pow(c,3,5). Generalize to R = ℤ[i]/𝔭:
   - 𝔭 = (3): R = 𝔽₉, elements (a + bi) with a,b ∈ ℤ/3, i² = −1.
     Implement as pairs with ring ops; inverse via (a−bi)/(a²+b²)
     (norm in 𝔽₃, nonzero for nonzero elements since −1 is a
     non-residue mod 3 ⇒ a²+b² = 0 ⟺ a = b = 0 ✓).
   - 𝔭 = (3+2i): R = 𝔽₁₃ prime field, i ↦ 8 (8² = 64 ≡ −1 mod 13).
     Same integer code path as 𝔽₅ with p = 13 and i-image 8.
   - ℙ¹(R): points (0:1) and (1:d), d ∈ R — N𝔭+1 points. Normalization:
     (c:d) with c invertible → (1 : d·c⁻¹).
   - Generator matrices mod 𝔭: T1 = (1,1;0,1); R_rot = (i,0;0,−i);
     T_iR = (1,i;0,1)·R_rot; S = (0,−1;1,0) — map entries through the
     ring embedding. T1 inverse = (1,−1;0,1); the other three are
     involutions in PSL (assert this in the ring: M² ≡ ±identity).
   - Keep every assert: permutation property, round-trip
     act(perm[n], δ) = n, cusp-class sizes {1, N𝔭} (class ∞ ⟺ c ≡ 0
     in R), connectivity of the glued dof graph.
2. **Copy count**: NC = N𝔭 + 1 threaded through assembly (now
   literal 6), and the √6 in σ_h → √(N𝔭+1) (see §3). NP constant
   parametrized (β₀, t₀ checks).
3. **Check values**: t₀(1) = N𝔭/2; 1ᵀM_ex1 = (N𝔭+1)·vol_w(K_h^ref)
   (≈ (N𝔭+1)·0.1417 at the 8×4×3 reference, ≈ ·0.1382 at 6×3×3).

No changes to: reference mesh builder, arb element enclosures, Lemma
G/S constants, window coefficient formulas (other than NP), Rump
routine.

## 5. Hard-won lessons (violate these and you will lose days)

- **Lemma D0 is load-bearing.** Do not reintroduce slab/τ bounds for
  t_α; they cost the certificate at 𝔭 = (2+i). The `tau` field still
  computed in `reference_arb` is legacy/diagnostic — unused in
  coefficients.
- **Never pc-sandwich main terms** — 38% mass inflation killed the
  first attempt. M_ex/Q_rem stay.
- **Parameter search**: scalar slack ANTI-correlates with pencil
  feasibility through ρ̃. The search must keep (a) the c_e ≥ d_e scalar
  conditions, (b) the constant-direction heuristic
  ρ̃(1−θ)ℒ(1)² ≥ D(1), and (c) the float sparse-LOBPCG pencil min-eig
  filter (≥ ~1.15 before attempting Rump). All three are in
  m3p_certify.certify already.
- **Memory**: dense A must be built as ONE Fortran-ordered array,
  filled blockwise from sparse, rank-ones added in 2048-row blocks,
  cholesky(..., overwrite_a=True, check_finite=False). A C-ordered
  array or a .toarray() of the full sum silently doubles memory.
- **Mesh generator invariants**: face-center nodes MUST use canonical
  keys shared between adjacent hexes (a per-hex key cracks the mesh —
  caused a 40-component bug); N₁ must be even (mirror symmetry); floor
  lift carries a (1+1e-9) rounding guard that the arb Lemma-G check
  needs.
- **LOBPCG matvec**: ravel the input (column-matrix × vector broadcast
  bug); tol governs the worst block vector, the reported pair sits
  ~2 orders below.
- Windows/parameters frozen at 𝔭=(2+i): ν* = 1.02, 8 windows,
  (θ,θ₂,α,θ₄,ρ̃) = (0.6, 0.9, 0.2, 0.85, 9.0) — good starting point,
  but let the search rerun per level (ρ̃ demand shrinks with N𝔭).

## 6. Execution order and budgets

Per level (do N𝔭 = 9 first, completely, before touching 13):
1. Generalize the residue ring; run the assert suite + a float margin
   sweep (M𝔭0) via congruence_prototype-style assembly with sparse
   LOBPCG. Decision numbers: expect min μ over λ∈(0,1) around 3–3.5
   for N𝔭 = 9 (vol 3.05) and 2.5–3 for N𝔭 = 13 (vol 4.27), by the
   sublinear trend (level 1: 7.3; N𝔭=5: 4.4). If μ < ~1.5 something is
   wrong — audit before proceeding (compare the assembled checks of §4).
   **STATUS (np9 worktree, 2026-07-11): DONE.** Residue ring
   generalized for (2+i)/(3)/(3+2i); all gluing asserts pass. M𝔭0 at
   6×3×3 for N𝔭=9: min μ ≈ **2.40**, dofs 26532, checks OK (t₀(1)=4.5,
   connected). Slightly below the 3–3.5 extrapolation but above the
   1.5 audit floor. Details: `PROGRESS_LADDER.md`.
2. Certify N𝔭 = 9 at the 6×3×3 reference (10 copies ≈ 26,700 dofs,
   dense A ≈ 5.7 GB — fits; ~4–6 min/window Cholesky, 8 windows).
   **STATUS (2026-07-11): DONE — Theorem 3.** All 8 windows PSD.
   Frozen (θ,θ₂,α,θ₄,ρ̃)=(0.5, 0.9, 0.15, 0.8, 1.5); ν*=1.02;
   26532 dofs; ~38 min. Note: default high-ρ̃ grid failed pencil
   (~0.90); lower ρ̃=1.5 gives pencil≈1.13. See PROOF.md / PROGRESS_LADDER.md.
3. Certify N𝔭 = 13: preferred 6×3×2 reference (≈ 25,000 dofs, 5 GB;
   viable since d_e is now only λ̃(1+1/α)γ² — verify scalar slack at the
   coarser reference first, γ grows with layer height). Fallbacks:
   6×3×3 at 11.2 GB with maximal memory hygiene, or out-of-core blocked
   right-looking Cholesky on np.memmap (admissible under Rump).
   **STATUS (2026-07-11, np9): DONE — Theorem 4.** M𝔭0 μ ≈ 1.45.
   Certify at 8×3×2 (33176 dofs), ν*=1.001, 16 windows, params
   (0.3, 0.98, 0.1, 0.85, 0.5). Requires power-of-two diagonal
   equilibration + per-row radius Rump (unscaled uniform extra fails
   even though A is SPD). All 16 windows PSD. See PROOF.md / PROGRESS.
4. On success: add "Theorem 3/4" sections to PROOF.md (mirror the
   Theorem 2 template — statement, frozen certificate, dependency delta
   = none), update CONGRUENCE.md §8 status, note the μ values.
   **Theorems 3 and 4 done.**

Reporting standard: every certified claim states its frozen parameters
and is reproducible by one command. Scalar checks by rigorous arb
comparison only (bool(x > y) on arb is rigorous); float appears only in
(a) parameter *choice* and (b) inside the Rump-covered Cholesky.

## 7. Open flanks you do NOT need to touch

- CAMWA doi:10.1016/j.camwa.2026.03.029 into refs/ (user access;
  non-load-bearing).
- Composite levels (e.g. (5) = (2+i)(2−i), >2 cusps), LaTeX writeup,
  N𝔭 ≥ 61 mode-bound ceiling — noted in CONGRUENCE.md §8, out of scope.
