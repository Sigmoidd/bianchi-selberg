# Theorem (certified): the Picard orbifold has no exceptional eigenvalues

**Theorem.** Let Γ = PSL(2, ℤ[i]) act on hyperbolic 3-space ℍ³. Then the
Laplace–Beltrami operator on L²(Γ\ℍ³) has **no eigenvalue in (0, 1)**.
Combined with the spectral decomposition (Friedman, Thm 3.8.1; see §5),
the spectrum in (0,1) is empty: **λ₁ ≥ 1**, the Selberg-type bound at
level 1.

This proof is *independent of the Selberg trace formula*: it shares no
input with the engine in the parent repository except the Humbert
fundamental domain. Date certified: **2026-07-10**. Status: complete
modulo the standard citations in §5 (each classical, none numerical).

## 1. Architecture

Three layers, each documented in this directory:

1. **Reduction** (`DESIGN.md` §§0–3, all proofs there): a Lax–Phillips
   cusp analysis reduces the theorem to positivity of an explicit
   quadratic form on H¹(K) for a compact region K ⊂ ℍ³ (truncated
   fundamental domain, Y = 1.25), for every spectral parameter
   s ∈ (0,1): the form 𝒜ₛ on the hyperplane {ℒₛ = 0}, equivalently the
   unconstrained pencil inequality (P_s) of `lower_bound_theory.md` §0.
   Only the cusp zero mode needs exact treatment; all other Fourier modes
   contribute nonnegatively; face identifications of the orbifold are
   relaxed away (testing over all of H¹(K) is a relaxation).

2. **Guaranteed finite-element lower bound** (`lower_bound_theory.md`,
   Theorem G1, all proofs there): a Crouzeix–Raviart discretization on an
   inner polyhedral mesh K_h ⊂ K (floor-lift Lemma G) certifies (P_s) on
   a window of s-values from three *finite* checks: two scalar
   inequalities (c_Σ ≥ 0, c_e > d_e) and one matrix positive-
   semidefiniteness check (N_h − D_h ⪰ 0 on the CR space). The
   interpolation constant is the self-contained Lemma I1
   (κ_sc = √(1/π²+1/15), arbitrary tetrahedra); the sliver between the
   polyhedral floor and the sphere is absorbed by closed-form
   column-mean estimates (Lemma S); the trace functional by the
   slab-mean identity (Lemma E).

3. **Verified certificate** (`m3_certify.py`): the finite checks of
   Theorem G1 executed rigorously — geometry admissibility and all
   constants in ball arithmetic (python-flint/Arb, 128-bit); the ℓ₀
   vector enclosed by a degree-5 Taylor expansion with remainder ball;
   scalar checks by rigorous arb comparison; the matrix check by the
   interval-to-midpoint reduction (Lemma R) plus a literal
   implementation of Rump's positive-definiteness verification
   (Thm 2.3 + Cor 2.4 + Lemma 2.5 of Rump, BIT 46 (2006); IEEE-754
   double, round-to-nearest, underflow included).

## 2. Frozen certificate parameters

    mesh:      12 × 6 × 6 mapped hex cells → 2592 tets, 5544 CR dofs
    Y = 1.25,  ρ = 55,  θ = 0.7,  θ₂ = θ' = α = 0.5,  ν* = 1.05
    s-grid:    8 uniform windows covering [0,1]
    ℓ₀ Taylor order p = 5,  arb precision 128
    κ:         κ_sc = √(1/π² + 1/15)   (Lemma I1; also passes with the
               sharper CZZ20 value √(1/π² + 1/120), kappa_mode="czz")

Certified run (2026-07-10), abridged:

    geometry admissible (Lemma G, arb-verified): True
    constants (arb, upper): gamma=0.09847 tau=0.29085 alpha_h=0.02240
      S_M=1.34e-01 S_Q=7.45e-03 S_S=5.71e-03 V_S=6.94e-03
    radii: ||l_rad||=6.6e-08 ; ||Qr||=2.5e-13 ; ||Mr||,||t_rad|| ≈ 1e-17
    all 8 windows: c_e>d_e True, PSD True (Rump shift c ≈ 2e-08,
      float diagnostic margins ≈ 1e-04)
    ==> M3 CERTIFIED (all windows): True

Reproduce:

    python cr_prototype.py     # float pipeline + box-model validation
    python m3_certify.py       # the certificate (both kappa modes)

## 3. Why each check suffices (one-paragraph summary)

An eigenfunction u with eigenvalue λ = 1−s² ∈ (0,1) is L² and
real-analytic. Its cusp zero mode is exactly c(y/Y)^{1−s}; Lemmas 1–2 of
DESIGN.md turn Q(u) − λM(u) = 0 into 𝒜ₛ(u|_K) ≤ 0, and orthogonality to
constants into ℒₛ(u|_K) = 0, while u|_K ≠ 0 by unique continuation. The
pencil inequality (P_s), verified for all s through the finite window
checks of Theorem G1, gives 𝒜ₛ > 0 on {ℒₛ = 0} — contradiction. The
spectral decomposition upgrades "no eigenvalues in (0,1)" to
"no spectrum in (0,1)".

## 4. What the certificate rests on (complete list of external inputs)

Mathematical:
- **Humbert fundamental domain** for PSL(2,ℤ[i]) (classical; shared with
  the parent repo, the only shared input).
- **Shimizu's lemma** (precisely invariant horoball) — one-line
  self-contained proof in DESIGN.md Lemma 0 (cf. Shimizu, Ann. of Math.
  77 (1963); Leutbecher, Math. Z. 100 (1967)).
- **Spectral decomposition of cofinite Kleinian groups** — Friedman,
  *The Selberg trace formula for PSL(2,O_K)*, arXiv:math/0612807,
  Theorem 3.8.1 (in-repo text: `../friedman_thesis.txt`), resting on
  [EGM98 §6.2]. Used only for: (i) spectrum ∩ [0,1) is atomic
  (Eisenstein part has eigenvalue 1+t² ≥ 1); (ii) the λ₁ ≥ 1 phrasing.
  The core statement "no L² eigenvalues in (0,1)" needs neither.
- **Real-analyticity of eigenfunctions** (Morrey–Nirenberg; supplies
  unique continuation).
- **Payne–Weinberger inequality** for convex domains (Arch. Rational
  Mech. Anal. 5 (1960); corrected proof Bebendorf, Z. Anal. Anwend. 22
  (2003)) — the only input to Lemma I1.
- **Rump, *Verification of positive definiteness*, BIT Numer. Math. 46
  (2006) 433–452** — Thm 2.3, Cor 2.4, Lemma 2.5 (local copy in
  `refs/`). Requires IEEE-754 semantics of the floating Cholesky; the
  paper: "any library routine can be used".

Everything else — cusp Fourier analysis, Hardy/ODE lemmas, Lemma T
(trace identity), Lemma G (floor inclusion), Lemma S (sliver), Lemma E
(functional errors), Lemma I1 (CR constant), Lemma R (interval
rank-one), Theorem G1 — is proved in this directory from scratch.

Software: python-flint (Arb ball arithmetic, correctness-critical),
NumPy/SciPy/LAPACK (dpotrf — covered by Rump's theorem; everything else
diagnostic), IEEE-754 double arithmetic.

Explicitly **not** used: the Selberg trace formula, any test function,
the scattering matrix / ζ_K, the torsion or conjugacy-class inventories,
Then's computed spectrum, and every result file of the parent repo.

## 5. Relation to the trace-formula theorem

The parent repo certifies B(ℤ[i]) ≤ 0.3199 < 1, hence the same
conclusion, via the trace formula. The two proofs share only the
fundamental domain. Agreement of two disjoint methods on λ₁ ≥ 1 is the
independent validation the parent README asked for — and this
architecture (torsion-free congruence covers, where the elliptic
bookkeeping vanishes and only the mesh grows) is reusable for levels
where the statement is new.

## Theorem 2 (level 𝔭 = (2+i)): the first congruence rung

**Theorem 2.** Let Γ₀(𝔭) ⊂ PSL(2,ℤ[i]) be the Hecke congruence subgroup
for 𝔭 = (2+i), N𝔭 = 5 (index 6, two cusps ∞ and 0, vol ≈ 1.832). Then
the Laplacian on L²(Γ₀(𝔭)\ℍ³) has **no eigenvalue in (0, 1)**.
Certified 2026-07-10. Unlike Theorem 1, this statement is not reachable
by the trace-formula positivity method (B ≈ 1.9 > 1 at this volume).

**Architecture** (documents in this directory):
- Multi-cusp criterion: `CONGRUENCE.md` §2 — per-cusp copies of the
  level-1 Lemmas 0–3, with per-cusp data (Λ_α, |T_α|, Y_α = Y); the
  cusp-0 chain is audited item-by-item in `0cuspchecks.md` (conjugated
  stabilizer, fold-compatible tiling, height-preserving chimney maps,
  dual-lattice mode bound, inherited precise invariance).
- Domain structure: exactly 6 isometric copies of the level-1 reference
  cell plus two product collars at the same height Y (Facts A–B,
  `CONGRUENCE.md` §1) — the Reference-Cell Principle: all interval data
  is level-1-sized; the index enters only through exact 𝔽₅ combinatorics
  and the size of the final linear algebra.
- Guaranteed lower bound: Theorem G1𝔭 (`CONGRUENCE.md` §7) = Theorem G1
  with three refinements, each stated and proved there: **Lemma D0**
  (the CR interpolant reproduces the cusp trace functionals *exactly* —
  no trace error terms, uninflated boundary blocks; also strengthens
  Theorem 1 retroactively, see the remark in `lower_bound_theory.md`),
  **exact main terms** (Taylor-enclosed exact-weight mass M_ex and
  stiffness remainder Q_rem; the weight sandwich survives only in error
  terms), and per-copy sliver constants with per-face column heights.
- Certificate: `m3p_certify.py` — same verification standard as
  Theorem 1 (arb geometry admissibility, degree-5 Taylor enclosures,
  rigorous arb scalar comparisons, literal Rump BIT 46 (2006) matrix
  certificate).

**Frozen certificate.** Reference mesh 8×4×3 (24-split; 2304 tets/copy),
28,400 global CR dofs; Y = 1.25; ν* = 1.02; (θ, θ₂, α, θ₄, ρ̃) =
(0.6, 0.9, 0.2, 0.85, 9.0); κ = κ_sc (Lemma I1); Taylor p = 5; arb
prec 128; 8 uniform s-windows. Result: every window passes with
c_e ≥ 0.847 vs d_e ≤ 0.138 and Rump-verified N_h − D_h ⪰ 0 (float
pencil diagnostic min-eig 1.44). Reproduce: `python m3p_certify.py`.

**Dependency delta vs Theorem 1: none.** The external-input list of §4
is unchanged — the additional ingredients (coset combinatorics, Facts
A–B, Lemma D0, C1–C6 of `0cuspchecks.md`) are all proved in-repo;
Friedman Thm 3.8.1 already covers cofinite Kleinian groups, hence
Γ₀(𝔭) with trivial character verbatim.

## Theorem 3 (level 𝔭 = (3)): the inert rung N𝔭 = 9

**Theorem 3.** Let Γ₀(𝔭) ⊂ PSL(2,ℤ[i]) be the Hecke congruence subgroup
for 𝔭 = (3) (inert, N𝔭 = 9, index 10, two cusps ∞ and 0). Then the
Laplacian on L²(Γ₀(𝔭)\ℍ³) has **no eigenvalue in (0, 1)**.
Certified 2026-07-11 (worktree np9).

**Architecture:** same as Theorem 2. Residue ring is ℤ[i]/(3) ≅ 𝔽₉
(pairs a+bi mod 3); Facts A–B give 10 isometric copies of the level-1
reference cell; Theorem G1𝔭 with Lemma D0 uses σ_h = √10 · α_h^ref.
Float margin M𝔭0 at the same mesh: min μ ≈ 2.40 on (0,1). Full log:
`PROGRESS_LADDER.md`.

**Frozen certificate.** Reference mesh 6×3×3 (24-split; 1296 tets/copy),
26 532 global CR dofs; Y = 1.25; ν* = 1.02; (θ, θ₂, α, θ₄, ρ̃) =
(0.5, 0.9, 0.15, 0.8, 1.5); κ = κ_sc (Lemma I1); Taylor p = 5; arb
prec 128; 8 uniform s-windows. Every window: c_S>0, c_e>d_e, Rump PSD.
Window extremes: c_e ∈ [0.836, 0.881], d_e ∈ [0.079, 0.222]. Float
pencil diagnostic min-eig ≈ 1.13 at the search windows. Wall time
≈ 38 min. Reproduce (from `independent_exclusion/`):
```text
python -u -c "from m3p_certify import certify; certify(6,3,3,params=(0.5,0.9,0.15,0.8,1.5),level='(3)')"
```
or after a parameter re-search: `python -u m3p_certify.py coarse "(3)"`.

**Dependency delta vs Theorem 2: none.** Only the residue ring and
copy count change; criterion, G1𝔭, Rump standard unchanged.

## Theorem 4 (level 𝔭 = (3+2i)): the split rung N𝔭 = 13

**Theorem 4.** Let Γ₀(𝔭) ⊂ PSL(2,ℤ[i]) be the Hecke congruence subgroup
for 𝔭 = (3+2i) (split, N𝔭 = 13, index 14, two cusps ∞ and 0). Then the
Laplacian on L²(Γ₀(𝔭)\ℍ³) has **no eigenvalue in (0, 1)**.
Certified 2026-07-11 (worktree np9).

**Architecture:** same as Theorems 2–3. Residue ring ℤ[i]/𝔭 ≅ 𝔽₁₃
with i ↦ 8; 14 isometric copies of the level-1 reference cell;
σ_h = √14 · α_h^ref (Lemma D0). Float margin M𝔭0: min μ ≈ 1.45 on
(0,1) (mesh-converged). Full log: `PROGRESS_LADDER.md`,
`theorem4postmortem.md`.

**Rump implementation note (load-bearing for this size).** Before the
BIT 46 test, the midpoint matrix is reduced by the exact power-of-two
diagonal congruence A ↦ SAS with
s_i = 2^{round(−½ log₂ a_ii)} (IEEE-exact; A ⪰ 0 ⟺ SAS ⪰ 0).
Radii are applied **per row** (r_i = Σ_j |S ΔA S|_ij subtracted from
diagonal i), not as a uniform max-row-sum shift — the latter exceeds
λ_min(A) ~ 10⁻⁵ at n ≈ 3·10⁴ while the matrix is still SPD. No change
to Theorem G1𝔭.

**Frozen certificate.** Reference mesh 8×3×2 (24-split; 1152 tets/copy),
33 176 global CR dofs; Y = 1.25; ν* = 1.001; (θ, θ₂, α, θ₄, ρ̃) =
(0.3, 0.98, 0.1, 0.85, 0.5); κ = κ_sc (Lemma I1); Taylor p = 5; arb
prec 128; **16** uniform s-windows. Every window: c_S>0, c_e>d_e,
Rump PSD after equilibration. Window extremes: c_e ∈ [0.870, 0.907],
d_e ∈ [0.074, 0.485]. Wall time ≈ 93 min. Reproduce
(from `independent_exclusion/`):
```text
python -u -c "import m3p_certify as m; m.NU_STAR=1.001; m.NWIN=16; m.certify(8,3,2,params=(0.3,0.98,0.1,0.85,0.5),level='(3+2i)')"
```

**Dependency delta vs Theorem 3: none** on the mathematical chain;
execution uses the equilibrated Rump path above (implementation-only).

## 6. Loose ends (tracked, non-load-bearing)

- CAMWA doi:10.1016/j.camwa.2026.03.029 (Carstensen–Dond–Maity–Nataraj)
  into `refs/` when accessible — motivational/primary-LEB citation only;
  no constant traces to it.
- Optional sharpenings: CZZ20's κ₁ (already wired in), per-element
  interpolation constants (Liu 2015), finer mesh for larger certified
  margins.
- Writeup for publication (LaTeX) — this file is the spine.
