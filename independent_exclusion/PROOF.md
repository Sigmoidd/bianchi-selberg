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

## 6. Loose ends (tracked, non-load-bearing)

- CAMWA doi:10.1016/j.camwa.2026.03.029 (Carstensen–Dond–Maity–Nataraj)
  into `refs/` when accessible — motivational/primary-LEB citation only;
  no constant traces to it.
- Optional sharpenings: CZZ20's κ₁ (already wired in), per-element
  interpolation constants (Liu 2015), finer mesh for larger certified
  margins.
- Writeup for publication (LaTeX) — this file is the spine.
