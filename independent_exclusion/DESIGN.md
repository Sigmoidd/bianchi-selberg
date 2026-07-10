# Independent exclusion of exceptional eigenvalues (Target B)

Goal: prove **λ₁(PSL(2,ℤ[i])\ℍ³) ≥ 1** by a method that never touches the
trace-formula engine. Architecture: Lax–Phillips cusp reduction + certified
finite-element lower bounds on the compact core.

**Inputs used:** the Humbert fundamental domain, Shimizu's horoball lemma,
standard spectral theory of cofinite Kleinian groups (continuous spectrum
[1,∞), discreteness below 1), a 1D ODE in the cusp, and verified linear
algebra. **Inputs NOT used:** the Selberg trace formula, any test function
h, the scattering matrix, the torsion/conjugacy-class inventory, Then's
spectrum, and every file in the parent repo. This is the independence the
README asks for.

## 0. Conventions

ℍ³ = {P = (z,y) : z = x₁+ix₂ ∈ ℂ, y > 0}, metric (|dz|²+dy²)/y²,
volume dV = y⁻³ dx₁dx₂dy. For u real,

    Q(u) = ∫ |∇u|²_hyp dV = ∫ (|∇_euc u|²) y⁻¹ dx₁dx₂dy,
    M(u) = ∫ u² dV.

Γ = PSL(2,ℤ[i]). Humbert fundamental domain
F = { (z,y) : x₁∈[−½,½], x₂∈[0,½], |z|²+y² ≥ 1 }.
Γ∞ ≅ ℤ[i] ⋊ ⟨z↦−z⟩; cusp cross-section orbifold T = coordinate rectangle
x₁∈[−½,½], x₂∈[0,½], Euclidean area |T| = ½.

Fix a truncation height **Y > 1** (strictly: this removes the mesh
degeneracy where the unit sphere touches y=1 at z=0). Split

    K = F ∩ {y ≤ Y}   (compact core: y_f(z) ≤ y ≤ Y, y_f = √(1−|z|²)),
    C = F ∩ {y > Y}   (cusp: T × (Y,∞)).

## 1. Lemmas

**Lemma 0 (precisely invariant horoball, Shimizu).** For γ = (a b; c d) ∈ Γ
with c ≠ 0: y(γP) = y / (|cz+d|² + |c|²y²) ≤ 1/(|c|²y) ≤ 1/y since
|c|² ≥ 1 for 0 ≠ c ∈ ℤ[i]. Hence {y > Y} with Y ≥ 1 meets its Γ-translates
only through Γ∞, and the cusp Fourier analysis below is legitimate on the
quotient. ∎

**Spectral facts (pinned citations).** The spectral decomposition of
L²(Γ\ℍ³) for cofinite Kleinian groups: Friedman, *The Selberg trace
formula for PSL(2,O_K)* (arXiv:math/0612807, in-repo
`friedman_thesis.txt`), **Theorem 3.8.1** — every f in the domain expands
into discrete eigenfunctions e_m plus Eisenstein integrals E_α(·, it);
the underlying eigenpacket theory is [EGM98, §6.2] (Elstrodt–Grunewald–
Mennicke, *Groups Acting on Hyperbolic Space*, Ch. 6). Since E_α(·, it)
has Laplace eigenvalue 1 + t² ≥ 1, the spectral measure on [0,1) is
purely atomic, supported on the discrete eigenvalues: **any spectral
point in (0,1) is an L² eigenvalue.** Eigenfunctions are real-analytic
(elliptic operator with real-analytic coefficients; Morrey–Nirenberg —
this also supplies unique continuation, in place of Aronszajn). λ = 0 has
eigenspace = constants: Δu = 0 with u ∈ L² on a connected finite-volume
quotient forces ∫|∇u|² = ⟨Δu, u⟩ = 0. Eigenfunctions with λ ≠ 0 are
orthogonal to constants (self-adjointness). Lemma 0 below is Shimizu's
lemma (cf. Shimizu, Ann. of Math. 77 (1963); Leutbecher, Math. Z. 100
(1967)); the one-line proof given is self-contained.

Let u be a real eigenfunction, Δu = λu, λ = 1−s², s ∈ (0,1), ‖u‖ = 1.
Work on the double cover of the cusp, C̃ = (ℝ²/ℤ²) × (Y,∞) (orbifold
integrals = ½ · cover integrals). Fourier modes û_μ(y), μ = (m,n) ∈ ℤ².

**Lemma 1 (nonzero modes are harmless).** For μ ≠ 0, the mode's
contribution to Q_C − λM_C is

    ½ ∫_Y^∞ ( û_μ′² + 4π²|μ|² û_μ² ) y⁻¹ dy − ½ λ ∫_Y^∞ û_μ² y⁻³ dy
      ≥ ½ (4π²Y² − λ) ∫_Y^∞ û_μ² y⁻³ dy ≥ 0,

using ∫ û² y⁻¹ ≥ Y² ∫ û² y⁻³ and 4π²Y² > 1 > λ. No boundary condition at
y = Y is needed. ∎

**Lemma 2 (exact zero-mode profile and defect).** û₀ satisfies
y²f″ − yf′ + λf = 0 on (Y,∞); the solutions are y^{1±s}, and L²(y⁻³dy)
forces û₀(y) = c·(y/Y)^{1−s} with c = û₀(Y) = (1/|T|)∫_T u(·,Y) dx.
Integration by parts with the ODE gives *exactly*

    ∫_Y^∞ (û₀′² − λ û₀² y⁻²) y⁻¹ dy = − û₀(Y) û₀′(Y)/Y = − (1−s) c² / Y²,

so the zero mode contributes −(1−s)c²/(2Y²) to Q_C − λM_C (½ from the
cover). ∎

**Lemma 3 (orthogonality constraint).** ⟨u,1⟩ = 0 and
∫_Y^∞ (y/Y)^{1−s} y⁻³ dy = Y⁻²/(1+s) give

    ∫_K u dV + c / ( 2(1+s) Y² ) = 0. ∎

## 2. The exclusion criterion

For v ∈ H¹(K) define t(v) = ∫_T v(·,Y) dx₁dx₂ (Euclidean trace integral on
the top face; c = 2t), and for s ∈ (0,1), λ = 1−s²:

    𝒜_s(v) = Q_K(v) − λ M_K(v) − 2(1−s) t(v)² / Y²,
    ℒ_s(v) = ∫_K v dV + t(v) / ( (1+s) Y² ).

**Theorem (criterion).** If for every s ∈ (0,1),

    𝒜_s(v) > 0   for all 0 ≠ v ∈ H¹(K) with ℒ_s(v) = 0,

then Γ\ℍ³ has no eigenvalue in (0,1); i.e. λ₁ ≥ 1.

*Proof.* Let u be an eigenfunction with λ = 1−s² ∈ (0,1), ‖u‖=1. Then
0 = Q(u) − λM(u) = [Q_K − λM_K](u) + (zero mode) + (μ≠0 modes)
≥ 𝒜_s(u|_K) by Lemmas 1–2, so 𝒜_s(u|_K) ≤ 0. Lemma 3 gives
ℒ_s(u|_K) = 0. Real-analyticity (unique continuation) gives u|_K ≠ 0.
Contradiction. ∎

Remarks.
- Testing over **all** of H¹(K) is a relaxation: no side/floor
  identifications, no orbifold structure, plain Neumann-type free
  boundaries. If positivity survives this relaxation (prototype says it
  does — see below), the certified computation never has to implement the
  face gluings.
- The criterion sees *every* L² eigenfunction — cuspidal or residual —
  exactly like the B < 1 positivity argument.
- Constants sanity check: v ≡ 1 has 𝒜_s < 0 but ℒ_s(1) ≠ 0, so the
  constraint is what excludes the λ=0 eigenfunction, as it must.
- Monotonicity: for fixed v, 𝒜_s(v) is decreasing in λ (both −λM_K and
  −2(1−s)t² decrease), so a finite λ-grid checked at right endpoints
  controls each window; the s-dependence of the *constraint* is handled by
  the interval-κ reduction below.

## 3. λ-sweep rigor (finite verification covers the continuum of s)

Write ℒ_s(v) = a(v) + κ t(v), a(v) = ∫_K v dV, κ = 1/((1+s)Y²) ∈ [κ₁,κ₂]
on an s-window. To certify 𝒜 ≻ 0 on ∪_{κ∈[κ₁,κ₂]} ker(a + κt):
prove 𝒜 ≻ ε on W = ker(a) ∩ ker(t) (independent of s), then bound the
2×2 Schur complement on a complement span{w₁,w₂} with interval κ. Both
steps are finite-dimensional after discretization and interval-friendly.

## 4. From criterion to proof: certified FEM lower bounds

The criterion is an infimum over infinite-dimensional H¹(K). A conforming
Galerkin computation bounds that infimum from ABOVE (Rayleigh–Ritz), so
**discrete positivity is evidence, not proof**. The certification stage
needs guaranteed lower bounds:

1. **Nonconforming elements with explicit constants** (Crouzeix–Raviart /
   enriched CR; Liu's framework [Liu 2015, *Appl. Math. Comput.*; also
   Carstensen–Gedicke, Cancès et al.]): guaranteed lower bounds
   λ_k ≥ λ_k^h / (1 + C_h² λ_k^h) with C_h from per-element interpolation
   constants. Needs the weighted forms: sandwich y⁻¹ ∈ [Y⁻¹, y_min⁻¹] and
   y⁻³ per element (weights are monotone in y ⇒ trivial interval bounds on
   each element's y-range), reducing to a piecewise-constant-coefficient
   problem that the framework handles. Fine meshes make the sandwich loss
   negligible.
2. **Interval assembly**: matrix entries as balls (the integrands are
   explicit; element maps are polynomial except y_f = √(1−|z|²), which is
   analytic and interval-evaluable; or remesh with flat-faced tetrahedra
   below a spherical cap handled by a monotone bound).
3. **Verified positive definiteness**: prove A − εM ≻ 0 by interval
   Cholesky/LDLᵀ (python-flint or a small dedicated verified routine).
   This replaces "smallest eigenvalue ≈ μ > 0" with a proof.
4. **Window sweep** over λ per §3.

### Flagged gaps / to-do (analogue of RIGOR_GAPS.md)

- [x] G1 ✅ (2026-07-10): weighted CR lower-bound theorem with explicit
      constants — **Theorem G1 in `lower_bound_theory.md`**, implemented in
      `cr_prototype.py`. Pencil formulation (P_s) removes the constraint;
      per-tet weight sandwiching; CR projection with κ₁ = √(1/π²+1/120)
      [CP22 §2.4]; slab-mean bound for the trace functional; column-mean
      sliver absorption. Float run: **all 8 s-windows pass at ν* = 1.05 on a
      12×6×6 mesh (2592 tets / 5544 CR dofs)** — c_e/d_e ≥ 1.88, PSD checks
      hold, sliver constants ≤ 0.14.
- [x] G3 ✅ folded into Lemma E (slab-mean form) of `lower_bound_theory.md`.
- [x] G5 ✅ folded into Theorem G1 (worst-case window data λ⁺, β⁺ + Δκ
      Young-shift; monotonicity does the rest).
- [ ] G2 (reduced): mesh geometry is now polyhedral-inner with proven
      inclusion (Lemma G: concavity + sag/8 lift); remaining work is only
      the *interval* evaluation of per-tet data and the ℓ₀ vector in M3.
- [ ] G4: citations to pin: EGM Ch. 6 for spectral facts; Aronszajn for
      unique continuation (or real-analyticity of eigenfunctions);
      Shimizu/Leutbecher for Lemma 0. NEW: [CITE-CHECK] κ₁ hypotheses
      (arbitrary tets) in CZZ20/CP20; fallback: per-element constants
      (Liu 2015). NEW: [TODO] CAMWA doi:10.1016/j.camwa.2026.03.029 PDF
      into `refs/` (paywalled; user access needed).

## 5. Prototype (this directory)

`fem_prototype.py` — trilinear Q1 conforming FEM on the mapped box
(ξ₁,ξ₂,ξ₃) ↦ (x₁, x₂, y_f(1−ξ₃) + Yξ₃), dense assembly of the weighted
stiffness S, mass M, volume functional a, top-trace functional t; for each
λ on a grid it computes μ(λ) = min of 𝒜_s on {ℒ_s = 0} as a constrained
generalized eigenvalue (Householder projection). Self-checks: S·1 = 0,
1ᵀM1 = vol(F) − 1/(4Y²) with vol(F) = 0.30532, t(1) = ½, refinement
stability. Purpose: measure the margin min_λ μ(λ) and decide Y; go/no-go
for the Neumann relaxation. It proves nothing by itself.

### M0 results (2026-07-09)

Margin is large and stable. With Y = 1.25 (meshes 16×8×6 and 24×12×9 agree
to 3 digits; vol(K) reproduced to 6 digits; S·1 = 0; t(1) = ½):

    mu(lambda) decreases from 8.23 (lambda=0.05) to 7.28 (lambda=0.999);
    min margin ~ 7.3.   Y = 1.5 gives ~ 8.9.

The minimizer is the first nonconstant Neumann mode of the core
(lambda_2^N ≈ 8.28), which has t(v) ≈ 0 and mean ≈ 0 by parity, so
mu(lambda) ≈ lambda_2^N − lambda. Consequences:
- The Neumann relaxation (no face identifications) is comfortably enough.
- The certification target is modest: any guaranteed lower bound
  ≥ lambda + epsilon (i.e. ≈ 1) for the constrained problem suffices —
  we have ~7 units of headroom for CR constants, weight sandwiching, and
  interval slop.
- Negative control: without the constraint the near-constant direction
  gives −3.39 (pure constant: −3.13), confirming boundary term and
  constraint are load-bearing and correctly signed.

## 6. Milestones

- **M0** ✅ (2026-07-09): prototype margin measured (~7.3), Y = 1.25 chosen
  (1.5 also viable), Neumann relaxation validated.
- **M2** ✅ (2026-07-10): guaranteed-lower-bound theory written and proved
  (`lower_bound_theory.md`, Theorem G1) and implemented (`cr_prototype.py`).
  Float pipeline passes all s-windows at ν* = 1.05 on 2592 tets; box model
  problem validates GLB ≤ exact.
- **M3** ✅ (2026-07-10): interval certification (`m3_certify.py`, arb,
  prec=128). Geometry admissibility verified in arb (Lemma G); all Lemma
  E/S constants as balls; ℓ₀ enclosed by degree-5 Taylor with rigorous
  remainder (‖ℓ_rad‖ ≈ 7·10⁻⁸); scalar window checks by rigorous arb
  comparison; matrix checks via Lemma-R factored rank-one reduction +
  verified float Cholesky (Higham/Rump error model, 10× guard). **All 8
  windows certified**; certified λ_min(A) > 0 with float margins ~10⁻⁴.
  ⇒ **λ₁(PSL(2,ℤ[i])\ℍ³) ≥ 1 certified**, modulo the citation-level
  residue below.
- **M3 hardening** ✅ (2026-07-10): both former [CITE-CHECK]s discharged.
  (i) κ dependency on CZZ20 removed: self-contained **Lemma I1**
  (lower_bound_theory.md), κ_sc = √(1/π²+1/15) ≤ 0.410, arbitrary tets,
  inputs Payne–Weinberger 1960 / Bebendorf 2003; certificate passes with
  both κ_sc (default) and the sharper CZZ20 value (optional).
  (ii) PSD certificate now a literal implementation of Rump, BIT 46
  (2006) Thm 2.3 + Cor 2.4 + Lemma 2.5 (`rump_psd_certificate`), incl.
  eps = 2⁻⁵³, eta = 2⁻¹⁰⁷⁴ underflow terms and the ϕ-trick diagonal;
  Rump op. cit.: "any library routine can be used" (LAPACK dpotrf, as in
  INTLAB isspd). Rump shift c ≈ 2·10⁻⁸ vs float margins ~10⁻⁴.
- **M4** (writeup) remaining residue:
  (iii) G4 spectral-theory citations (EGM Ch. 6, Shimizu, Aronszajn).
  (iv) CAMWA 2026 paper into refs/ (user access), reconcile as primary
  LEB citation (optional — no longer load-bearing for any constant).
- **M1**: mesh refinement study; polyhedral floor treatment (G2 decision).
- **M2**: weighted lower-bound theory (G1, G3) written with proofs.
- **M3**: interval assembly + verified factorization; full λ sweep (G5).
- **M4**: writeup. Result: λ₁ ≥ 1 with no trace-formula input — the
  independent validation of the main theorem, and reusable machinery for
  torsion-free congruence covers (which would ALSO be new theorems).
