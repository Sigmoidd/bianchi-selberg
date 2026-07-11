# Congruence levels: the open problem, by reference cells

Target: **Γ₀(𝔭) ⊂ PSL(2, ℤ[i]), 𝔭 = (2+i), N𝔭 = 5** — index 6, two cusps
(∞ and 0), vol = 6·0.30532 = 1.832. Goal: certify **no eigenvalues in
(0,1)** for Γ₀(𝔭)\ℍ³ with the Target-B machinery. Unlike level 1 this is
open/publishable territory (parent README roadmap). The trace-formula
route degrades linearly in volume (B ≈ vol-linear); this route's margin
tracks the actual spectral gap and only the *linear algebra* grows.

## 1. The Reference-Cell Principle (scale protection)

**Principle.** In the level-𝔭 certificate, every interval-enclosed
quantity is either (a) a per-cell quantity of the *level-1 reference cell*
K (element matrices, ℓ₀ Taylor enclosures, Lemma E/G/S constants —
computed once, radii independent of the index), or (b) *exact integer
combinatorics* (coset action, gluing table, cusp classes). The only
computation that grows with the index is the verified positive-
definiteness check (Rump), whose cost is linear algebra, not error.

This holds because of two exact structural facts, proved below for
𝔭 = (2+i) (both generalize):

**Fact A (exact tiling by isometric copies).** With right-coset
representatives {γᵢ} of Γ₀(𝔭)\Γ ≅ ℙ¹(𝔽₅) (six points; 𝔽₅ = ℤ[i]/(2+i),
i ↦ 3), D = ⋃ᵢ γᵢF is a fundamental domain for Γ₀(𝔭), and each γᵢF is
isometric to F. Assembling each copy in its own pulled-back reference
coordinates makes all per-copy data equal to the reference data; only
the boundary identifications know the embedding.

**Fact B (collars assemble at the same height).** Truncate every copy at
its reference height Y (= 1.25). The cusp-∞ collar is T × (Y,∞) exactly
as at level 1 ({y>1} is precisely invariant already under Γ∞ ⊂ Γ₀).
For cusp 0 take σ₀ = (0,−1;1,0) and reps γ_k = (0,−1;1,k), k = 0..4
(bottom rows (1:k)): then

    σ₀⁻¹ γ_k = (1, k; 0, 1),

a *height-preserving* translation. Hence in cusp-0 normalized
coordinates the five chimneys γ_k({y>Y}×T) are (T+k) × (Y,∞), and since
{0,…,4} represent ℤ[i]/(2+i), the translates T+k tile the cusp-0
cross-section T₀ (fundamental domain of (2+i)ℤ[i] ⋊ ⟨−1⟩, area 5/2).
So the assembled domain is **exactly**

    D_Y = [6 isometric copies of K]  ∪  T×(Y,∞)  ∪  T₀×(Y,∞).

Volume check: 6(0.30532 − 1/(4Y²)) + (1/2 + 5/2)/(2Y²) = 6·0.30532 ✓.

## 2. Multi-cusp criterion (constants for Γ₀(2+i))

Both cusp stabilizers contain the order-2 rotation (diag(i,−i) is
c ≡ 0, and its σ₀-conjugate is again diagonal), so both cross-sections
are folded, exactly as at level 1. Per-cusp data:

| cusp α | lattice Λ_α | \|T_α\| | shortest dual \|μ\| | mode bound (need > λ) | Y_α |
|---|---|---|---|---|---|
| ∞ | ℤ[i] | 1/2 | 1 | 4π²Y² ≈ 61.7 | 1.25 |
| 0 | (2+i)ℤ[i] | 5/2 | 1/√5 | 4π²Y²/5 ≈ 12.3 | 1.25 |

With t_α(v) = ∫_{T_α} v(·,Y) dx (coordinate integral; c_α = t_α/|T_α|),
the level-1 Lemmas 1–3 of DESIGN.md apply verbatim per cusp:

- zero-mode defect of an eigenfunction: −Σ_α |T_α|(1−s) c_α²/Y²
  = −(1−s)/Y² · [ 2 t_∞² + (2/5) t₀² ]; so β_∞(s) = 2(1−s)/Y²,
  β₀(s) = (2/5)(1−s)/Y²;
- nonzero modes ≥ 0 by the table's mode bounds (both ≫ 1);
- orthogonality to constants: ℒ_s(v) = a(v) + κ_c(s)(t_∞(v) + t₀(v)),
  κ_c(s) = 1/((1+s)Y²)  — the same κ_c for both cusps since Y_∞ = Y₀,
  so the constraint stays rank-one with the *combined* vector t_∞ + t₀.

**Criterion (level 𝔭).** If for all s ∈ (0,1) and all v ∈ H¹(D_Y-core):
Q(v) − (1−s²)M(v) − β_∞ t_∞(v)² − β₀ t₀(v)² > 0 subject to ℒ_s(v) = 0,
then Γ₀(𝔭)\ℍ³ has no eigenvalue in (0,1). Pencil form, Theorem-G1
window machinery, Lemma R: identical, with the rank-one boundary term
replaced by rank-two (per-vector treatment; both t-vectors are exact in
the same sense as level 1: face areas at y = Y).

Slab bounds (E-t) for t₀: the cusp-0 slab pulls back to the *reference*
top slabs of the five class-0 copies (Fact B's maps are height-
preserving), so τ₀ uses level-1 per-cell constants with a √5-area factor
in the Cauchy–Schwarz — per-cell data unchanged (Reference-Cell
Principle).

## 3. Gluing combinatorics (exact, 𝔽₅)

Reference-cell side pairings of F and their face maps:

| δ | matrix | faces | reference face map |
|---|---|---|---|
| T₁ | (1,1;0,1) | x₁=−½ → x₁=+½ | identity in (x₂,y) |
| R | (i,0;0,−i) ≡ (3,0;0,2) | x₂=0 → itself | x₁ → −x₁ |
| T_iR | (1,i;0,1)·R ≡ (3,1;0,2) | x₂=½ → itself | x₁ → −x₁ |
| S | (0,−1;1,0) ≡ (0,4;1,0) | floor → itself | (x₁,x₂) → (−x₁,x₂) |

Cosets = ℙ¹(𝔽₅) via bottom row (c:d); right action (c:d)·M. Copy c's
δ-source face glues to copy j := (c:d)·δ⁻¹'s δ-target face via the
reference face map (derivation: γ_c x ~ γ_j (δx) iff γ_c δ⁻¹ γ_j⁻¹ ∈
Γ₀(𝔭)). If j = c the identification is internal to one copy and is
**relaxed** (Neumann), exactly as at level 1 — valid, it enlarges the
test space. If j ≠ c the interface is *interior to the domain* and must
be **conforming**: relaxing it would hand each copy its own constant
function and the criterion would (correctly) fail. Cusp class of coset
(c:d): ∞ iff c ≡ 0; else 0. Sizes 1 and 5 ✓.

## 4. Mesh: the symmetric 24-split

Conformity across interfaces requires the reference boundary
triangulation to be compatible with the face maps: identity/translation
across the x₁-faces AND the x₁-mirror on x₂=0, x₂=½ and the floor. Kuhn
splits (uniform diagonals) violate the mirror maps. Resolution: the
**24-split** — each mapped hex gets a body-center node and six
face-center nodes; each quad face → 4 triangles; each triangle + body
center → tet. This split is invariant under all face maps in the table
(the 4-way face split has no diagonal convention) and automatically
face-to-face across translated neighbors. Cost: 4× the tets of Kuhn;
compensated by a coarser base grid. Lemma G (floor lift) reproves
verbatim on the quarter-triangles (smaller sag, if anything); Lemmas
E/S/I1/T and Theorem G1 are mesh-split-agnostic. Base grid needs N₁
even (mirror-symmetric grid lines).

## 5. Verification-scaling plan

- Per-cell interval data: level-1 sized, computed once (Principle).
- Gluing: exact integer, asserts (permutation property, involutivity,
  class sizes, connectivity of the assembled dof graph, 1ᵀM1 = 6·vol_w,
  t_∞ᵀ1 = ½, t₀ᵀ1 = 5/2).
- PSD checks: matrices are sparse + rank ≤ 3 (ℓ, t_∞, t₀). Rump Thm 2.3
  is sparsity-aware (Δ(A) built from the symbolic factor); path: sparse
  Cholesky with fill-reducing ordering + the same Cor 2.4/Lemma 2.5
  test, or dense on a downsized certified mesh if margins allow. This is
  the only index-growing cost.
- Expected margin: smaller than level-1's ≈ 7 (bigger core ⇒ smaller
  first nonconstant "Neumann" eigenvalue; the cusp-difference direction
  t_∞ vs t₀ is the new low-energy candidate that the assembled pencil
  must beat). Measured by `congruence_prototype.py` (M0-analog): go/no-go
  below.

## 6. Status

- [x] **M𝔭0 ✅ (2026-07-10): margin measured and mesh-converged —
      μ(λ) ≥ 4.39 on (0,1).** Coarse 4×2×2 (4816 dofs, dense):
      min μ = 4.457; refined 6×3×3 (16,020 dofs, low-rank/sparse solve,
      run by user): min μ = 4.393 — a −1.4% drift, conforming
      upper-bias as expected, continuum value ≈ 4.3–4.4. μ decreases
      6.66 → 4.39 as λ → 1; second eigenvalue ≈ 10.6 well separated.
      All structural checks pass at both resolutions: gluing graph
      connected, t_∞(1) = ½ and t₀(1) = 5/2 exact, Q·1 = 0,
      weighted-volume deficit vs 6·vol_w(K) shrinking O(h²)
      (0.111 → 0.043), consistent with the sliver analysis. The repo
      file now carries a sparse LOBPCG path (constraint via the
      M-orthogonal complement, boundary terms as low-rank operator
      corrections), cross-validated against dense on the coarse mesh.
      **Scaling observation:** level-1 margin 7.3 → level-𝔭 margin 4.5
      at 6× the volume — sublinear degradation, versus the trace-formula
      B ≈ 0.32·6 ≈ 1.9 > 1 which fails outright here. This quotient is
      only reachable by the Target-B machinery.
- [x] M𝔭1 ✅ (2026-07-10): Theorem G1𝔭 below (§7).
- [x] **M𝔭2 ✅ (2026-07-10): CERTIFIED — Γ₀(2+i)\ℍ³ has no eigenvalue
      in (0,1).** All 8 s-windows pass in `m3p_certify.py`: scalar
      checks by rigorous arb comparison (c_e ≥ 0.847 vs d_e ≤ 0.138 —
      slack ≥ 6×), matrix checks by the literal Rump certificate.
      **Frozen certificate:** reference mesh 8×4×3 (24-split, 2304
      tets/copy), 28,400 global CR dofs, Y = 1.25, ν* = 1.02,
      (θ, θ₂, α, θ₄, ρ̃) = (0.6, 0.9, 0.2, 0.85, 9.0), κ = Lemma I1,
      Taylor p = 5, arb prec 128. Float pencil min-eig 1.44 (diagnostic).
      The decisive step was **Lemma D0** (exact trace reproduction by
      I_CR): it removed τ and the boundary Young inflations, taking the
      hardest window from PSD-fail/scalar-knife-edge to 6× slack. Also
      load-bearing: exact-main-term refinement (M_ex, Q_rem — the pc
      sandwich alone cost 38% mass inflation), and the 8×4×3 mesh
      (halving δ̂ to get S_M = 0.27).
      Precise invariance at both cusps at height Y is inherited from
      level 1: {y > 1} is precisely invariant under the full Γ (Lemma 0),
      and both scaling matrices σ_α ∈ Γ, so the horoballs at ∞ and 0 are
      Γ-images of the level-1 horoball; Γ₀(𝔭) ⊂ Γ needs nothing new.
      Residue: same citation-level items as level 1 (PROOF.md §4 —
      Friedman Thm 3.8.1 covers cofinite Kleinian groups, so Γ₀(𝔭) with
      trivial character verbatim), plus the M𝔭 writeup itself.

## 7. Theorem G1𝔭 (guaranteed lower bound, two cusps)

Everything is as in `lower_bound_theory.md` (whose lemmas apply per
reference cell verbatim) with these changes.

**Reference constants** (on the 24-split reference mesh K_h, computed in
arb): γ, α_h^ref, τ (the level-1 slab constant τ_∞), S_M, S_Q, S_Σ, V_Σ
(per-face-optimized column heights H_F allowed in Lemma S — each column
is independent, so H may vary per floor face). Lemma G is re-verified on
the quarter-triangle floor faces (their diameters are ≤ the cell
diagonal, so the same sag bound applies; the center node carries the
cell lift).

**Lemma D (exact trace reproduction + disjoint copies).** The cusp
cross-sections at height Y are exactly tiled by boundary faces of the
mesh (top faces of the top-layer tets, per copy), and the CR
interpolant matches the mean of v on *every* mesh face. Hence, with
e = v − I_CR v,

    (D0)  t_∞(e) = t₀(e) = 0,   i.e.  t_α(I_CR v) = t_α(v)  exactly.

Consequently the trace functionals need no error terms at all: no slab
constants, no Young split of t(v)² (the boundary terms enter D_h with
coefficient β̃_α, uninflated), and

    (D3)  |a(e)| ≤ √6 α_h^ref √(Q_pc(e))     (six disjoint copies),

so σ_h = √6 α_h^ref and d_e has no trace part. ∎ (The level-1 files
carry a slab-based bound τ_h > 0 for the same quantity; that is valid
but strictly lossier — the certified level-1 run stands a fortiori.)

**Theorem G1𝔭 statement.** Fix a window W = [s⁻,s⁺] ∋ s₀, parameters
θ, θ₂, θ', α ∈ (0,1), ρ > 0, ν* = 1+ε₀ > 1. With λ⁺ = 1−(s⁻)²,
β_∞⁺ = 2(1−s⁻)/Y², β₀⁺ = β_∞⁺/5, Δκ = max_W |κ_c(s)−κ_c(s₀)|, define

    d₂ := ρ(1/θ₂ − 1) = ρ̃/θ₂,   ρ̃ := ρ(1−θ₂),   ω := 12 d₂ V_Σ,
    c_Q := 1 − (ω + ν*λ⁺) S_Q,        c_Σ := 1 − (ω + ν*λ⁺) S_Σ,
    λ̃  := ν*λ⁺(1 + S_M) + ω S_M,     β̃_α := ν*β_α⁺ + 4 d₂ Δκ²,
    σ_h := √6 α_h^ref                  (Lemma D0: no trace error),
    c_e := c_Q (1 − (1/θ₄−1)ρ_w) − ρ̃(1/θ − 1) σ_h²,
    d_e := λ̃(1+1/α) γ².

(The ℒ-shift error uses E² ≤ 2·(6V_Σ)M_Σ + 4Δκ²(t_∞²+t₀²), whence the
12 in ω and the 4 in β̃; a_Σ over the six slivers has ∫y⁻³ ≤ 6V_Σ.)

**Exact-main-term refinement.** The Young split v = Iv + e is pointwise,
so the *main* terms need no weight freezing — only the error steps
(E-M/E-t/E-a, which convert ‖e‖-type norms to ‖∇e‖) use sandwiched
weights. Define per tet the nonnegative stiffness remainder

    Q_rem := Σ_T ( ∫_T y⁻¹ − w_T^Q |T| ) ∇φ∇φᵀ  ⪰ 0   (∫_T y⁻¹ Taylor-enclosed),

the exact-weight mass M_ex (entries ∫_T φ_aφ_b y⁻³, Taylor-enclosed like
ℓ₀), and ρ_w := max_T (y_T⁺/y_T⁻ − 1). Then for any θ₄ ∈ (0,1):

    Q(v) ≥ [Q_pc + (1−θ₄) Q_rem](Iv) + (1 − (1/θ₄−1)ρ_w) Q_pc(e),
    M(v) ≤ (1+α) M_ex(Iv) + (1+1/α) γ² Q_pc(e),

(first line: Pythagoras on Q_pc plus the Young split of Q_rem(v) with
Q_rem(e) ≤ ρ_w Q_pc(e); second line: pointwise Young, exact weights on
the main term). Accordingly c_e is replaced by

    c_e := c_Q (1 − (1/θ₄−1) ρ_w) − ρ̃(1/θ − 1) σ_h².

Discrete data on the glued global CR space (exact 𝔽₅ scatter of the
reference element enclosures): Q_pc, Q_rem, M_ex, a, t_∞, t₀,
ℓ₀ = a + κ_c(s₀)(t_∞ + t₀);

    N_h := c_Q [Q_pc + (1−θ₄) Q_rem] + ρ̃(1−θ) ℓ₀ℓ₀ᵀ,
    D_h := λ̃(1+α) M_ex + β̃_∞ t_∞t_∞ᵀ + β̃₀ t₀t₀ᵀ,

with the boundary terms *uninflated*, by Lemma D0 (t_α(Iv) = t_α(v)
exactly, so no Young step is needed for the trace quadratics).

If c_Σ ≥ 0, c_e ≥ d_e and N_h − D_h ⪰ 0 on a family of windows covering
(0,1), then Γ₀(2+i)\ℍ³ has **no eigenvalue in (0,1)**.

*Proof.* Identical to Theorem G1 with: the multi-cusp criterion of §2 in
place of the level-1 criterion; the six per-copy slivers absorbed by the
reference Lemma S constants (per-copy maxima; disjoint columns); the
window shift of ℒ carried out with the combined error E above; the CR
projection step using (D1)–(D3) for the per-cusp trace and volume
functionals; Lemma R applied per vector (rank one positive for ℓ₀, two
rank-one negatives for t_∞, t₀) in the interval-to-midpoint reduction;
and the Rump certificate for the matrix inequality. Testing over all of
H¹ of the glued core (conforming across distinct-copy interfaces, free
on self-identifications) is a relaxation of the true form domain, as
before. ∎

## 8. The ladder: N𝔭 = 9 and 13 (plan, 2026-07-10)

Targets: 𝔭 = (3) (inert, N𝔭 = 9, index 10) and 𝔭 = (3+2i) (split,
N𝔭 = 13, index 14). Both prime level ⇒ two cusps; expected results:
"no eigenvalues in (0,1)" for each — new theorems, same machinery.

**What generalizes verbatim (prime-independent, proofs already written):**
Facts A–B hold for any prime 𝔭 with reps γ_k = (0,−1;1,k), k over
Gaussian-integer lifts of ℤ[i]/𝔭 — σ₀⁻¹γ_k = (1,k;0,1) is the same
computation, and the C2 fold-tiling argument of `0cuspchecks.md` only
needs residues closed under negation. Cusp classes: (0:1) alone at ∞,
N𝔭 cosets at 0. Per-cusp data: |T₀| = N𝔭/2, β₀ = 2(1−s)/(N𝔭 Y²),
t₀(1) = N𝔭/2, shared κ_c (equal Y). Nonzero-mode bound:
4π²Y²/N𝔭 > 1 ⟺ N𝔭 < 4π²Y² ≈ 61.7 — comfortable through the whole
near-term ladder (first bites at N𝔭 ≥ 61). Constant-direction scaling
*improves* with N𝔭: D(1)'s boundary part grows ~N𝔭/2 but ℒ(1)² grows
~N𝔭², so the ρ̃ demand shrinks. Lemma D0, exact main terms, Lemma I1,
Rump: untouched.

**What must change (small, localized):**
1. Residue arithmetic: ℤ[i]/(3) ≅ 𝔽₉ is NOT a prime field — replace the
   mod-5 integer tricks in `build_gluing` by generic ring ops on pairs
   (a+bi mod 3, i² = −1); for (3+2i), ℤ[i]/𝔭 ≅ 𝔽₁₃ with i ↦ 8
   (8² ≡ −1). Parametrize NP, the ℙ¹ point list, and the generator
   matrices mod 𝔭. The assert suite (permutations, round-trip, class
   sizes {1, N𝔭}, connectivity, t₀(1) = N𝔭/2, volume) transfers and is
   the error net.
2. Copy count NC = N𝔭 + 1 threaded through assembly (currently 6).

**Sizes and memory strategy (the only scaling cost is the Cholesky):**
- N𝔭 = 9: 10 copies × 6×3×3 reference (1296 tets) ≈ 26,700 dofs →
  dense A ≈ 5.7 GB: current machinery as-is.
- N𝔭 = 13: 14 copies: 6×3×3 → ≈ 37,400 dofs → 11.2 GB (borderline on a
  16 GB machine). Options in order: (a) leaner reference 6×3×2
  (≈ 25,000 dofs, 5 GB; Lemma D0 makes the coarser mesh viable — d_e is
  now only the γ² term, estimated ≈ 0.3 vs c_e ≈ 0.8); (b) memory-
  hygiene attempt at 11.2 GB; (c) out-of-core blocked right-looking
  Cholesky on a memmap — Rump's Thm 2.3 covers any standard-model
  execution order ("any library routine"), so a blocked own
  implementation is admissible; hours, not days.

**Execution order per level:** (i) M𝔭0 margin sweep (sparse LOBPCG,
cheap) — decision numbers: expect μ ≈ 3–3.5 at vol 3.05 (N=9) and
μ ≈ 2.5–3 at vol 4.27 (N=13) by sublinear-decay extrapolation
(7.3 → 4.4 at 6×); (ii) parameter search + certify. Watch-items: if μ
dips toward ~1.5, tighten ν* (1.01) and revisit window count; the S_M
and γ constants are reference-level and do not change with N𝔭.

**Status (worktree np9, 2026-07-11):** Residue ring generalized
(`GaussianResidueRing` / `set_level`). **N𝔭=9 CERTIFIED** (Theorem 3).
**N𝔭=13 CERTIFIED** (Theorem 4): mesh 8×3×2, ν*=1.001, 16 windows,
params (0.3, 0.98, 0.1, 0.85, 0.5); Rump with power-of-two diagonal
equilibration + per-row radii. M𝔭0 μ ≈ 1.45. Details: `PROOF.md`,
`PROGRESS_LADDER.md`.

**Beyond the ladder (noted, not planned):** composite/split levels
(e.g. (5) = (2+i)(2−i)) have > 2 cusps — Facts A–B need the general
cusp-width bookkeeping (still finitely many height-preserving
translation classes, but widths vary per cusp); the mode-bound ceiling
N𝔭 < 4π²Y² eventually forces more s-windows or a taller Y.
