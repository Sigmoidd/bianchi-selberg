# Cusp-0 audit for Γ₀(2+i) (this file was found empty; filled 2026-07-10)

Every cusp-0 ingredient of the certified theorem, re-derived from scratch
and checked. Verdict up front: **all claims audit clean**; the two
genuinely subtle points are C2 (fold-compatibility of the chimney tiling)
and C6 (inherited precise invariance), both proved below. One convention
pitfall is flagged in C8.

## C1. Conjugated stabilizer at cusp 0 (lattice and fold)

σ₀ = S = (0,−1;1,0) ∈ Γ, cusp 0 = σ₀(∞). γ = (a,0;c,d) ∈ Γ₀(𝔭) fixes 0
(b = 0), ad = 1, c ∈ 𝔭. Direct computation:

    σ₀⁻¹ (a,0;c,d) σ₀ = (d, −c; 0, a).

So the conjugated stabilizer is upper triangular with unit diagonal
pairs (d,a), ad = 1: translations z ↦ z − c·a-unit with c ranging over
𝔭 ⇒ **translation lattice Λ₀ = 𝔭 = (2+i)ℤ[i]** ✓; the unit choices
a = ±i give the rotation z ↦ −z ⇒ **the ± fold is present at cusp 0** ✓.
Hence |T₀| = N𝔭/2 = 5/2 (folded cross-section area) ✓.

## C2. Fold-compatible tiling T₀ = ⊔ₖ (T + k)

Claim (Fact B): T₀ := ⊔_{k=0}^{4} (T + k) is a fundamental domain for
Λ₀ ⋊ ⟨−1⟩, where T is one for ℤ[i] ⋊ ⟨−1⟩. This needs the coset
decomposition ℤ[i]⋊± = ⊔ₖ (𝔭⋊±)·τ_k (τ_k = translation by k), which is
NOT automatic in a semidirect product. Proof: for (λ, ε) (z ↦ εz + λ),
(λ', ε)τ_k maps z ↦ ε(z+k) + λ' = εz + (εk + λ'), so (λ,ε) ∈ (𝔭⋊±)τ_k
⟺ λ − εk ∈ 𝔭 ⟺ k ≡ ελ (mod 𝔭). Since ε = ±1 and {0,…,4} represent
ℤ[i]/𝔭, exactly one k works for each (λ,ε) ✓. (The fold is compatible
because −1 acts on ℤ[i]/𝔭 and permutes the residues.) ∎

## C3. Chimney heights (Fact A/B computation)

σ₀⁻¹γ_k with γ_k = (0,−1;1,k): σ₀⁻¹ = (0,1;−1,0), and
(0,1;−1,0)(0,−1;1,k) = (1,k;0,1) — a *height-preserving* translation.
So in cusp-0 coordinates the five chimneys are (T+k)×(Y,∞), assembling
exactly to T₀×(Y,∞) at the same truncation height Y ✓. det γ_k = 1 ✓;
bottom rows (1:k) plus identity's (0:1) exhaust ℙ¹(𝔽₅) ✓.

## C4. Dual lattice / nonzero-mode bound at cusp 0

Dual of Λ₀ = 𝔭 under the pairing e^{2πi Re(μ̄z)}: Re(μ̄(2+i)m) ∈ ℤ for
all m ∈ ℤ[i] ⟺ μ̄(2+i) ∈ ℤ[i] (self-duality of ℤ[i]) ⟺
μ ∈ (2−i)⁻¹ℤ[i], shortest |μ| = 1/|2−i| = 1/√5 ✓. Mode bound:
(2π/√5)²Y² = 4π²Y²/5 ≈ 12.34 > λ⁺ = 1 ✓ (Lemma 1 applies per cusp;
folding only restricts to an invariant subspace, so the bound is
unaffected ✓).

## C5. Zero-mode defect and constraint bookkeeping (β₀, κ_c)

Level-1 Lemma 2 with general cross-section area |T_α|: the cusp-α
zero-mode contributes −|T_α|(1−s)c_α²/Y² with c_α the cross-section
mean; with t_α := ∫_{T_α} v dx = |T_α| c_α this is −(1−s)t_α²/(|T_α|Y²),
so β_∞ = 2(1−s)/Y² (|T_∞| = ½) and β₀ = (2/5)(1−s)/Y² (|T₀| = 5/2) —
matching the code's 2(1−s)/(N𝔭·Y²) ✓. Constraint: ∫_{C_α} u =
|T_α| c_α Y⁻²/(1+s) = t_α Y⁻²/(1+s), same κ_c(s) = 1/((1+s)Y²) for both
cusps because Y_∞ = Y₀ = Y ⇒ the constraint is rank-one in t_∞ + t₀ ✓.

## C6. Precise invariance of both collars at height Y

{y > Y}, Y ≥ 1, is precisely invariant under the FULL Γ = PSL(2,ℤ[i])
(level-1 Lemma 0: y(γP) ≤ 1/(|c|²y) ≤ 1/y). For cusp 0 and γ ∈ Γ₀(𝔭):
γ·σ₀H ∩ σ₀H ≠ ∅ ⟺ (σ₀⁻¹γσ₀)H ∩ H ≠ ∅ ⟹ σ₀⁻¹γσ₀ ∈ Γ_∞ (Γ-level
precise invariance; σ₀ ∈ Γ) ⟹ γ stabilizes 0 ✓. Nothing about Γ₀(𝔭)
is used beyond Γ₀(𝔭) ⊂ Γ. Same argument at ∞ with σ = id ✓.

## C7. Coset/cusp-class/gluing conventions in the code

- Class of coset Γ₀γ, γ = (a,b;c,d): copy γF's cusp end sits at
  γ(∞) = a/c (lowest terms since (a,c) is half of a unimodular pair);
  a/c ~_{Γ₀} ∞ ⟺ c ≡ 0 (mod 𝔭). Code tests c ≡ 0 on the bottom-row
  label (c:d) ✓. Class sizes 1 and 5 asserted ✓.
- Gluing: copy c's δ-source face glues to copy j with Γ₀γ_c δ⁻¹ = Γ₀γ_j,
  i.e. j = act((c:d), δ⁻¹) — code uses inv[name] with the involutions
  R, T_iR, S self-inverse and T1 inverted explicitly ✓; round-trip
  assert act(perm[n], δ) = n present ✓.
- t₀ assembly: the copy → collar coordinate change is the translation
  (1,k;0,1) (C3), measure- and height-preserving, so summing reference
  top-face integrals over the five class-0 copies computes ∫_{T₀} v dx
  exactly ✓.

## C8. Convention pitfall (flagged, not an error)

We do NOT normalize σ₀ to make Λ₀ unimodular (no det-scaling). All
formulas carry |T_α| explicitly, so this is consistent — but any future
comparison with literature that uses unit-covolume scaling matrices
(width factors in Eisenstein/scattering normalizations) must convert.
The mode bound, β₀, κ_c above are stated in the UNSCALED coordinates the
code uses.

## C9. Numerical cross-checks that would catch cusp-0 errors

- t₀(1) = 5/2 exactly (assembled) ✓ — wrong lattice/fold would break it.
- 1ᵀM_ex1 = vol_w(6 copies) ✓ against 6·vol_w(K_h).
- Gluing graph connected (1 component) ✓ — wrong permutation direction
  typically disconnects or double-glues.
- The certified window checks themselves: a wrong β₀ (e.g. missing the
  1/N𝔭) would move D(1) by ~2× and the parameter search / pencil
  diagnostics would have shown a very different profile.
