# Methodology review request — certified Selberg trace-formula bound on λ₁ for Bianchi groups

You are reviewing the **mathematical methodology** of a computational project (do
not edit code; this is a read-only critical review). Repo: `bianchi-selberg/`.
Files of interest: `picard_stf.py` (engine), `verify_group_data.py`,
`verify_matthies.py`, `verify_eisenstein.py`, `README.md`.

## Goal
Rigorously (interval-arithmetic certified, via python-flint/Arb) exclude
exceptional Laplace eigenvalues 0 < λ < 1 on Γ\ℍ³ for Bianchi groups
Γ = PSL(2,O_K), K imaginary quadratic — the archimedean analogue of Selberg's
conjecture. Best known unconditional general bound is λ₁ ≥ 975/1024
(Blomer–Brumley 2011). We attack specific small groups by direct computation.

## Core method (level 1, Γ = PSL(2,ℤ[i]) done)

Spectral parametrization: λ = 1 + r², so λ<1 ⇔ r = iσ, σ∈(0,1] (complementary
series). λ₀=0 ⇔ r₀=i.

Selberg trace formula for cofinite Kleinian groups (Friedman arXiv:math/0612807
Thm 4.1.1; Balkanova et al arXiv:1712.00880 Thm 2.2):
    Σ_j h(r_j) = I + NCE + CE + PAR + (scattering terms),
h even, holomorphic in a strip; g = Fourier transform of h.

**Test function:** h(r) = sinc^{2k}(δr) = (sin δr/δr)^{2k} ≥ 0 on ℝ, with
g = ĥ supported in [−2kδ, 2kδ] (a cardinal B-spline). Choose 2kδ ≤ ℓ₀ = systole
= log((3+√5)/2). Then **every loxodromic term vanishes** (g(logN(T))=0 since
logN(T) ≥ ℓ₀ > support). This is the key trick; it also works at every
congruence level since lengths only grow.

**Positivity criterion.** Define B := (geometric side) − h(i). On the spectral
side, λ₀=0 gives h(i); every exceptional eigenvalue λ=1−σ² gives
h(iσ) = (sinh δσ/δσ)^{2k} > 1; every other discrete eigenvalue gives h(r_j) ≥ 0
(h ≥ 0). Hence
    Σ_{exceptional} h(iσ_j) = B − Σ_{tempered} h(r_j) ≤ B.
If **B < 1** (certified), no exceptional eigenvalue can exist (each would
contribute > 1), so **λ₁ ≥ 1**. Result: B ≤ 0.319954 < 1 for ℤ[i], k=2,
δ=0.240365.

**Claim to scrutinize:** is "B<1 ⇒ no exceptional eigenvalue" airtight? In
particular: (a) are ALL non-λ₀ discrete spectral terms h(r_j) ≥ 0 (including any
residual/small eigenvalues in (0,1) other than complementary series, and the
possibility r_j on the imaginary axis)? (b) Is the continuous-spectrum /
scattering contribution correctly moved to the geometric side with the right
sign, so that B genuinely upper-bounds the exceptional sum? (c) Is dropping the
tempered terms (h(r_j) ≥ 0) valid — could any tempered term be negative because
sinc^{2k} with even power is ≥ 0, yes, but verify the exponent parity argument.

## Terms (ℤ[i]), all Arb-certified
- I = vol/(4π²)∫h(r)r²dr, vol = 2ζ_K(2)/π² = 0.30532 (Humbert; matches Then).
- NCE (non-cuspidal elliptic) = g(0)·log(7+4√3)/9  [one order-3 class,
  N(T₀)=(2+√3)², |E(R)|=3 (PSL); R ~ R⁻¹ so it is a single conjugacy class,
  sin²(π/3)=3/4, giving 1/(4·3·(3/4)) = 1/9].
- CE (cuspidal elliptic) = (5/16 log2) g(0) + (1/4)∫₀^∞ g tanh(x/2) dx
  [4 order-2 classes].
- PAR/scattering: uses η(ℤ[i]) = 2γ+2log2+3logπ−4logΓ(¼), and
  φ(s)=π ζ_K(s−1)/((s−1)ζ_K(s)); the φ′/φ integral is done by contour shift to
  Re s=2 through the entire ξ~(s)=s(s−1)π^{−s}Γ(s)ζ_K(s); only prime term n=2
  survives (log 3 > support).
- Infinite tails all bounded by explicit closed forms; B-spline integrals split
  at knots so integrands are polynomials (acb_calc needs holomorphic integrands).

## Independent validation (the part I'm most confident in)
`verify_matthies.py`: apply the trace formula to the counting-limit test function
h_k = 1_{|r|<k} (so g(0)=k/π, h(0)=1). Each term yields a k-asymptotic; we
reconstruct Matthies' independently-derived Weyl-law coefficients (reported in
Aurich–Steiner–Then, gr-qc/0404020):
    a₂ = −3/(2π),
    a₃ = (1/π)[13/16 log2 + 7/4 logπ − logΓ(¼) + 2/9 log(2+√3) + 3/2].
All of a₂ and every a₃ constant match EXACTLY (Arb). The single 2/9 log(2+√3)
proves the non-cuspidal elliptic inventory is complete (one order-3 class).

**Scrutinize:** is matching the Weyl-law a₂,a₃ a *sound* proof that the geometric
side is correct, or could errors cancel? Which terms does a₃ actually constrain,
and which remain unconstrained by this check (e.g. does the identity/vol term or
the h(i) subtraction get tested)? Is the counting-limit derivation of each term's
k-coefficient rigorous (I did the sinc→indicator limit informally)?

## Extension to ℤ[ω] (Eisenstein–Picard), level 1 — in progress
Certified mechanical constants (`verify_eisenstein.py`): vol=0.16916,
systole ℓ₀=0.86255, L(1,χ₋₃)=π/(3√3), L′(1,χ₋₃)=0.222663,
η(ℤ[ω])=(V/π)·6(γL(1,χ₋₃)+L′(1,χ₋₃))=0.94550 [V=√3/2], [Γ∞:Γ′∞]=3.
Derived identity: η_Λ=(V/π)c₀ (c₀=const term of Σ|μ|⁻²ˢ at s=1; no −w term
because T(1⁻)=0), validated on ℤ[i].

**Gated:** the O₃ non-cuspidal + cuspidal elliptic inventory and the scattering
normalization for K=ℚ(√−3). Definitive source EGM book Ch.7 (not online);
alternative is direct computation validated against the ℤ[i] elliptic value.

## Specific questions for your review
1. Soundness of "B<1 ⇒ λ₁≥1": any spectral-side term that could be negative or
   any exceptional contribution < the h(iσ)>1 bound I claim?
2. Is the Matthies a₂/a₃ reconstruction a genuine independent check of the
   geometric side, and precisely which terms does it leave untested?
3. The η_Λ=(V/π)c₀ derivation — is the partial-summation argument correct
   (boundary term handling)?
4. For ℤ[ω]: is validating a direct elliptic-term computation against the ℤ[i]
   value 2/9 log(2+√3) a sound way to trust it for the new field, or is that a
   non-sequitur (different group)?
5. Any hidden normalization/convention pitfalls (factors of 2, |D_K|^{...} in the
   scattering determinant, PSL vs SL, orientation) you'd flag before trusting a
   ℤ[ω] number.

Give a critical, specific review. Flag anything unsound. Do not rubber-stamp.
