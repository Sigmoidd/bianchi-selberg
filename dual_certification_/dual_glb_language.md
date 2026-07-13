# Dual language around GLB ≈ 6.45 — what is and is not certified

**Hard map (unchanged):**  
`eta_le_eta0=false`, `width_lt_tol=false`, `counting_certified=false`  
⇒ `pipeline_ok=true`, **`rung4_certified=false`**.

This note freezes honest dual-interval language so GLB numbers are not
over-read as a certified first eigenvalue.

---

## 1. Certified pieces (machine / theorem)

| Claim | Status | Source |
|-------|--------|--------|
| No Laplace eigenvalue in **(0,1)** on Γ₀(2+i) | **Certified** | independent_exclusion Theorem 2 / m3p |
| F5 gluing of 6 copies | **Exact** | congruence combinatorics |
| Lemma K tail (Fourier) | **Arb-enclosed** | lemma_K.py |
| Theorem D(K) defect formula under H,A,S | **Paper + Arb C1** | theorem_DK.tex; needs η≤η₀ |

FEM exclusion ⇒ continuum  
\[
\lambda_1(\Gamma_0(2+i)\backslash\mathbb{H}^3)\;\ge\;1.
\]
That is the only **hard** lower bound currently dual-safe.

---

## 2. Engineering discrete spectrum on \(V_h^{P1,\mathrm{per}}\)

Multi-copy CR with **cross-copy** PAIRINGS identified (conforming-periodic
flavor), self-faces free:

| Quantity | Value (mesh 4×2×2, 6 copies) | Flavor |
|----------|------------------------------:|--------|
| first positive discrete Neumann \(\lambda_h^+\) | **≈ 6.825** | **upper-flavor** for the *CR discrete* operator (Rayleigh) |
| CR GLB sketch \(\lambda_h^+/(1+\kappa_1^2 h_{\max}^2\lambda_h^+)\) | **≈ 6.451** | **lower-flavor** for continuum *if* full CR hypotheses hold on this space |
| Arb reference tet volumes (multi-copy mesh) | **384/384 GREEN** | `multi_copy_cr_arb_glb.py` |
| multi-copy κ₁ theory | **YELLOW** | not re-proved on glued space |
| \(J_h^{\mathrm{cross}}\) | **0** | by construction on \(V_h^{P1,\mathrm{per}}\) |

**Do not write:** “certified \(\lambda_1\in[6.45,6.83]\).”  
**May write:** “engineering CR bracket sketch \(\lambda_h^+\approx 6.83\), GLB sketch \(\approx 6.45\), under unverified multi-copy CR interpolation hypotheses; continuum lower bound still only \(\lambda_1\ge 1\) from FEM exclusion.”

---

## 3. What GLB is *not*

1. **Not** a dual upper bound (GLB is a lower bound when valid).  
2. **Not** automatic for Γ₀ multi-copy without replaying m3/m3p-style Arb checks (κ₁, h_T, weighted forms, BC).  
3. **Not** a substitute for Hejhal+D(K) residual certification (η≤η₀).  
4. **Not** a counting certificate on (1, λ₁).  

Free single-cell **Neumann without pairings** remains **lower-bound flavor only** for the quotient and must not be used as a dual upper bound.

---

## 4. Dual interval status (honest)

```
continuum lower (certified):     λ₁ ≥ 1
engineering CR lower sketch:     ≳ 6.45   (if CR theory applies)
engineering CR Rayleigh:         ≈ 6.83   (discrete upper flavor)
Hejhal+D(K) existence interval:  empty/huge until η ≲ η₀ ~ 6e-8
counting on (1, λ₁):             not certified (N candidate YELLOW)
```

Until η and counting move, **rung4_certified stays false**.

---

## 5. Allowed phrases for papers / certificates

| Allowed | Forbidden |
|---------|-----------|
| “FEM certifies exclusion of (0,1)” | “λ₁ = 6.45…” |
| “engineering first positive CR mode ≈ 6.83” | “certified dual interval [6.45, 6.83]” |
| “GLB sketch 6.45 under CR hypotheses” | “GLB proves continuum lower bound 6.45” |
| “pairing-conforming trial space \(V_h^{P1,\mathrm{per}}\)” | “Neumann free P1 upper-bounds the quotient” |
| “pipeline_ok; rung4_certified=false” | flipping hard flags without numbers |

---

## 6. Reproduce numbers

```bash
python cr_glb_p1_per.py --N1 4 --N2 2 --N3 2
python larger_poincare_residual.py --skip-poincare   # conforming block
python v_h_p1_per_spectrum.py --N1 4 --N2 2 --N3 2
```

---

*Maintained with dual_certification_ / morningbrief.md. Do not green hard map from this file alone.*
