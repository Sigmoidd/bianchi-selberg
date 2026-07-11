# Brief: Bianchi congruence exclusion & ℙ¹(ℤ[i]/𝔫)

**Audience:** collaborators · **Worktree:** `np9` · **Date:** 2026-07-11  
**Method:** verified FEM (Crouzeix–Raviart + Rump PSD), **no** Selberg trace formula.

---

## Certified theorems (full statements)

All for the Laplacian on L²(Γ\ℍ³), trivial character; **no eigenvalue in (0,1)**.

| Thm | Group | Level | Index | Status |
|-----|--------|-------|-------|--------|
| **1** | PSL(2,ℤ[i]) | 1 | 1 | Certified (prior) |
| **2** | Γ₀(𝔭), 𝔭=(2+i) | N𝔭=5 | 6 | Certified |
| **3** | Γ₀(𝔭), 𝔭=(3) inert | N𝔭=9 | 10 | Certified (np9) |
| **4** | Γ₀(𝔭), 𝔭=(3+2i) split | N𝔭=13 | 14 | Certified (np9) |

**Theorem 3.** For 𝔭=(3) ⊂ ℤ[i], the Laplacian on L²(Γ₀(𝔭)\ℍ³) has no eigenvalue in (0,1).  
*Certificate:* mesh 6×3×3, 26 532 CR dofs; Y=1.25; ν\*=1.02; (θ,θ₂,α,θ₄,ρ̃)=(0.5,0.9,0.15,0.8,1.5); 8 s-windows; Rump PSD all windows. Float margin min μ≈2.40.

**Theorem 4.** For 𝔭=(3+2i) ⊂ ℤ[i], the Laplacian on L²(Γ₀(𝔭)\ℍ³) has no eigenvalue in (0,1).  
*Certificate:* mesh 8×3×2, 33 176 dofs; ν\*=1.001; (θ,θ₂,α,θ₄,ρ̃)=(0.3,0.98,0.1,0.85,0.5); **16** s-windows; Rump with **power-of-two diagonal equilibration** SAS (s_i=2^{round(−½log₂ a_ii)}) and **per-row** radius shifts. Float margin min μ≈1.45.  
*Note:* Unscaled uniform max-row radius exceeded λ_min(A) while A was still SPD — a false negative, not a spectral obstruction.

**Common architecture (Thms 2–4).** Multi-cusp Lax–Phillips criterion (two cusps at prime level) → Theorem G1𝔭 (CR lower bound, Lemma D0: exact cusp traces, σ_h=√(index)·α_h^ref) → arb enclosures + Rump BIT 46. Cosets ↔ ℙ¹(ℤ[i]/𝔭); right action (c:d)·M=(c,d)M drives face gluing. Master writeup: `PROOF.md`.

---

## Combinatorics: right action on ℙ¹(ℤ[i]/𝔫)

| Layer | Status |
|--------|--------|
| Action formula (c:d)·M | **General** (any R) |
| Field ℙ¹ (prime 𝔭, 𝔽₉) | **Done** — ladder |
| Square-free CRT ℙ¹ ≅ ∏ ℙ¹(R_j) | **Done** — e.g. 𝔫=(5)=(2+i)(2−i), \|ℙ¹\|=36, gluing asserts green |
| Arbitrary square-free 𝔫 | **Near** — structure exists; factor 𝔫 → product of local fields |
| 𝔭^k (local non-field rings) | **Open** — unimodular-row ℙ¹ |
| Multi-cusp geometry + G1𝔫 + certify | **Open** — >2 cusps, widths, RAM |

**Index (square-free 𝔫):** N(𝔫)∏_{𝔭\|𝔫}(1+1/N𝔭)=∏(N𝔭+1)=\|ℙ¹\|.  
**Code:** `congruence_prototype.py`, note `P1_ACTION.md`. Prime regression for Thms 2–4 preserved.

---

## Takeaways for collaborators

1. **Three congruence rungs certified** (N𝔭=5,9,13) plus level 1 — all in the verified-FEM pipeline.  
2. **N=13 was knife-edge on Rump conditioning**, not on the exclusion form (μ≈1.45>0); equilibration + per-row radii close it on a 16 GB machine.  
3. **ℙ¹ right action is not the bottleneck for “general 𝔫”**; multi-cusp analysis and certificates are.  
4. **Main repo isolation:** this work lives in worktree `np9`; merge to main is a separate step.

**Reproduce (from `independent_exclusion/`):**  
Thm 3: `certify(6,3,3,params=(0.5,0.9,0.15,0.8,1.5),level='(3)')`  
Thm 4: `NU_STAR=1.001; NWIN=16; certify(8,3,2,params=(0.3,0.98,0.1,0.85,0.5),level='(3+2i)')`  
Gluing: `python -u congruence_prototype.py gluing`
