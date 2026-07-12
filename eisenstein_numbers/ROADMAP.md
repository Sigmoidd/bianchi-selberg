# FEM independent exclusion — Eisenstein–Picard roadmap

```
Geometry → Reference cell → Float M0 → Interpolation constants
  → CR assembly → Window theorem → Interval certification → Congruence ladder
```

| Stage | Artifact | Status |
|-------|----------|--------|
| 1 Geometry | `GEOMETRY.md` + `geometry_fund.py` | **FROZEN** (EGM F_3, pairings, exact vol) |
| 2 Reference cell | `reference_cell.py` | **P_3 mesh** (default); R_comp legacy |
| 3 Float M0 | `fem_omega_m0.py` | **PASS** (Q1 on R_comp; margin healthy) |
| 4 Interpolation constants | `constants.py` | κ_sc, κ₁, β(s), **g1_coeffs** |
| 5 CR assembly | `cr_omega.py` | **PASS** on P_3; t(1)=\|T\|, vol≈vol(K_Y) |
| 6 Window theorem | `cr_omega.py` | **FLOAT PASS 8/8** (rho=55, N_tri=6,N3=3) |
| 7 Interval certification | `cert_omega.py` | **CERTIFIED 8/8** Arb+Rump (N_tri=6,N3=3,rho=55) |
| 8 Congruence ladder | `cert_omega_p.py` | **CERTIFIED Γ₀(N=3)** 8/8; N=7+ open |

**Run all:** `python -u run_roadmap.py`  
**Geometry:** `python -u geometry_fund.py`  
**G1 float:** `python -u cr_omega.py`  
**Interval cert:** `python -u cert_omega.py [N_tri N3]`  
**Γ₀ cert:** `python -u cert_omega_p.py 3 6 3`  
**Pairing matrices:** `python -u pairing_matrices.py`

**Level-1 group:** Γ = PSL(2, ℤ[ω]), K=ℚ(√−3).

**Certified:**
```
Theorem A (2026-07-11): λ₁(PSL(2,ℤ[ω])\ℍ³) ≥ 1  — cert_omega.py 6 3
Theorem B (2026-07-12): no eigenvalue in (0,1) on Γ₀(1−ω)\ℍ³ — cert_omega_p.py 3 6 3
Pairing matrices: pairing_matrices.py PASS
```

**Writeup:** [`PROOF.md`](PROOF.md) (Theorems A–B + G4 citations §7).  
**Pairings:** [`PAIRING_MATRICES.md`](PAIRING_MATRICES.md).  
**Progress log:** [`PROGRESS.md`](PROGRESS.md).

**Optionals done:** non-Neumann FE (`non_neumann_omega.py` 8/8, κ=I1+CZZ);
journal draft `papers/paper3_eisenstein.tex`.

**Remaining:** ladder float/cert N=7,13.

**Not equivalent to Thms 1–4** (ℤ[i] congruence FEM). This is the ω FEM path
(level 1 + first congruence rung).
