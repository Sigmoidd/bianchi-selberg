# Stage 1 — Geometry freeze (PSL(2, ℤ[ω]) level 1)

**Status: FROZEN** for level-1 FEM input.  
Runnable check: `python -u eisenstein_numbers/geometry_fund.py`.

---

## 1. Group and cusp

| Object | Value |
|--------|--------|
| Γ | PSL(2, ℤ[ω]), ω = e^{2πi/3} = −1/2 + i√3/2 |
| O_K | ℤ[ω], K = ℚ(√−3), disc D = −3, class number 1 |
| Cusp | single cusp ∞ (class number 1) |
| Γ_∞ | Λ ⋊ U, Λ ≅ ℤ[ω] (translations), U = ⟨z ↦ ωz⟩ order 3 |
| [Γ_∞ : Γ'_∞] | 3 (Friedman Cor. 5.3.3; `verify_eisenstein.py`) |
| covol(Λ) | V_Λ = √3/2 |
| **\|T\|** (orbifold cross-section) | V_Λ / 3 = **√3/6** |

---

## 2. Fundamental polyhedron F_3 (EGM)

**Primary citation.** Elstrodt–Grunewald–Mennicke, *Groups Acting on Hyperbolic Space*, Springer 1998, **§7.2–7.3**: construction of  
B_d (Ford region), P_d (planar domain for Γ_∞ on ℂ), and  

```text
F_d = { (z,t) ∈ B_d : z ∈ P_d }.
```

**Theorem** ([EGM98, §7.3]; restated as [DP20, Thm 2.3]): F_d is a fundamental domain for Γ_d = PSL(2, O_d) on ℍ³.

**Explicit P_3** ([DP20, §2.3], matching EGM for d = 3):

```text
P_3 = T_up ∪ T_low,

T_up  = { x+iy | 0 ≤ x,   x/√3 ≤ y ≤ (1−x)/√3 }     (⇒ x ≤ 1/2),
T_low = { x+iy | 0 ≤ x ≤ 1/2,  −x/√3 ≤ y ≤ x/√3 }.
```

Euclidean area(P_3) = **1/(2√3) = √3/6 = |T|**.

**B_3 for d = 3.** The only exterior Ford sphere cutting the prism P_3 × ℝ_{>0} is the unit hemisphere |z|² + t² = 1. Moreover max_{P_3}|z| = 1/√3 < 1, so

```text
F_3 = { (z,t) : z ∈ P_3,  t ≥ √(1−|z|²) }.
```

Code: `geometry_fund.in_P3`, `in_B3`, `in_F3`.

---

## 3. Truncated core and reference cell

### 3.1 Canonical core (EGM truncation)

```text
K_Y = F_3 ∩ { t ≤ Y },   Y > 1 fixed (default Y = 1.25).
```

Faces: top P_3×{Y}, floor (unit hemisphere over P_3), vertical walls over ∂P_3.

### 3.2 Computational reference cell (FEM mesh)

**Default (CR / G1 / paper path):** EGM section P_3 — `build_P3_cell` / `build_reference_cell(..., domain="P3")`.

```text
K_Y = { (z,y) : z ∈ P_3,  y_f(z) ≤ y ≤ Y },  y_f = √(1−|z|²)
+ optional Lemma-G floor lift (δ̄ ~ O(h²)).
```

| Property | P_3 / K_Y (default mesh) | R_comp (legacy Q1 M0 only) |
|----------|--------------------------|----------------------------|
| planar area \|T\| | √3/6 | √3/6 |
| max \|z\| | 1/√3 ≈ 0.577 | 1.0 (degenerate y_f→0 at z=1) |
| min y_f | √(2/3) ≈ 0.816 | 0 |
| hyp. volume | vol(F)−\|T\|/(2Y²) ≈ 0.0768 | larger (~0.18); not EGM |
| CR / G1 | **float PASS 8/8** | do not use (S_Q blows up) |

Prism→tet extrusion sorts base vertices by index so shared faces get a
consistent diagonal (conforming mesh; required for ker(Q)=constants).

---

## 4. Side-pairing maps (replace Neumann)

Generators ([Swan71]; EGM Ch. 7), implemented in `geometry_fund.GEN`:

| Generator | Matrix | Action on ℍ³ |
|-----------|--------|----------------|
| T₁ | `[[1,1],[0,1]]` | z ↦ z+1 |
| T_ω | `[[1,ω],[0,1]]` | z ↦ z+ω |
| U | `[[ω²,0],[0,ω]]` | z ↦ ωz (order 3) |
| S | `[[0,−1],[1,0]]` | inversion; pairs floor |

**Face map for F_3**

| Face | Pairing |
|------|---------|
| floor \|z|²+t²=1 | **S** (self-paired) |
| vertical walls over ∂P_3 | elements of ⟨T₁, T_ω, U⟩ = Γ_∞ |
| top y=Y | truncation (Lax–Phillips); matched by t-functional |

**Operational maps on the P_3 mesh** (`face_pairings_p3.py`)

| Generator | Edges of ∂P_3 paired |
|-----------|----------------------|
| T₁ | RIGHT (x=½) ↔ LEFT (x=0) |
| T_ω | LOW ↔ UP (slanted sides) |
| U | vertical faces via z↦ωz |
| S | FLOOR ↔ FLOOR (sphere inversion) |
| top y=Y | free (cusp ODE / β) |

Used by multi-copy gluing (`congruence_omega_proto.py`).  
Checks: `geometry_fund.verify_side_pairings()`; `python -u face_pairings_p3.py`.

**Neumann relaxation (level-1 cert).** Side/floor IDs are **not** imposed at level 1 (Picard `DESIGN.md`). Congruence multi-copy **does** use the dictionary above for cross-copy faces.

---

## 5. Exact reference-cell volume

### 5.1 EGM truncated core K_Y (paper domain)

**Humbert formula** for vol(Γ\ℍ³):

```text
vol(F) = |D|^{3/2} ζ_K(2) / (4 π²)
       = 3√3 · ζ_K(2) / (4 π²),   D = −3,
ζ_K(2) = ζ(2) L(2, χ_{−3}).
```

Numeric (arb): **vol(F) ≈ 0.169156934402** (matches `verify_eisenstein.py`).

**Cusp tail** (product structure of F ∩ {t > Y} ≅ P_3 × (Y,∞), hyp. measure t^{−3} dt dx dy):

```text
vol(F ∩ {t > Y}) = |T| / (2 Y²),   |T| = √3/6.
```

**Closed formula:**

```text
vol(K_Y) = vol(F) − |T| / (2 Y²).
```

At default Y = 1.25:

| quantity | value |
|----------|-------|
| vol(F) | 0.169156934402 |
| tail | 0.092376043070 |
| **vol(K_Y)** | **0.076780891331** |

Cross-check: midpoint quadrature of ∫_{P_3} ½(yf⁻² − Y⁻²) dx dy reproduces the closed form to relative error ≲ 10⁻³ at n=600.

### 5.2 Mesh volume self-check

On the default **P_3** mesh, `1ᵀ M 1` must match `vol(K_Y)` (up to lift O(δ̄)).  
Observed: 1ᵀM1 ≈ 0.07666 vs exact 0.07678 at N_tri=6, N3=3.

Legacy R_comp (M0 only): `vol_K_comp_quad(Y)` = ∫_{R_comp} ½(yf⁻²−Y⁻²).

API: `vol_F_exact()`, `vol_KY_exact(Y)`, `cusp_tail_volume(Y)`, `vol_K_comp_quad(Y)`.

---

## 6. Load-bearing inequalities (unchanged)

| Fact | Status |
|------|--------|
| Shimizu: y(γP) ≤ 1/(\|c\|² y) for c≠0 | min N(c)=1 for 0≠c∈ℤ[ω] ✓ |
| Nonzero modes ≥ 0 | dual Λ^* min length 2/√3; 4π²Y²/\|μ\|² ≫ 1 at Y=1.25 ✓ |
| Zero-mode ODE | β = (1−s)/(\|T\| Y²) ✓ |
| Spectral atoms on (0,1) | Friedman Thm 3.8.1 ✓ |

---

## 7. Checklist (geometry tasks)

- [x] Cite EGM/DP20 for complete fundamental polyhedron F_3
- [x] Identify computational reference cell vs EGM K_Y
- [x] Explicit side-pairing generators + face table (for future non-Neumann)
- [x] Exact vol(K_Y) formula; wire into M0 volume self-check
- [x] Mesh true P_3 (not R_comp) for CR/G1; conforming prism split
- [x] Float G1: all 8 windows PASS at rho=55, N_tri=6, N3=3

**Interval cert (Stage 7):** `python -u cert_omega.py` — **PASS 8/8** at
N_tri=6, N3=3, rho=55 (Lemma G arb + exact-weight Taylor enclosures +
Rump SAS/per-row). Claims λ₁ ≥ 1 under Neumann relaxation (DESIGN.md).
