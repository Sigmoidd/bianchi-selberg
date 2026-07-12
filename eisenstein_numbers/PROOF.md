# Theorems (certified): Eisenstein–Picard level 1 and first congruence rung

**Date (machine):** level 1 — 2026-07-11; Γ₀(N=3) — 2026-07-12.  
**Method:** independent exclusion (Lax–Phillips + Theorem G1 / G1𝔭 + Arb + Rump BIT 46).  
**Not** the Selberg trace formula; **not** Theorems 1–4 of the ℤ[i] congruence ladder.

---

## Theorem A (level 1)

Let Γ = PSL(2, ℤ[ω]) with ω = e^{2πi/3}, acting on hyperbolic 3-space ℍ³.
Then the Laplace–Beltrami operator on L²(Γ\ℍ³) has **no eigenvalue in (0, 1)**.
Combined with the spectral decomposition for cofinite Kleinian groups
(Friedman, Thm 3.8.1), the spectrum in (0,1) is empty: **λ₁(Γ\ℍ³) ≥ 1**.

**Scope / relaxation.** Face pairings of the orbifold are **not** imposed
(Neumann free H¹(K)). As in `independent_exclusion/DESIGN.md`, testing over
the larger free space is a valid relaxation: positivity on free H¹(K)
implies the same on the subspace of functions that descend to the quotient.

**Reproduce:**

```text
python -u cert_omega.py 6 3
```

---

## Theorem B (congruence rung 1)

Let O = ℤ[ω], 𝔭 = (1−ω) (so N(𝔭) = 3), and

```text
Γ₀(𝔭) = { [a b; c d] ∈ PSL(2,O) : c ∈ 𝔭 }.
```

Then [Γ : Γ₀(𝔭)] = N(𝔭)+1 = 4, the quotient Γ₀(𝔭)\ℍ³ has exactly two
cusps (∞ and 0), and the Laplacian on L²(Γ₀(𝔭)\ℍ³) has **no eigenvalue in
(0, 1)**.

**Reproduce:**

```text
python -u cert_omega_p.py 3 6 3
```

---

## 1. Architecture (both theorems)

Three layers (same blueprint as Picard `independent_exclusion/PROOF.md`):

### 1.1 Reduction (analytic criterion)

**Level 1.** Lax–Phillips cusp analysis reduces λ ∈ (0,1) to positivity of
𝒜ₛ on {ℒₛ = 0} for the truncated core K_Y ⊂ F_3, every s ∈ (0,1)
(`DESIGN.md` §§0–3, with |T| = √3/6):

```text
|T| = area(P_3) = √3/6,
β(s) = (1−s)/(|T| Y²),   κ_c(s) = 1/((1+s) Y²).
```

**Congruence.** Two cusps ∞, 0 with cross-section areas
(`T0_AREA.md`, proved for field primes):

```text
|T_∞| = √3/6,   |T_0| = N(𝔭)·|T_∞| = √3/2   (N=3),
β_∞(s) = (1−s)/(|T_∞| Y²),   β_0(s) = (1−s)/(|T_0| Y²),
κ_c(s) = 1/((1+s) Y²)   (same Y for both cusps).
```

Continuous spectrum of cofinite Kleinian groups starts at 1
(Friedman 3.8.1); atoms in (0,1) are L² eigenvalues.

### 1.2 Guaranteed FEM lower bound

**Theorem G1** (level 1): Crouzeix–Raviart on an inner polyhedral mesh
K_h ⊂ K_Y; finite checks c_Σ ≥ 0, c_e > d_e, and N_h − D_h ⪰ 0 per
s-window. Constants: Lemma I1 (κ_sc), Lemma E (slab τ), Lemma S (floor
lift), Lemma G (floor inclusion). Full proofs:
`independent_exclusion/lower_bound_theory.md`.

**Theorem G1𝔭** (congruence): Reference-Cell Principle — NC = index
isometric copies of the P₃ cell; conforming glue on cross-copy faces;
Neumann on self-IDs. **Lemma D0:** σ_h = √(index)·α_h^ref (no τ in d_e)
when multi-copy tops stay free at y=Y. Coefficients:
`CERT_OMEGA_P.md` §3; Gaussian parallel `CONGRUENCE.md` §7.

### 1.3 Verified certificate

| Theorem | Script | Arithmetic |
|---------|--------|------------|
| A | `cert_omega.py` | Arb enclosures + Rump BIT 46 |
| B | `cert_omega_p.py` | Same stack; multi-copy scatter + two-cusp pencil |

Element Q, M, a, t by Taylor enclosures (`weighted_moments` from
`m3p_certify.py`); scalar window checks by rigorous arb comparison;
matrix PSD by Rump BIT 46 (2006) with power-of-two diagonal equilibration
and per-row radii (`rump_certify_inplace`).

---

## 2. Geometry (EGM fundamental domain)

**Group.** Γ = PSL(2, O_{−3}) = PSL(2, ℤ[ω]), disc D = −3, class number 1.

**Fundamental polyhedron** ([EGM98, §7.3]; [DP20, Thm 2.3 / §2.3]):

```text
F_3 = { (z,t) ∈ B_3 : z ∈ P_3 },
P_3 = T_up ∪ T_low   (explicit inequalities: geometry_fund.py / GEOMETRY.md),
K_Y = F_3 ∩ { t ≤ Y },   Y = 1.25.
```

**Volume.**

```text
vol(F)   = 3√3 · ζ_K(2) / (4π²)  ≈ 0.169156934402   (Humbert)
vol(K_Y) = vol(F) − |T|/(2 Y²)   ≈ 0.076780891331   @ Y=1.25
|T|      = √3/6.
```

**Computational mesh.** Tet mesh of K_Y over **P_3** (`build_P3_cell`).
Prism→tet splits order base vertices by index for conformity.

**Side-pairing generators** (Swan / EGM; matrices in `geometry_fund.GEN`
and `PAIRING_MATRICES.md`):

| δ | Matrix in SL(2,O) | Möbius |
|---|-------------------|--------|
| T₁ | `[[1, 1], [0, 1]]` | z ↦ z+1 |
| T_ω | `[[1, ω], [0, 1]]` | z ↦ z+ω |
| U | `[[ω², 0], [0, ω]]` | z ↦ ωz (order 3) |
| S | `[[0, −1], [1, 0]]` | inversion; floor self-paired |

Ideal-vertex cycles of **inverses** (polyhedron section; user ground truth,
audited by `check_generators.py`):

```text
T₁⁻¹ : (∞)(0 1)(ω ω²)
T_ω⁻¹: (∞)(0 ω)(1 ω²)
U⁻¹  : (∞)(0)(1 ω ω²)   # equals our forward U; same ⟨U⟩
S⁻¹  : (∞ 0)(1)(ω ω²)
```

---

## 3. Frozen certificate parameters

### 3.1 Theorem A (level 1)

```text
mesh:      P_3, N_tri=6, N3=3  →  648 tets, 1440 CR dofs
Y = 1.25,  ρ = 55,  θ = 0.7,  θ₂ = θ' = α = 0.5,  ν* = 1.05
s-grid:    8 uniform windows covering [0,1]
Taylor p = 5,  arb precision 128
κ:         κ_sc = √(1/π² + 1/15)   (Lemma I1)
```

Certified run (2026-07-11), abridged:

```text
Lemma G: 72 floor faces, admissible=True
1'Mm1 ≈ 0.076661  (≈ vol(K_Y)),  t(1) = √3/6
all 8 windows: c_e>d_e True, Rump PSD True
==> P_3 / PSL(2,Z[omega]) M3 CERTIFIED (all windows): True
```

### 3.2 Theorem B (Γ₀(1−ω))

```text
mesh:      same P_3 6×3 reference; NC = 4 copies
Y = 1.25,  θ = 0.5,  θ₂ = 0.9,  α = 0.2,  θ₄ = 0.5,  ρ̃ = 9,  ν* = 1.05
s-grid:    8 windows; Taylor p = 5; arb prec 128
σ_h:       √4 · α_h^ref   (Lemma D0; tops free)
n ≈ 5404 dofs
```

Certified run (2026-07-12), abridged:

```text
Lemma G arb: OK; t_∞(1)=|T_∞|, t_0(1)=|T_0|
min c_e/d_e ≳ 11; Rump PSD 8/8
==> Γ₀(1−ω) CERTIFIED (all windows): True
```

Float precursor: `congruence_omega_proto.py 3 6 3` (min μ ≈ 1.12 > 0).

---

## 4. Congruence ingredients (Theorem B)

| Ingredient | Artifact | Status |
|------------|----------|--------|
| Cosets ≅ ℙ¹(O/𝔭) | `residue_omega.py` | exact; N=3 has 4 points |
| Gluing perms (right action by δ⁻¹) | same | S²=id, U³=id, round-trips |
| Cusp areas \|T_∞\|, \|T_0\| | `T0_AREA.md`, `cusps_omega.py` | proved |
| Face dictionary on P₃ | `face_pairings_p3.py` | RIGHT↔LEFT, LOW↔UP, U, S |
| Pairing matrices g∈Γ | `PAIRING_MATRICES.md` | frozen SL(2,O) + checks |
| Lemma D0 (tops free) | face dict excludes TOP | verified in cert |

**N=3 specialities.** ω ≡ 1 in 𝔽₃ ⇒ T₁ ≡ T_ω as residue translations;
U ≡ id on ℙ¹(𝔽₃). Cross-copy glue uses T₁/T_ω and S; geometric U-maps
remain internal / residue-trivial.

---

## 5. Why each check suffices (summary)

An eigenfunction u with eigenvalue λ = 1−s² ∈ (0,1) is L² and
real-analytic. Cusp zero modes are exactly of the form c(y/Y)^{1−s}
(per cusp). The Lax–Phillips reduction turns Q(u)−λM(u)=0 into
𝒜ₛ(u|_K) ≤ 0 with ℒₛ(u|_K)=0, while u|_K ≠ 0 by unique continuation.
The pencil inequality, verified for all s through the finite window
checks of Theorem G1 / G1𝔭, gives 𝒜ₛ > 0 on {ℒₛ = 0} — contradiction.
The spectral decomposition upgrades “no eigenvalues in (0,1)” to
“no spectrum in (0,1)”.

---

## 6. Relation to other results in this repo

| Result | Group | Status |
|--------|--------|--------|
| Picard level-1 λ₁ ≥ 1 | PSL(2,ℤ[i]) | `independent_exclusion/PROOF.md` |
| Thms 1–4 congruence | Γ₀(𝔭) ≤ PSL(2,ℤ[i]), N𝔭=5,9,13 | `independent_exclusion/` |
| **Theorem A** | **PSL(2,ℤ[ω]) level 1** | **this file / `cert_omega.py`** |
| **Theorem B** | **Γ₀(1−ω) ≤ PSL(2,ℤ[ω])** | **this file / `cert_omega_p.py`** |
| STF / CE gate B&lt;1 for ω | same group (trace-formula path) | `framework.py`, `AUDIT.md` |
| Ladder rungs N=7,13 | Γ₀(π) over ℤ[ω] | combinatorics only |

The STF/CE path and the FEM path are **independent** checks of the same
level-1 spectral bound. Theorem B has no STF counterpart in this repo.

---

## 7. Citations required in a formal paper writeup (G4)

These are classical and not re-proved numerically:

1. **Spectral decomposition / continuous spectrum ≥ 1** — Friedman,
   *The Selberg trace formula for PSL(2,O_K)*, arXiv:math/0612807,
   Theorem 3.8.1 (cf. [EGM98 §6.2]). Used for the λ₁ ≥ 1 phrasing;
   “no L² eigenvalues in (0,1)” is the machine-certified core.
2. **EGM fundamental domain F_d** — Elstrodt–Grunewald–Mennicke,
   *Groups Acting on Hyperbolic Space*, Springer 1998, §7.2–7.3;
   P_3 formula [DP20, §2.3] (or equivalent Swan/EGM presentation).
3. **Shimizu horoball / precisely invariant cusp** — DESIGN.md Lemma 0
   (Shimizu, Ann. of Math. 77 (1963); Leutbecher, Math. Z. 100 (1967));
   same argument with N(c) ≥ 1 for O = ℤ[ω].
4. **Real-analyticity of eigenfunctions / unique continuation**
   (Morrey–Nirenberg; elliptic regularity).
5. **Rump, BIT Numer. Math. 46 (2006)** — Thm 2.3 / Cor 2.4 / Lemma 2.5
   (verification of positive definiteness).
6. **CR interpolation constant** — Lemma I1 (Payne–Weinberger /
   Bebendorf) or CZZ20/CP22 for the sharper κ.
7. **Side-pairing generators of Γ** — Swan; EGM Ch. 7; matrices and
   mesh dictionary in `PAIRING_MATRICES.md` / `GENERATOR_CYCLES.md`.

Software correctness-critical stack: python-flint (Arb), IEEE-754 double
for Rump’s floating Cholesky (as permitted by Rump: any LAPACK dpotrf).

---

## 8. Optional: non-Neumann FE space and κ=CZZ

The certified level-1 space uses Neumann free H¹(K) (relaxation).
Imposing the pairing dictionary on CR face dofs yields the geometric
orbifold FE space (`non_neumann_omega.py`):

```text
python -u non_neumann_omega.py 6 3
# raw 1440 → paired 1342 dofs; t(1)=|T|; ker(Q)≃const
# G1 float 8/8 for κ_sc (Lemma I1) and κ₁ (CZZ)
# min c_e/d_e ≈ 1.53 (I1), ≈ 2.09 (CZZ)
```

Positivity on free H¹ already implies positivity on this subspace;
the optional run is a geometric check + sharper-κ probe, not a second
interval cert.

**Journal packaging:** `papers/paper3_eisenstein.tex` (Paper III outline:
Theorems A–B, geometry, pairings, criterion, G1/G1𝔭, certificates).

## 9. What is *not* claimed

- No claim for Γ₀(𝔭) at N(𝔭) ∈ {7,13} (combinatorics only so far).
- No re-proof of Theorems 1–4 on ℤ[i].
- Congruence multi-copy glues cross-copy faces via the frozen dictionary
  + residue action; self-IDs remain Neumann.
- The STF constant B(Q(ω)) is a separate software gate; these certs do
  not use the trace formula or conjugacy-class inventory.

---

## 10. Commands (reproduce both theorems)

```text
# Geometry freeze
python -u geometry_fund.py
python -u check_generators.py
python -u pairing_matrices.py
python -u non_neumann_omega.py 6 3   # optional paired FE

# Level 1
python -u cr_omega.py
python -u cert_omega.py 6 3

# Congruence rung 1
python -u residue_omega.py
python -u cusps_omega.py
python -u face_pairings_p3.py
python -u congruence_omega_proto.py 3 6 3
python -u cert_omega_p.py 3 6 3
```
