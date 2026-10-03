# bianchi-selberg — Repo Structure Guide

## TL;DR: Four Certified Proofs that λ₁ ≥ 1

| Method | Path | DOFs | Result |
|--------|------|------|--------|
| **Selberg Trace Formula** | `picard_stf.py`, `bianchi_omega_arb.py` | — | B ∈ [0.305, 0.320] < 1 ✅ |
| **Weyl Law Reconstruction** | `verify_matthies.py` | — | All a₂, a₃ coefficients match exactly ✅ |
| **FEM Level 1** | `independent_exclusion/m3_certify.py` | 5,544 CR | λ₁ ≥ 1 (no eigenvalues in (0,1)) ✅ |
| **FEM Congruence N𝔭=5** | `independent_exclusion/m3p_certify.py` | 28,400 CR | Γ₀((2+i)): no eigenvalues in (0,1) ✅ |
| **FEM Congruence N𝔭=9** | `independent_exclusion/m3p_certify.py` | 26,532 CR | Γ₀((3)): no eigenvalues in (0,1) ✅ |
| **FEM Congruence N𝔭=13** | `independent_exclusion/m3p_certify.py` | 33,176 CR | Γ₀((3+2i)): no eigenvalues in (0,1) ✅ |

---

## Directory Map

### 📊 Level 1: Picard Group PSL(2, ℤ[i]) via Trace Formula

**Goal:** Prove λ₁ ≥ 1 on the Picard orbifold using the Selberg trace formula.

| File | Role |
|------|------|
| `picard_stf.py` | **Main engine.** All trace-formula terms (identity, elliptic, cuspidal, scattering) with Arb-certified bounds. |
| `verify_group_data.py` | Independent verification: systole, centralizers, Friedman identity, lattice Euler constant η. All PASS. |
| `verify_matthies.py` | **KEY VALIDATION.** Reconstructs Matthies' Weyl-law coefficients a₂ and a₃ from engine constants; all match exactly (Arb-certified). |
| `elliptic_inventory.py` | Direct O_K arithmetic for non-cuspidal elliptic terms; validates 1/9·log(2+√3). |
| `cuspidal_ce.py` | Cuspidal elliptic enumeration + Friedman identity check. |
| `bianchi_omega.py` | High-precision mpmath assembly (B ≈ 0.3105 Picard, B ≈ 0.534 Eisenstein). |
| `bianchi_omega_arb.py` | **Certified version.** B(ℤ[i]) ∈ [0.3046, 0.3179], B(ℤ[ω]) ∈ [0.5252, 0.5437], both < 1. |
| `scan.py` | Robustness check: B(k, δ) over parameter space. |
| `final_run.py` | Headline certified run. |

**Status:** ✅ **DONE.** B < 1 certified. Picard level 1 proven and cross-validated.

---

### 🌐 Eisenstein–Picard Group PSL(2, ℤ[ω]) via Trace Formula

**Goal:** Same for the Eisenstein–Picard orbifold.

| File | Role |
|------|------|
| `verify_eisenstein.py` | Mechanical constants (vol, systole, L(1,χ₋₃), L′(1,χ₋₃), η(ℤ[ω])). All Arb-certified. |
| `bianchi_omega_arb.py` | Eisenstein trace-formula assembly: B(ℤ[ω]) ∈ [0.525, 0.544] < 1. |
| `eisenstein_numbers/` | Independent FEM proof (see below). |

**Status:** ✅ **DONE.** B < 1 certified (trace formula), independent FEM cert also done.

---

### 🏗️ Independent Exclusion: FEM + Rump Verification

**Path:** `independent_exclusion/`

**Goal:** Prove λ₁ ≥ 1 without the Selberg trace formula. Uses Lax–Phillips cusp reduction + Crouzeix–Raviart FEM + rigorous Rump matrix verification.

| File | Role |
|------|------|
| **`PROOF.md`** | Master document: Theorems 1–4, architecture, frozen certificate parameters. |
| **`DESIGN.md`** | Lax–Phillips reduction to pencil inequality on truncated domain (Lemma 0–3). |
| **`lower_bound_theory.md`** | Theorem G1: guaranteed CR lower bound via three finite checks (scalar inequalities + PSD matrix). Lemma I1 (CR constant), Lemma S (sliver), Lemma E (trace functional). |
| **`m3_certify.py`** | **Level 1 certificate (Theorem 1).** 2592 tets, 5544 CR dofs, 8 s-windows, arb prec 128, Rump BIT 46. Reproduces: `python m3_certify.py`. |
| **`m3p_certify.py`** | **Congruence certificates (Theorems 2–4).** N𝔭=5, 9, 13; per-row Rump equilibration. |
| `CONGRUENCE.md` | Multi-cusp criterion: Lemma D0 (trace exactness), Lemma R (interval rank-one), Facts A–B (coset structure). |
| `0cuspchecks.md` | Cusp-0 audit trail for N𝔭=5: conjugated stabilizer, fold-compatible tiling, height-preserving chimney maps. |
| `cr_prototype.py` | Float CR pipeline + box-model validation. |
| `fem_prototype.py` | FEM setup and assembly. |
| `congruence_prototype.py` | Coset combinatorics and index arithmetic. |
| `theorem4postmortem.md` | Equilibration recipe for 33k-DOF systems (per-row radii, diagonal SAS preconditioning). |
| `PROGRESS_LADDER.md` | Live log: float margins, mesh convergence, per-level timings. |
| `HANDOFF_LADDER.md` | Roadmap for N𝔭 = 13, 17, 25, … |

**Theorems Certified:**
1. **Theorem 1 (Level 1, Picard):** λ₁(PSL(2,ℤ[i])\ℍ³) ≥ 1. Mesh 12×6×6 → 5544 dofs.
2. **Theorem 2 (N𝔭=5):** λ₁(Γ₀((2+i))\ℍ³) ≥ 1. Index 6, vol ≈ 1.83, **beyond trace-formula reach** (B ≈ 1.9 > 1). Mesh 8×4×3 × 6 copies → 28,400 dofs.
3. **Theorem 3 (N𝔭=9):** λ₁(Γ₀((3))\ℍ³) ≥ 1. Index 10, torsion-free, 10 copies → 26,532 dofs.
4. **Theorem 4 (N𝔭=13):** λ₁(Γ₀((3+2i))\ℍ³) ≥ 1. Index 14, 14 copies → 33,176 dofs (per-row Rump equilibration required).

**Status:** ✅ **DONE.** All four theorems certified and reproducible.

---

### 🔷 Dual Certification: Cusp-Defect Analysis (Hejhal Method)

**Path:** `dual_certification_/`

**Goal:** Approximate eigenvalues and eigenfunctions via truncated cusp expansions, with rigorous defect bounds.

| File | Role |
|------|------|
| **`theorem_DK.tex`** | Full paper: Assumption H (Hecke growth), Lemma K (Bessel upper bound, C_K=1), Theorem D(K) (defect bound with C₁, C₂). |
| **`lemma_K.py`** | Arb-enclosable majorant S_{M,∞} for Fourier tails under Hecke growth. `python lemma_K.py --test --bench --constants`. |
| `constants_DK.md` | Reference numerics: C₁ ≈ 1.3·10⁶, C₂ ≈ 11.7 (Picard); tail decay ~ 10⁻⁸³ at M=400. |
| `plot_kappa_vs_M.py` | Conditioning-number diagnostic: κ_proxy vs truncation M. |

**Rungs (status per `ROADMAP.md`):**
- Rung 0–1: **DONE** — Lemma K and Theorem D(K), elementary C_K proof.
- Rung 2–3: **DONE** — FEM half-line, two-cusp Hejhal when η small, Krawczyk N=5.
- Rung 4: **Infra done** — N=5 dual pipeline; equation 13 (full cert) open.

---

### 🌀 Eisenstein Numbers: FEM for PSL(2, ℤ[ω])

**Path:** `eisenstein_numbers/`

**Goal:** Independent FEM proof for the Eisenstein–Picard group (parallel to `independent_exclusion/` for Picard).

| File | Role |
|------|------|
| **`PROOF.md`** | Theorems A, B₁, B₂ with FEM architecture. |
| **`cert_omega.py`** | **Level 1 certificate.** 6×3 mesh (P₃ FE), 8 s-windows, `python -u cert_omega.py 6 3`. |
| `geometry_fund.py` | Fundamental domain assembly and validity checks. |
| `cr_omega.py` | Float CR pipeline for ω. |
| `GEOMETRY.md` | EGM P₃ cell coordinates and orientation. |
| `framework.py` | STF + CE gates for ℤ[ω] (see `AUDIT.md` for load-bearing math). |
| `smoke_test.py` | Ring and coset checks. |

**Status:** ✅ **DONE.** Theorem A certified; Theorems B₁–B₂ (congruence rungs on ω) float-ready, cert open.

---

## Reference materials and shared trace engine

Stable citations and source roles are in [docs/REFERENCES.md](docs/REFERENCES.md).
Local derivations are in [docs/ANALYTIC_DERIVATIONS.md](docs/ANALYTIC_DERIVATIONS.md)
and [docs/SYSTOLES.md](docs/SYSTOLES.md). Private OCR files are not required.

| Path | Role |
|---|---|
| `groups/` | Explicit group/class data, exact ring arithmetic, systole and class-number checks |
| `fields/` | Characters, field constants, and prime splitting |
| `core/` | Exact B-spline derivative, Arb quadrature, assembly, and proof gates |
| `scripts/feasibility/` | The four portable attached screening/evidence scripts |
| `examples/field_screen.py` | Arb mechanical screen; no spectral conclusion |
| `certificates/` | Frozen historical and current regression reports |
| `tests/test_trace_core.py` | Analytic regression, prime splitting, proof-gate, and portability checks |
| `RIGOR_GAPS.md` | Current proof ledger, including open new-field inventories |

## 🎯 How to Navigate

### **Just the headlines?**
→ Read `README.md` (main theorem and proof dependency graph).

### **I want the trace-formula proof.**
→ Start with `picard_stf.py` (the engine). Check: `verify_group_data.py` (constants), `verify_matthies.py` (Weyl validation), `bianchi_omega_arb.py` (certified B < 1).

### **I want the FEM proof (independent of trace formula).**
→ Start with `independent_exclusion/PROOF.md`. Reproduce with `python independent_exclusion/m3_certify.py` (level 1) or `python -c "from independent_exclusion.m3p_certify import certify; certify(...)"` (congruence).

### **I want the Hejhal eigenvalue approximation.**
→ See `dual_certification_/` and its `ROADMAP.md`. Lemma K bounds are in `lemma_K.py`.

### **I want Eisenstein–Picard (ℤ[ω]).**
→ Trace formula: `bianchi_omega_arb.py`. FEM: `eisenstein_numbers/cert_omega.py`.

### **I want to understand the congruence ladder.**
→ Read `independent_exclusion/PROGRESS_LADDER.md` (status), then `CONGRUENCE.md` (architecture), then `m3p_certify.py` (implementation). Theorems 2–4 reproduce with one-liners.

---

## 🏆 Key Validation Strategies

1. **Matthies fingerprinting** (`verify_matthies.py`): Reconstruct every Weyl-law coefficient from engine constants; match exactly. Confirms elliptic sector is complete.
2. **Dual-method agreement**: Trace formula (B < 1) + independent FEM (matrix PSD) both certify λ₁ ≥ 1 on level 1. Only shared input: fundamental domain.
3. **Per-row Rump equilibration** (`theorem4postmortem.md`): Scale matrix diagonally, apply radii per row, not as uniform shift. Enables 33k-DOF systems.
4. **Reference-cell principle** (`CONGRUENCE.md` Facts A–B): Congruence covers = 6, 10, or 14 isometric copies of level-1 mesh. Index enters via exact 𝔽_N combinatorics, not floating-point magic.

---

## 📋 Reproducibility Checklist

- [ ] Trace formula (Picard): `python picard_stf.py`
- [ ] Weyl validation: `python verify_matthies.py`
- [ ] Trace formula (Eisenstein): `python bianchi_omega_arb.py`
- [ ] FEM level 1: `python independent_exclusion/m3_certify.py`
- [ ] FEM N𝔭=5: `python -c "from independent_exclusion.m3p_certify import certify; certify(8,4,3,level='(2+i)')"` (takes ~15 min)
- [ ] FEM N𝔭=9: `python -c "import independent_exclusion.m3p_certify as m; m.certify(6,3,3,level='(3)')"` (takes ~40 min)
- [ ] FEM N𝔭=13: `python -c "import independent_exclusion.m3p_certify as m; m.NU_STAR=1.001; m.NWIN=16; m.certify(8,3,2,level='(3+2i)')"` (takes ~90 min)
- [ ] Eisenstein FEM: `cd eisenstein_numbers && python cert_omega.py 6 3`

---

## 🚀 What's Next (Roadmap)

- **Trace formula:** LP-optimized test functions (Booker–Strömbergsson), Eisenstein–Picard at higher levels.
- **FEM:** Congruence ladder on ω (N𝔭=3, 7, …), mesh refinement for tighter margins.
- **Hejhal:** Full η-dependent certified eigenvalue enclosures (rung 4 of `dual_certification_/`).

---

**Last Updated:** 2026-07-12  
**Certified by:** Arb + Rump (IEEE-754, per-row equilibration)  
**Author:** @Sigmoidd
