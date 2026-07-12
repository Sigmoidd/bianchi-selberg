# Does independent exclusion FEM generalize to PSL(2,ℤ[ω])?

**Date:** 2026-07-11 · **Status:** M0 float probe — **promising, not certified**  
**Code:** `fem_omega_m0.py` · **Picard analogue:** `independent_exclusion/fem_prototype.py`

---

## Short answer

| Layer | Generalizes? | Notes |
|-------|----------------|-------|
| Criterion (Lax–Phillips + zero mode) | **Yes, abstractly** | Same lemmas; plug in \|T\|, dual lattice, min \|c\| |
| Level-1 float margin (M0) | **Yes, empirically** | min μ ≈ **1.62** at 12×12×6, Y=1.25 (Picard ~7.3) |
| CR + Theorem G1 + Rump cert | **Not yet** | Needs mesh admissibility (Lemma G), CR assembly, windows |
| Congruence Γ₀(𝔫) ⊂ PSL(2,ℤ[ω]) | **Not started** | Needs ℙ¹(ℤ[ω]/𝔫), multi-cusp Facts A–B for Eisenstein |

We did **not** reproduce Theorems 1–4. We ran the **same first experimental step** (M0) that justified the Picard FEM ladder: discrete positivity margin of the exclusion form under Neumann relaxation.

---

## What must change vs Picard (ℤ[i])

| Input | Picard (done) | Eisenstein (this probe) |
|-------|----------------|-------------------------|
| Ring | ℤ[i] | ℤ[ω] |
| Cusp lattice area | \|T\|=1/2 | \|T\|=√3/6 (= covol Λ / 3) |
| Base of K | rectangle | parallelogram slice z=u+vω |
| Floor | √(1−\|z\|²) | same (sphere) |
| Mode bound | 4π²Y² > 1 | 4π²Y²/(2/√3)² ≈ 46 > 1 at Y=1.25 |
| Shimizu | min \|c\|²≥1 | same for 0≠c∈ℤ[ω] |
| Pairings / 24-split | T1,R,TiR,S for ℤ[i] | **Different** face group for full cert |
| β in 𝒜_s | 2(1−s)/Y² | (1−s)/(\|T\| Y²) (same rule) |

**Unchanged:** spectral atomicity on (0,1) (Friedman 3.8.1 for general cofinite Kleinian), ODE in the cusp, constraint ℒ_s, CR theory G1 (field-agnostic once geometry is admissible).

---

## M0 results (float, conforming Q1)

```
Y = 1.25, lams in (0,1)
mesh  8×8×4:  min mu ≈ 1.82
mesh 12×12×6: min mu ≈ 1.62   (refinement lowers upper bias)
checks: t(1)=|T| exact, |S·1|~0, mode-bound OK
```

| Comparison | Picard M0 | Eisenstein M0 |
|------------|-----------|-----------------|
| min μ | ~7.3 | ~1.6 |
| Verdict | comfortable | **positive, thinner** |

Same pattern as congruence ladder: thinner margin ⇒ harder Rump/CR, not automatic failure. Picard M0 ≫ 1 enabled easy certification; here headroom is closer to the N𝔭=13 situation.

**Caveats (honest):**

1. Domain is the **EGM truncated core** K_Y over P_3 (GEOMETRY.md; `build_P3_cell`). Neumann relaxation is intentional (as in DESIGN.md); side-pairing maps to drop it are recorded in GEOMETRY.md §4. (Legacy R_comp parallelogram is M0-only; it hits \|z\|=1.)
2. vol(K) self-check uses **vol(K_Y)=vol(F)−\|T\|/(2Y²)**; CR float G1 and **interval cert** both pass 8/8 at rho=55 (N_tri=6, N3=3) via `cert_omega.py` (m3/m3p port).
3. No CR, no arb, no Rump — **evidence only**.

---

## Roadmap if we continue (same steps as Picard)

1. **M0** ✅ float margin (this file’s experiment).  
2. **Geometry note** — ~~freeze EGM domain~~ done (`GEOMETRY.md` / `geometry_fund.py`); keep Lemma 0 + fold in the paper writeup.  
3. **CR prototype** — port `cr_prototype` / 24-split or hex-compatible mesh; check Lemma G on floor.  
4. **m3-style certify** — windowed G1 + equilibrated Rump (reuse Thm 4 Rump path).  
5. **Congruence (later)** — ℙ¹(ℤ[ω]/𝔫), multi-cusp; analogue of Thms 2–4.

---

## Relation to B < 1 (trace formula)

| Method | Eisenstein level 1 status |
|--------|---------------------------|
| STF positivity B < 1 | Done in software (`framework.py` step 5): B ∈ [0.525, 0.544] |
| FEM independent exclusion | **M0 only** — margin ~1.6, cert not built |

They are **independent** routes to λ₁ ≥ 1. FEM is what makes Thms 2–4 “trace-formula free”; extending that route to ω is exactly this program.

---

## Run

```powershell
cd eisenstein_numbers
python -u fem_omega_m0.py
```
