# Morning brief — Rung 4 residual path

**Updated:** larger Poincaré + conforming residual; Taylor+Rump Route A; dual GLB language  
**Hard map:** `pipeline_ok=true`, `rung4_certified=false` (unchanged).

---

## 0. Snapshot

| Item | Result |
|------|--------|
| C1 1.1+1.2+1.6 | C1≈1.98e3, η₀≈6.4e-8, gap ~3.9 orders |
| **Larger Poincaré** | H-only best **δ≈0.055**; bigger word groups **worse** |
| **Conforming \(V_h^{P1,\mathrm{per}}\)** | λ₊=**6.825**, GLB=**6.451**, \(J_h^{\mathrm{cross}}=0\) |
| **Taylor Arb weights** | I1 rel rad **7e-4** (was 0.17 min/max) |
| **Rump Dirichlet** | **PSD PASS** for t up to ~first D (~18.6) |
| Dual language | **`dual_glb_language.md`** — no false [6.45,6.83] cert |

---

## 1. Larger Poincaré + conforming residual

**Code:** `larger_poincare_residual.py`

### Poincaré group size (M=32, hybrid multi, periodize)

| group | \|F\| | δ_aut | vs H-only |
|-------|------:|------:|----------:|
| **H_only** (id,T1±,R,TiR) | 5 | **0.055** | 1× |
| words ≤1 (+S) | 6 | 0.46 | 0.12× worse |
| words ≤2 | 23 | 5.5 | much worse |

**Conclusion:** naive larger finite Poincaré (without group closure / height control) **hurts**. Best remains short horizontal periodization. Residual still ~0.05 ≫ η₀~6e-8.

### Conforming CR trial

| quantity | value |
|----------|------:|
| first pos λ_h | **6.825** |
| Rayleigh check | 6.825 |
| GLB sketch | **6.451** |
| \(J_h^{\mathrm{cross}}\) | **0** |

Correct dual **upper-flavor** discrete trial space (not free Neumann).

```bash
python larger_poincare_residual.py --M 32
```

---

## 2. Taylor-moment Arb + Rump (Route A)

**Code:** `route_A_rump.py`

| | min/max height (old) | **Taylor p=4** |
|--|---------------------:|---------------:|
| worst rel I1 | 0.17 | **7.4e-4** |
| worst rel I3 | 0.48 | **2.3e-2** |
| ‖Q_rad‖_F | 1.63 | **3.1e-3** |

**Rump PSD (BIT 46):**

| pencil | t | rump |
|--------|--:|:----:|
| Neumann Q | 0 | FAIL (constant null — expected) |
| Dirichlet Q | 0 | **PASS** |
| Dirichlet Q−tM | up to ~λ₁^D | **PASS** |

Dirichlet Rump succeeds near the first Dirichlet eigenvalue on the float mid matrix. Still **YELLOW**: radii not folded into interval eigenproof; truncation Δ open; **not** counting_certified.

```bash
python route_A_rump.py --N1 4 --N2 2 --N3 2 --taylor-p 4
```

---

## 3. Dual language around GLB ≈ 6.45

**Doc:** `dual_glb_language.md`

| Allowed | Forbidden |
|---------|-----------|
| FEM: λ₁ ≥ 1 certified | “certified λ₁ ∈ [6.45, 6.83]” |
| engineering λ_h⁺ ≈ 6.83 | “GLB proves continuum lower 6.45” |
| GLB sketch 6.45 under CR hypotheses | free Neumann as dual upper bound |
| rung4_certified=false | flipping hard flags |

---

## 4. Checklist

### Done
- [x] C1 + residual APIs + V_h spectrum / J_h / GLB
- [x] Eigen-r drive (flat)
- [x] **Larger Poincaré compare** (H-only wins)
- [x] **Conforming residual language**
- [x] **Taylor weights + Rump Dirichlet PASS**
- [x] **dual_glb_language.md**

### Open (dual GREEN)
- [ ] η ≲ 6e-8 (best δ still ~0.05)
- [ ] Interval eigen / Rump with radii + Δ_trunc
- [ ] counting_certified
- [ ] rung4_certified

---

## 5. Next

| Pri | Action |
|----:|--------|
| 1 | Production Hejhal iterate (not finite Poincaré toys) |
| 2 | Fold Q_rad into interval residual + Rump-with-radius |
| 3 | Multi-copy CR Arb checks (m3p style) for GLB on Γ₀ |
| 4 | Keep hard map |

---

## 6. Reproduce

```bash
cd dual_certification_
python larger_poincare_residual.py --M 32
python route_A_rump.py --N1 4 --N2 2 --N3 2 --taylor-p 4
# dual language: dual_glb_language.md
python lemma_K.py --test
```

---

## 7. Bottom line

**Poincaré:** short horizontal group best (δ~0.055); larger words worse.  
**Conforming space:** λ_h=6.825, GLB sketch=6.451, jumps 0 — dual language frozen.  
**Route A:** Taylor weights ~200× tighter radii; **Dirichlet Rump PSD passes**.  
Still **not dual-certified** (η, full interval spectrum, counting). Hard map unchanged.

*Path: `dual_certification_/morningbrief.md`*
