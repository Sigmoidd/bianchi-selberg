# Morning brief — Rung 4 residual path

**Updated:** Hejhal iterate, Rump-with-radius, multi-copy Arb/GLB  
**Hard map:** `pipeline_ok=true`, `rung4_certified=false` (unchanged).

---

## 0. Snapshot

| Item | Result |
|------|--------|
| C1 production | ~1.98e3, η₀~6.4e-8, gap ~3.9 orders |
| **Hejhal iterate** | best **δ≈0.036**, η≈0.037 (M=28, periodize+reproject) |
| **Rump+radius** | Dirichlet **PASS** for t≲0.5 λ₁^D; fails near λ₁^D |
| **Multi-copy Arb GLB** | 384/384 vol>0; λ_h=6.825, GLB=**6.451** |
| Dual language | `dual_glb_language.md` — no false cert interval |

---

## 1. Production Hejhal iterate

**Code:** `hejhal_iterate.py`  
Outer loop: golden r-refine on collocation rel → two-cusp near-kernel → optional periodize L² reproject → PAIRINGS δ_aut.

| iter | r | δ_aut | η |
|-----:|--:|------:|--:|
| 0 | 6.74 | **0.0361** | **0.0369** |
| 1 | 7.26 | 0.0361 | 0.0370 |
| 2 | 7.63 | 0.0361 | 0.0370 |

**Best η ≈ 0.037** (was ~0.16 two-cusp single shot, ~0.055 H-only periodize).  
Still **~5.8 orders** above η₀~6e-8. Iterate improves vs static residual but plateaus.

```bash
python hejhal_iterate.py --M 28 --iters 3 --r0 6.0 --periodize
```

---

## 2. Rump-with-radius (fold Q_rad)

**Code:** `route_A_rump_radius.py`  
Extra = ‖R_Q + |t| R_M‖_F; Rump certifies λ_min(A_mid) > extra ⇒ interval hull PSD.

| Dirichlet t | rump+radius | margin |
|------------:|:-----------:|-------:|
| 0 | **PASS** | 1.1e-2 |
| 1 | **PASS** | 1.0e-2 |
| 0.5 λ₁^D | **PASS** | 3.2e-3 |
| 0.9 λ₁^D | **FAIL** | −3e-3 |

Taylor I1 rel ~7e-4. Interval residual on first D mode dominated by radius (~0.036).  
**YELLOW:** not counting_certified; near-eigenvalue Rump still open.

```bash
python route_A_rump_radius.py --N1 4 --N2 2 --N3 2 --taylor-p 4
```

---

## 3. Multi-copy CR Arb + GLB (Γ₀)

**Code:** `multi_copy_cr_arb_glb.py`

| check | status |
|-------|--------|
| Reference tet vol>0 (Arb, oriented) | **384/384 GREEN** |
| h_max^U Arb | 0.278 |
| Multi-copy Q1, cusp areas | **GREEN float** |
| first pos λ_h / GLB | **6.825 / 6.451** |
| multi-copy κ₁ theory | **YELLOW** |
| continuum GLB as λ₁ lower | **RED** (only FEM ≥1 certified) |

```bash
python multi_copy_cr_arb_glb.py --N1 4 --N2 2 --N3 2
```

---

## 4. Checklist

### Done
- [x] Hejhal outer iterate (δ~0.036)
- [x] Rump-with-radius on Dirichlet pencils
- [x] Multi-copy Arb geometry + GLB sketch + theory labels
- [x] dual_glb_language.md

### Open (dual GREEN)
- [ ] η ≲ 6e-8 (best still ~0.037)
- [ ] Rump+radius near spectral edge / full interval eigen
- [ ] multi-copy CR κ certified → continuum GLB
- [ ] counting_certified, rung4_certified

---

## 5. Next

| Pri | Action |
|----:|--------|
| 1 | Deeper Hejhal (true pullback S, larger M, multi-start r) |
| 2 | Radius-aware eigenvalue enclosures (not just PSD) |
| 3 | m3p-style Arb κ on multi-copy (long) |
| 4 | Keep hard map |

---

## 6. Reproduce

```bash
cd dual_certification_
python hejhal_iterate.py --M 28 --iters 3
python route_A_rump_radius.py --N1 4 --N2 2 --N3 2
python multi_copy_cr_arb_glb.py --N1 4 --N2 2 --N3 2
```

---

## 7. Bottom line

**Iterate** cut residual to **δ~0.036** (still far from η₀). **Rump+radius** certifies Dirichlet PSD well below λ₁^D. **Multi-copy Arb** greens reference geometry; GLB 6.45 stays engineering. Dual not certified. Hard map unchanged.

*Path: `dual_certification_/morningbrief.md`*
