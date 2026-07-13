# Morning brief — Rung 4 residual path

**Updated:** next three steps — \(V_h^{P1,\mathrm{per}}\) spectrum, production two-cusp residual, Route A Arb scaffold  
**Hard map (do not change):**  
`F5_gluing=true, FEM_exclusion_01=true, two_cusp_hejhal_scan=true, defect_bound_applied=true,`  
`eta_le_eta0=false, width_lt_tol=false, counting_certified=false`  
⇒ **`pipeline_ok=true`, `rung4_certified=false`** ← still correct.

---

## 0. Snapshot

| Item | Status |
|------|--------|
| C1 1.1+1.2+1.6 | **C1≈1.98e3**, η₀≈6.4e-8, gap **~3.9** orders |
| Hybrid Fourier δ | best ~0.09 (periodize M=48) |
| \(J_h\) / periodic ID | PASS |
| **Q/M on \(V_h^{P1,\mathrm{per}}\)** | **shipped** — first pos Neumann ~**6.83** |
| **Production two-cusp residual** | **shipped** — δ_∞~0.17, rel~1.6e-4 (M=40) |
| **Route A Arb scaffold** | **shipped** — N(1) candidate **[0,1]** YELLOW |

---

## 1. \(V_h^{P1,\mathrm{per}}\) spectrum (step 1)

**Code:** `v_h_p1_per_spectrum.py`  
Reuses `assemble_level_p` (cross-copy pairings, free self/top).

| mesh | dofs | λ₁^N | first pos N | λ₁^D |
|------|-----:|-----:|------------:|-----:|
| 4×2×2, 6 copies | 4816 | ~0 | **6.825** | 19.58 |

Checks: |Q1|~0, t∞(1)=1/2, t0(1)=5/2, connected.  
**Language:** engineering CR spectrum on conforming-periodic space — **not** certified dual upper bound (needs CR GLB / residual / counting).

```bash
python v_h_p1_per_spectrum.py --N1 4 --N2 2 --N3 2
```

---

## 2. Production two-cusp residual (step 2)

**Code:** `production_hejhal_residual.py`  
Two-cusp near-kernel + true H³ PAIRINGS δ_aut (+ optional Poincaré).

| quantity (M=40, r=6, Y0=0.8, periodize) | value |
|----------------------------------------|------:|
| two-cusp rel σ_min/σ_max | 1.58e-4 |
| κ_eq | ~6e3 (ill-conditioned; model) |
| δ_aut(∞) | **0.173** |
| δ_aut(0) | 0.131 |
| τ_proxy | 9.6e-4 |
| η_proxy | ~0.174 |

vs single-cusp periodize M=48 (δ~0.09): two-cusp model residual is **same order**, slightly worse δ, much better collocation rel. Still **~6 orders** short of η₀~6e-8.

```bash
python production_hejhal_residual.py --M 40 --r 6.0
# optional: --r-scan 5.5:7.0:7
```

---

## 3. Route A Arb / D–N scaffold (step 3)

**Code:** `route_A_arb_scaffold.py`

| deliverable | status |
|-------------|--------|
| Float Picard D/N | done |
| Relative mid/rad on Q,M | done (placeholder for Arb) |
| Interval matvec residual | done (~1e-10 on first pos mode) |
| KAPPA1 GLB sketch | done |
| Candidate N(λ) integer interval shape | **YELLOW** e.g. N(1)=**[0,1]** |
| Truncation Δ constants | scaffold only |
| Certified counting | **still false** |

```bash
python route_A_arb_scaffold.py --N1 6 --N2 3 --N3 3
```

---

## 4. Master checklist

### Done
- [x] C1 geometry + §1.6 rewrite
- [x] Hybrid multi / periodize / hybrid_scan_M
- [x] \(J_h\) + \(V_h^{P1,\mathrm{per}}\) identification
- [x] **Spectrum on \(V_h^{P1,\mathrm{per}}\)**
- [x] **Production two-cusp residual API**
- [x] **Route A Arb/D–N scaffold + candidate N(λ)**

### Open (dual GREEN blockers)
- [ ] η ≲ η₀ (~6e-8) — residual still ~0.1
- [ ] width &lt; 0.1 with production residual
- [ ] counting_certified (Arb D/N + proved Δ_trunc → N=[0,0] on (1,λ₁))
- [ ] CR GLB / Rump for continuum-safe bounds
- [ ] rung4_certified

---

## 5. Next priorities

| Pri | Action |
|----:|--------|
| 1 | Drive residual: true Hejhal iterate / denser periodization / eigen-r search |
| 2 | CR GLB on \(V_h^{P1,\mathrm{per}}\) first positive mode |
| 3 | Arb tet assembly for Route A (replace relative mid/rad) |
| 4 | Keep hard map |

---

## 6. Reproduce

```bash
cd dual_certification_
python v_h_p1_per_spectrum.py --N1 4 --N2 2 --N3 2
python production_hejhal_residual.py --M 40 --r 6.0
python route_A_arb_scaffold.py --N1 6 --N2 3 --N3 3
python p1_per_jh_bridge.py
python validate_gap_numbers.py
python lemma_K.py --test
```

---

## 7. Bottom line

All three next steps landed as **engineering infrastructure**:
1. Conforming-periodic CR spectrum (λ₊≈6.83)  
2. Two-cusp production-shaped residual (δ~0.17, still far from η₀)  
3. Route A candidate counting enclosure shape (N(1)=[0,1] YELLOW)  

Dual remains **not certified**. Hard map unchanged.

*Path: `dual_certification_/morningbrief.md`*
