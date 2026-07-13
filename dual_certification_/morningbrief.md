# Morning brief — Rung 4 residual path

**Updated:** next three — eigen-r residual drive, CR GLB on \(V_h^{P1,\mathrm{per}}\), Arb tet assembly  
**Hard map:** `pipeline_ok=true`, `rung4_certified=false` (unchanged).

---

## 0. Snapshot

| Step | Result |
|------|--------|
| C1 1.1+1.2+1.6 | C1≈1.98e3, η₀≈6.4e-8, gap ~3.9 orders |
| Best Fourier δ | ~0.09–0.16 (periodize / two-cusp) |
| **Eigen-r drive** | residual **flat in r** (η≈0.162, best r≈5.5) |
| **CR GLB on \(V_h^{P1,\mathrm{per}}\)** | λ₊≈6.825 → **GLB≈6.451** (factor 0.945) |
| **Arb tet Route A** | 96/96 vol>0 certified; Q/M mid+rad assembled |

---

## 1. Eigen-r residual drive

**Code:** `eigen_r_residual_drive.py`  
Grid + golden refine on two-cusp η = δ_aut + τ_proxy.

| (M=28, periodize, r∈[5.5,7.5]) | value |
|--------------------------------|------:|
| δ_aut | **0.161** (essentially flat vs r) |
| best η | **0.162** at r≈5.5 |
| factor vs grid start | ~1.00× |

**Honest conclusion:** for the present model operator, **r-search does not buy orders**. Residual is automorphy-limited (PAIRINGS), not a sharp eigen-r miss. Need a different trial (true Hejhal iterate / conforming periodization), not denser r-grid.

```bash
python eigen_r_residual_drive.py --M 28 --n-grid 5 --refine-iters 4
```

---

## 2. CR GLB on \(V_h^{P1,\mathrm{per}}\)

**Code:** `cr_glb_p1_per.py`  
\[
\lambda \ge \frac{\lambda_h}{1+\kappa_1^2 h_{\max}^2\lambda_h},\quad
\kappa_1=\mathrm{KAPPA1}\approx 0.331.
\]

| mode | λ_h | GLB |
|------|----:|----:|
| first pos Neumann | **6.825** | **6.451** |
| λ₁^D | 19.58 | 16.79 |

GLB/λ_h ≈ **0.945** at this mesh (h_max≈0.278).  
**Language:** engineering GLB shape on multi-copy CR — full CR hypotheses for Γ₀ multi-copy not re-certified here (cf. m3p for exclusion). Not dual GREEN.

```bash
python cr_glb_p1_per.py --N1 4 --N2 2 --N3 2
```

---

## 3. Arb tet assembly (Route A)

**Code:** `route_A_arb_tet.py`  
Per-tet Arb: vol>0, h_T, y^{-1}/y^{-3} integral enclosures → Q/M mid+rad.

| (mesh 4×2×2, Y=1.25) | value |
|----------------------|------:|
| tets vol>0 certified | **96/96** |
| h_max | 0.421 |
| worst rel rad wQ / wM | 0.17 / 0.48 |
| ‖Q_rad‖_F | 1.63 |
| first pos λ^N (mid) | 7.69 |
| status | **YELLOW** |

Min/max height bounds are conservative (large weight radii). Next: tighter Taylor moments (as in m3_certify.ell_entries) + Rump on shifted pencils.

```bash
python route_A_arb_tet.py --N1 4 --N2 2 --N3 2
```

---

## 4. Checklist

### Done
- [x] C1 geometry + §1.6
- [x] Hybrid / periodize / two-cusp residual APIs
- [x] \(V_h^{P1,\mathrm{per}}\) ID, spectrum, **CR GLB sketch**
- [x] **Eigen-r residual drive** (flat — documented)
- [x] **Arb tet Q/M mid+rad** for Route A core

### Open (dual GREEN)
- [ ] Residual η ≲ 6e-8 (still ~0.16)
- [ ] Tighter Arb weights + Rump D/N
- [ ] counting_certified with proved Δ_trunc
- [ ] rung4_certified

---

## 5. Next

| Pri | Action |
|----:|--------|
| 1 | True production Hejhal / Poincaré over larger group / conforming FEM trial residual |
| 2 | Taylor-moment Arb weights (m3 style) + Rump PSD path for Route A |
| 3 | Feed GLB(6.45) into dual language carefully (not as certified λ₁) |
| 4 | Keep hard map |

---

## 6. Reproduce

```bash
cd dual_certification_
python eigen_r_residual_drive.py --M 28 --n-grid 5 --refine-iters 4
python cr_glb_p1_per.py --N1 4 --N2 2 --N3 2
python route_A_arb_tet.py --N1 4 --N2 2 --N3 2
python lemma_K.py --test
```

---

## 7. Bottom line

Three more engineering layers closed: **r-drive** shows residual is flat (not eigen-r limited); **CR GLB** lowers first positive discrete mode 6.83→6.45; **Arb tet assembly** certifies volumes and builds explicit Q/M radii (YELLOW). Dual still blocked by **η~0.16 ≫ η₀~6e-8** and uncertified counting. Hard map unchanged.

*Path: `dual_certification_/morningbrief.md`*
