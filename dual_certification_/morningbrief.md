# Morning brief — Rung 4 residual path

**Updated:** §1.6 **accepted** (proof rewrite) + larger-M hybrid scan  
**Hard map (do not change):**  
`F5_gluing=true, FEM_exclusion_01=true, two_cusp_hejhal_scan=true, defect_bound_applied=true,`  
`eta_le_eta0=false, width_lt_tol=false, counting_certified=false`  
⇒ **`pipeline_ok=true`, `rung4_certified=false`** ← still correct.

---

## 0. Snapshot

| GREEN / shipped | Note |
|-----------------|------|
| F5 + FEM (0,1) + model Hejhal scan | pipeline_ok |
| Conditioning \(b_{\mathrm{eq}}\approx 0.626\) | κ_eq 42–155 |
| C1 1.0–1.2 geometry | y_min=1/√2, box D_Y |
| **C1 §1.6** | **ACCEPTED** — λ-free \(A_{\mathrm{bdry}}\) |
| Hybrid multi + periodize δ_aut | best ~0.09 at M=48 |
| Residual / hybrid M-tables | see §2–3 |

**Production C1 (sharp, r=6):** **C1 ≈ 1.98e3**, **η₀ ≈ 6.37e-8**, gap @ η=5e-4 ≈ **3.9 orders**.

| Goal | Need | Gap @ η=5e-4 |
|------|------|--------------|
| width &lt; 0.1 (Y0=1.5) | η ≲ **5e-5** | ~1 order (width only) |
| η ≤ η₀ (production) | η ≲ **6e-8** | **~3.9 orders** |
| η ≤ η₀ (conservative) | η ≲ **2e-13** | **~9.3 orders** |

---

## 1. C1 geometry — including §1.6 **proved**

\[
A_{\mathrm{bdry}}=C_{\mathrm{tr}}C_{\mathrm{Sob}}A_{\mathrm{met}}
\quad\text{(λ-free)},
\qquad
A_{\mathrm{res}}=(1+\lambda)(1+C_P A_{\mathrm{met}}).
\]

| # | Change | Factor | Status |
|---|--------|-------:|:------:|
| 1.0 | \(C_K=1\), \(C_{\mathrm{tr}}=3/h_{\min}\) | ~14× vs cons. | **Done** |
| 1.1 | \(y_{\min}=1/\sqrt2\) | 2.39× | **Done** |
| 1.2 | Box \(D_Y\) | 2.29×; **5.24×** w/ 1.1 | **Done** |
| **1.6** | Drop \((1+\sqrt\lambda)\) from \(A_{\mathrm{bdry}}\) | **7.083×** | **Done (proved)** |
| 1.3–1.5 | mesh / directional / Sobolev | 1.2–4× each | **Open** |

**Proof (theorem_DK.tex Remark `rem:no-double-sqrt`):**  
Step 1 \(H^1\) corrector is λ-free. Step 2 puts all λ in \(A_{\mathrm{ell}}=1+\lambda\subset A_{\mathrm{res}}\). Earlier \((1+\sqrt\lambda)\) in \(A_{\mathrm{bdry}}\) double-counted.  
Code default: `bdry_sqrt_lam=False`. Legacy: `True` for regression.

**Stack 1.0–1.2+1.6 vs original cons C1~1e6:** roughly **~500×**. Residual still dominates.

---

## 2. δ_aut — hybrid multi + larger M

**Shipped:** `delta_aut_pairing.py`, `hybrid_scan_M.py`.

### Best single-point (M=48, jump_w=2, n_face=16)

| mode | δ_aut | factor vs col |
|------|------:|--------------:|
| collocation | 3.39 | 1× |
| multi hybrid | 0.443 | 7.6× |
| **periodize** | **0.088** | **38×** |

### Larger M (jump_w=2, n_face=16) — absolute δ **worsens**

| M | col δ | multi δ | periodize δ | fac multi | fac per |
|--:|------:|--------:|------------:|----------:|--------:|
| 48 | 3.4 | 0.44 | **0.088** | 7.6 | 38 |
| 64 | 16 | 0.86 | 0.17 | 19 | 93 |
| 100 | 26 | 1.7 | 0.34 | 15 | 75 |
| 200 | 149 | 20 | 4.1 | 7.4 | 37 |

**Honest conclusion:** hybrid/periodize still beat collocation by **~10–90×**, but **raising M alone does not drive δ→0**. Extra modes enrich the null-space of the collocation pin; face-jump constraints do not scale with \(n_{\mathrm{modes}}\) unless face sampling and jump weight grow carefully. Extra n_face / jump_weight trials (36–64 faces, w=5–10) did **not** restore the M=48 periodize floor.

Still **O(10⁻¹)** best vs targets **10⁻⁶–10⁻⁹**.

| # | Status |
|---|:------:|
| 2.1 hybrid multi | **Done** |
| 2.2b Poincaré H | **Done** |
| 2.3 true σ₀ / eigen-r | **Open** |
| 2.5 \(V_h^{P1,\mathrm{per}}\) | **Open** |

---

## 3–4. τ scan / η₀ / Route A

- Residual M-scan (rel): M:100→800 ~3.7×; τ_svd headline not SVD-kernel zero.
- `U_norm` API done; inventing U banned.
- Route A: trunc scaffold RED; counting not certified.

---

## 5. Checklist

### Done
- [x] Hard map honest; pipeline_ok / not certified
- [x] 1.1+1.2+**1.6 proved** → C1 **1.98e3**, η₀ **6.4e-8**
- [x] Hybrid multi + periodize; hybrid_scan_M table
- [x] Larger-M hybrid scan (absolute δ worsens with M)
- [x] Route A trunc scaffold

### Open
- [ ] δ_aut ≲ 1e-5 (need true near-automorphic trial / conforming space)
- [ ] counting_certified
- [ ] eta_le_eta0 / width_lt_tol / rung4_certified

---

## 6. Next priorities

| Pri | Action |
|----:|--------|
| 1 | Conforming \(V_h^{P1,\mathrm{per}}\) or CR \(J_h\) trial (not free Neumann) |
| 2 | True two-cusp σ₀ Hejhal iterate at eigen-r (not model collocation pin) |
| 3 | Optional 1.3–1.5 C1; Route A Arb D/N |
| 4 | Keep hard map |

---

## 7. Reproduce

```bash
cd dual_certification_
python c1_audit_16.py
python validate_gap_numbers.py
python hybrid_scan_M.py --M 48,64,100 --jump-weight 2
python delta_aut_pairing.py --M 48 --compare --jump-weight 2
python lemma_K.py --test
```

---

## 8. Bottom line

§1.6 is **closed**: \(A_{\mathrm{bdry}}\) is λ-free; production C1 **1.98e3**, theorem gap **~3.9 orders** at η=5e-4. Hybrid/periodize remains the best residual path (**δ~0.09** at M=48, ~38× vs collocation), but **larger M does not improve absolute δ** with the present operator. Rung 4 still needs a true automorphic trial, not bigger Fourier truncations of a non-automorphic near-kernel. Hard map unchanged.

*Path: `dual_certification_/morningbrief.md`*
