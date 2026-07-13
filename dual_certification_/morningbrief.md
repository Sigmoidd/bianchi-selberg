# Morning brief — Rung 4 residual path

**Updated:** continued session — hybrid multi-pairing δ_aut, C1 §1.6 audit, Route A trunc scaffold  
**Hard map (do not change):**  
`F5_gluing=true, FEM_exclusion_01=true, two_cusp_hejhal_scan=true, defect_bound_applied=true,`  
`eta_le_eta0=false, width_lt_tol=false, counting_certified=false`  
⇒ **`pipeline_ok=true`, `rung4_certified=false`** ← still correct.

---

## 0. What you already have (do not re-litigate)

| GREEN | Note |
|-------|------|
| F5 + 6-copy gluing | exact combinatorics |
| FEM exclusion (0,1) only | Thm 2; upper near first eig = engineering |
| Two-cusp Hejhal **scan** | model residual, not Maass residual |
| Defect bound applied | honest huge interval when η large |
| Conditioning mitigated | \(b_{\mathrm{eq}}\approx 0.626 < 4\) at M=100…800; κ_eq ≈ 42…154 |
| C1 geometry 1.1+1.2 | **shipped** (see §1) |
| δ_aut hybrid multi + periodize | **shipped** (see §2) — **~8–38×** vs collocation |
| C1 §1.6 audit | **shipped** (see §1.6) — ~7× hypothetical only |
| Residual scan M≤800 | **shipped** (see §3) |
| Route A trunc scaffold | **shipped** (RED formulas, not cert) |

**Frozen scan** (`rung4_result.json`, model residual):  
\(r\approx 6.001\), \(\sigma_{\min}\approx 7\cdot 10^{-4}\), \(\delta\approx 1.4\cdot 10^{-4}\), \(\tau\approx 8.3\cdot 10^{-4}\), \(\eta\approx 9.6\cdot 10^{-4}\).  
Gap analysis uses **η = 5.04×10⁻⁴** as validation baseline.

**Independent validation** (`python validate_gap_numbers.py`):

| Regime | C1 | η₀ | η/η₀ @ 5e-4 | orders |
|--------|---:|---:|------------:|-------:|
| Conservative | \(1.02\cdot 10^{6}\) | \(2.41\cdot 10^{-13}\) | \(2.1\cdot 10^{9}\) | **9.3** |
| Sharp legacy (\(y_{\min}=1/2\), \(D_Y=3\)) | \(7.34\cdot 10^{4}\) | \(4.64\cdot 10^{-11}\) | \(1.1\cdot 10^{7}\) | **7.0** |
| **Sharp default (1.1+1.2)** | \(1.40\cdot 10^{4}\) | \(1.27\cdot 10^{-9}\) | \(4.0\cdot 10^{5}\) | **5.6** |
| Sharp default + \(U=2\) | \(1.40\cdot 10^{4}\) | \(5.08\cdot 10^{-9}\) | \(9.9\cdot 10^{4}\) | **5.0** |

Legacy quotes match within ≤1%. New default C1 is **5.24×** smaller than legacy sharp.

Targets from current default (C1≈1.4e4, η≈5e-4):

| Goal | Need roughly | Gap |
|------|--------------|-----|
| width &lt; 0.1 at Y0=1.5 | η ≲ **7e-6** | **~1.9 orders** |
| η ≤ η₀ (theorem, new default) | η ≲ **1.3e-9** | **~5.6 orders** |
| η ≤ η₀ (conservative) | η ≲ **2e-13** | **~9.3 orders** |

Primary blocker: still **η ≤ η₀** (production residual), not κ.

---

## 1. Reduce C1 further (geometry) — no circular η₀

**Formula (code):**  
\(C_1 = 2\,A_{\mathrm{bdry}}A_{\mathrm{res}}(1+A_{\mathrm{cusp}})\),  
\(A_{\mathrm{bdry}}=C_{\mathrm{tr}}C_{\mathrm{Sob}}A_{\mathrm{met}}(1+\sqrt\lambda)\),  
\(A_{\mathrm{res}}=(1+\lambda)(1+C_P A_{\mathrm{met}})\),  
\(C_{\mathrm{tr}}=3/h_{\min}\), \(C_K=1\).

**Default sharp geometry (now):**  
\(y_{\min}=1/\sqrt{2}\),  
\(D_Y=\sqrt{1^2+(1/2)^2+(Y-y_{\min})^2}\) ≈ 1.243 at Y=1.25.

| # | Formula change | Exp. factor on C1 | Status |
|---|----------------|------------------:|:------:|
| 1.0 | Sharp \(C_K=1\); \(C_{\mathrm{tr}}=3/h_{\min}\) | ~14× vs cons. | **Done** |
| 1.1 | \(y_{\min}:\ 1/2\to 1/\sqrt 2\) | **2.39×** | **Done** |
| 1.2 | Box \(D_Y\) (not crude 3) | **2.29×** alone; **5.24×** w/ 1.1 | **Done** |
| 1.3 | Mesh \(C_{\mathrm{tr}}\) / larger \(h_{\min}\) | ~1.2–4× | **Open** |
| 1.4 | Directional \(A_{\mathrm{met}}\) | ~1.5–2.5× | **Open** |
| 1.5 | Sobolev mean split | ~1.3–2× | **Open** |
| 1.6 | Drop \((1+\sqrt\lambda)\) / rewrite \(A_{\mathrm{ell}}\) | **7.08×** measured (`bdry_sqrt_lam=False`) | **Audited** — not production default |
| 1.7 | Stack remaining 1.3–1.5 + proved 1.6 | realistic **5–15×** more | **Open** |

**§1.6 audit result:** factor equals \(1+\sqrt{37}\) exactly on C1. Production keeps `bdry_sqrt_lam=True` until theorem_DK.tex is rewritten. See `C1_AUDIT_1_6.md`.

**Even after 1.1+1.2:** η₀ ~ 1e-9 vs η ~ 5e-4 → **5.6 orders** still open. Geometry alone does not green Rung 4.

**Circular η₀ ban:** never define \(y_{\min}/D_Y/h_{\min}\) from measured η; never invent \(U_{\mathrm{norm}}\).

---

## 2. Reduce δ (automorphy defect)

**Shipped:** `delta_aut_pairing.py` — true SL(2,ℂ) action, hybrid multi-pairing operator, optional Poincaré periodization.

**Measured (M=48, r=6, Y0=0.8, jump_weight=2):**

| mode | δ_aut | factor vs collocation |
|------|------:|----------------------:|
| collocation (baseline) | 3.39 | 1× |
| **multi (hybrid col+jumps)** | **0.443** | **7.6×** |
| **periodize (hybrid + Poincaré H)** | **0.088** | **38×** |

Per-gen (periodize): T1 ~0; R/TiR ~0.088; S ~0.011.

Still **O(10⁻¹)** vs target **10⁻⁶–10⁻⁹** — concrete progress, not green. Pure jump-only SVD was underdetermined (null-space garbage); hybrid is required.

| # | Change | Exp. factor | Status |
|---|--------|------------:|:------:|
| 2.0 | Honest language: model residual ≠ Maass residual | — | **Done** |
| 2.1 | PAIRINGS face jumps + hybrid multi operator | **7.6×** measured | **Done** (partial close) |
| 2.2 | True height pullback | hygiene | **Done** |
| 2.2b | Finite Poincaré on horizontal gens | **~38×** total vs col | **Done** (engineering) |
| 2.3 | Production \(S\) / σ₀ near true eigen-r | 10×–100× more if eigen | **Open** |
| 2.4 | Drop diagonal reg from Maass residual claims | removes ~1e-8 floor | **Open** |
| 2.5 | \(V_h^{P1,\mathrm{per}}\) conforming (not Neumann free) | correct trial space | **Open** |
| 2.6 | Companion \(J_h\) CR face-mean jumps | engineering | **Open** |

**Warning:** P1 Neumann with dropped pairing is **lower-bound flavor** only. Need \(V_h^{P1,\mathrm{per}}\subset H^1(\Gamma\backslash\mathbb{H}^3)\) or CR \(J_h\).

---

## 3. Reduce τ (L² residual) — M / n_pts / Arb / equilibration

**Shipped:** `residual_scan_M.py` + `residual_scan_M.csv` / `.md`.

**Measured (single-cusp model, r=6, Y0=0.8, amp+Sinkhorn):**

| M | n_pts | n_modes | rel | τ_svd | κ_eq |
|--:|------:|--------:|----:|------:|-----:|
| 100 | 484 | 316 | 2.38e-2 | 0.635 | 42 |
| 200 | 900 | 632 | 1.59e-2 | 0.597 | 63 |
| 400 | 1764 | 1256 | 1.03e-2 | 0.564 | 97 |
| 800 | 3364 | 2520 | 6.46e-3 | 0.509 | 155 |

log-log \(b(\mathrm{rel})\approx -0.63\). M:100→800 gives **~3.7×** on rel (not 10×).  
Path (3) alone: expect η proxy ~1e-2 → few×10⁻³ at M=1200; **not** 7 orders.

Note: `τ_disc=||Va||` for the exact SVD near-kernel is machine-tiny (by construction). Headline residual is **rel** / **τ_svd**.

| # | Change | Exp. factor | Status |
|---|--------|------------:|:------:|
| 3.0 | κ table M=100…800; n_pts=g² | — | **Done** |
| 3.1 | \(b_{\mathrm{eq}}\approx 0.626\); Dr V Dc | enables M↑ | **Done** |
| 3.2 | Residual scan M=100…800 (+ τ_svd) | **~3.7×** rel | **Done** |
| 3.2b | M=1000,1200 | ~1.5–2× more if slope holds | **Open** |
| 3.3 | Honest τ beyond SVD kernel | 1–2× | **Partial** |
| 3.4 | Arb 256-bit residual balls | cert enclosure | **Open** |
| 3.5 | Oversample collocation 2× | ~2–5× aliasing | **Open** |
| 3.6 | Y0↑ for collar / width | width win | **Partial** |

---

## 4. Relax η₀ via explicit ‖U‖₂

\[
\eta_0
=\min\Bigl(
\frac{1}{4(1+A_{\mathrm{bdry}})^2},\;
\frac{U_{\mathrm{norm}}^2}{4C_1^2}
\Bigr).
\]

| # | Change | Factor | Status |
|---|--------|-------:|:------:|
| 4.0 | A priori η₀ | — | **Done** |
| 4.1 | `U_norm` API | \(U=2\Rightarrow 4\times\) on 2nd term | **Done** |
| 4.2 | Measure certified \(\|U\|_2\) | only then U&gt;1 | **Open** |
| 4.3 | Paper sync if relaxation claimed | — | **Open** |

With new C1, U=2 gives η₀≈5e-9 — still **~5 orders** under η=5e-4. Do not invent U.

---

## 5. Master checklist (Done vs Open)

### Done
- [x] Hard map semantics; pipeline_ok true / rung4_certified false
- [x] FEM exclusion (0,1); F5 gluing; two-cusp model scan; Rung 3 Krawczyk
- [x] Conditioning \(b_{\mathrm{eq}}<4\); κ_eq 42–155 for M 100–800
- [x] Sharp \(C_K=1\); sharp_geom defaults
- [x] **1.1+1.2** y_min=1/√2 + box D_Y → C1 **1.40e4**, η₀ **1.27e-9**
- [x] η₀ non-circular + `U_norm` API
- [x] Gap docs + `validate_gap_numbers.py` (legacy + new default)
- [x] Overnight research swarm on (1)–(4)
- [x] **`delta_aut_pairing.py`** hybrid multi + periodize (**7.6× / 38×**)
- [x] **`residual_scan_M.py`** M=100…800 table
- [x] **§1.6 audit** `bdry_sqrt_lam` + `C1_AUDIT_1_6.md` (7.08× hypothetical)
- [x] **Route A** `truncation_constants_scaffold()` (formulas, still RED)

### Open — width &lt; 0.1 (η ≲ 7e-6 at new C1)
- [ ] Drive δ_aut from ~0.09 → ≲1e-5 (still ~3–4 orders)
- [ ] 3.2b M=1000–1200; 3.4–3.5 Arb / oversample
- [ ] Optional 1.3–1.5 + proved 1.6

### Open — η ≤ η₀ (~1e-9 new default)
- [ ] All of above + ~3–4 more residual orders beyond width target
- [ ] 4.2 real U_norm only if justified
- [ ] **Never** invent U_norm / C1 to green η_le_eta0

### Open — dual GREEN / rung4_certified
- [ ] eta_le_eta0 true
- [ ] width_lt_tol true (ε&lt;0.1)
- [ ] counting_certified (Route A/B)
- [ ] Milestone 1 reproducibility

---

## 6. Priority for next work

| Pri | Action | Why |
|----:|--------|-----|
| 1 | Larger M hybrid/periodize + true σ₀ coupling at eigen-r | δ still ~0.09 |
| 2 | Conforming \(V_h^{P1,\mathrm{per}}\) trial (not Neumann free) | correct space |
| 3 | Prove/reject 1.6 drop in theorem_DK.tex | 7× C1 if proved |
| 4 | Route A: Arb D/N + use trunc scaffold | counting_certified |
| 5 | Keep hard map; only flip when numbers say so | honesty |

---

## 7. Reproduce

```bash
cd dual_certification_
python validate_gap_numbers.py
python delta_aut_pairing.py --M 48 --compare --jump-weight 2
python c1_audit_16.py
python residual_scan_M.py --M 100,200,400,800
python lemma_K.py --test
python residue_F5_tests.py
# optional long:
# python residual_scan_M.py --M 1000,1200
# python rung4_N5_dual.py --M 32 --Y0 1.25 --Y0-defect 1.5
```

---

## 8. Bottom line

Pipeline still **runs**, does **not** certify. C1 geometry **5.24×** (orders 7.0→5.6); hybrid multi-pairing + periodize cuts δ_aut by **~8× / ~38×** (3.4 → 0.44 → 0.088) but remains **~5–7 orders** short of η₀. §1.6 can buy another **7×** on C1 only after a proof rewrite. Route A still RED on counting; trunc scaffold documents the missing constants. Hard map unchanged.

---

*Path: `dual_certification_/morningbrief.md`*
