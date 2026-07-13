# Rung 4 certificate — Γ₀(2+i), N𝔭 = 5 dual pipeline

**Source of truth:** `DualCertification_AgentReady.md` §6 Rung 4 + §13  
**Date:** 2026-07-12  
**Runner:** `python rung4_N5_dual.py`

---

## AgentReady checklist

```
[x] FEM lower bound already certified (μ≈4.4 pencil diagnostic; Thm 2:
    no eigenvalue in (0,1)) is re-used / gluing re-verified.
      → independent_exclusion/PROOF.md Theorem 2; fem_N5_status()
      → frozen mesh 8×4×3, n≈28400, Rump 8/8 (2026-07-10)

[x] Two-cusp Hejhal uses the coupling analysis of Rung 3.
      → two_cusp_hejhal_N5.build_block_system + residual scan in rung4_N5_dual.py
      → σ₀ block form from two_cusp_coupling.md / Rung 3

[ ] Overlap + counting succeed with ε<0.1.
      → OVERLAP: FEM [1,∞) ∩ Hejhal+D(K) interval computed (status YELLOW/RED
        depending on residual size — see results).
      → COUNTING: YELLOW — Route A prototype only (route_A_status.md).
      → Full GREEN blocked until (i) Hejhal residual η ≲ 1e-8–1e-10 with
        production automorphy defect, (ii) certified counting N(λ)=0 on (1,λ₁).

[x] Residue-ring arithmetic for F5 verified by independent unit-test suite.
      → residue_F5_tests.py  OVERALL PASS

[x] Gluing of the 6 copies is exact (combinatorial identity, not numerical).
      → residue_F5_tests + congruence_prototype.build_gluing('(2+i)')
      → perms bijective; cusp_class (1 ∞ + 5 finite); round-trips OK
```

---

## Stopping condition (AgentReady)

> Certified interval \([\lambda_1-\varepsilon,\lambda_1+\varepsilon]\) with \(\varepsilon<0.1\),  
> disjoint from \((0,1)\), counting certifies it is the first eigenvalue.

| Requirement | Status |
|-------------|--------|
| Disjoint from (0,1) | FEM certifies empty (0,1); Hejhal interval must sit above 1 |
| \(\varepsilon<0.1\) | Needs \(C_1\eta + C_2 e^{-2\pi Y_0}<0.1\) |
| Counting first | **Open** (Route A/B not certified) |
| **Rung 4 §13 CERTIFIED** | **NO** — **correct** for current residuals |

**`pipeline_ok = true`, `rung4_certified = false` is the right output.**

### Quantitative gap (see also `RUNG4_GAP.md`)

| | Conservative geom | **Sharp geom (default now)** |
|--|------------------:|-----------------------------:|
| \(C_1\) (r=6, Y0=0.8) | \(\sim 1.02\cdot 10^{6}\) | \(\sim 7.3\cdot 10^{4}\) |
| \(\eta_0\) | \(\sim 2.4\cdot 10^{-13}\) | \(\sim 4.6\cdot 10^{-11}\) |
| Current \(\eta\) | \(\sim 5\cdot 10^{-4}\) | same |
| \(\eta/\eta_0\) | \(\sim 2\cdot 10^{9}\) (~9 orders) | \(\sim 10^{7}\) (~7 orders) |

With \(\eta\approx 5\cdot 10^{-4}\), \(C_1\sim 10^{6}\) (old): \(\mathrm{d}\lambda\sim 515\), Hejhal interval \([-478,552]\), overlap width \(\sim 550\gg 0.1\).

**Primary blocker shifted:** conditioning (\(b_{\mathrm{eq}}\approx 0.63\)) is no longer first; **\(\eta\le\eta_0\)** (mathematical residual certification) is. Still need \(\delta,\tau\sim 10^{-10}\)–\(10^{-13}\) even after sharp \(C_1\).

Hard map (do not soften):
```
F5_gluing true, FEM_exclusion_01 true, two_cusp_hejhal_scan true,
defect_bound_applied true, eta_le_eta0 false,
overlap_nonempty_disjoint_unit false, width_lt_tol false,
counting_certified false  ⇒  rung4_certified false
```


---

## What is GREEN now

| Piece | Evidence |
|-------|----------|
| F₅ arithmetic + 6-copy gluing | `python residue_F5_tests.py` → PASS |
| FEM exclusion (0,1) | Theorem 2 provenance + glue re-verify |
| Two-cusp residual scan | Rung 3 operator, r-grid + golden refine |
| D(K) plug-in | `defect_to_lambda_error` |
| Dual overlap arithmetic | FEM half-line ∩ Hejhal interval |
| Synthetic tight-η demo | Shows width&lt;0.1 path when η=1e-10, Y0=1.5 |

---

## Blockers to full §13 certification

1. **Production Hejhal residual.** Discrete collocation relative residual is still \(\eta\sim 10^{-3}\)–\(10^{-4}\) (model operator does not yet vanish at a true eigenvalue). With \(C_1\sim 10^6\), need \(\eta\lesssim 10^{-8}\) for \(\mathrm{d}\lambda<0.1\) even at large \(Y_0\).
2. **Counting certificate.** Route A/B must return \(N(\lambda)=[0,0]\) on \((1,\lambda_1-\varepsilon)\).
3. **Reproducibility suite (Milestone 1)** — not run for a dual eigenvalue interval.
4. **Optional:** sharpen \(C_1\) in Theorem D(K) to reduce residual demand.

---

## Reproduce

```bash
cd dual_certification_
python residue_F5_tests.py
python rung4_N5_dual.py --M 32 --Y0 1.25 --r-grid 6.0:14.0:21
python two_cusp_hejhal_N5.py --M 48          # Rung 3 Krawczyk
```

---

## Hard items (last run)

See `rung4_result.json` for numbers. Typical:

| Item | Typical |
|------|---------|
| F5_gluing | YES |
| FEM_exclusion_01 | YES |
| two_cusp_hejhal_scan | YES |
| defect_bound_applied | YES |
| eta_le_eta0 | NO (until residual drops) |
| width_lt_tol | NO |
| counting_certified | NO |
| **pipeline_ok** | **YES** |
| **rung4_certified** | **NO** |

---

## Language (AgentReady §12)

- FEM exclusion of (0,1): **certified** (Theorem 2).  
- First eigenvalue near ~45: **engineering target** until Hejhal residual + counting close.  
- Dual loop: **addresses** the logical gap; **not** a full §13 certificate in this folder yet.  
- Then \(N=5\) is the first congruence level where TF \(B\gtrsim 1\); dual certification is the proposed route.

---

## Next actions (ordered)

1. Production single-/two-cusp Hejhal iterate driving \(\delta_{\mathrm{aut}},\tau_{\mathrm{tail}}\) below \(10^{-8}\).  
2. Route A/B certified counting on \((1,\lambda_\star)\).  
3. Re-run `rung4_N5_dual.py` until `rung4_certified=true`.  
4. Milestone 1 reproducibility on the final interval.
