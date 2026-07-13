# Referee report — Round 2 (post-revision)

**Date:** 2026-07-12  
**Scope:** Revisions along Luke → Lemma K → tail → trace → residual → Poincaré → D(K)  
**Evidence:** `REVISION_NOTES.md`, re-run of `python lemma_K.py --test --bench --constants`

---

## Verdict

### **CONDITIONAL ACCEPT** (analysis core for Rung 0–1)

The revision **closes all blocking items** from Round 1 that can be closed on paper + Arb majorant code. Lemma K is now an elementary, checkable upper-bound argument. Theorem D\((K)\) has named constants, non-circular \(\eta_0\), metric comparison, and no double-counted \(C_K\).

Remaining weaknesses are **engineering sharpness** (conservative \(C_{\mathrm{trace}}=16\), \(D_Y=10\), soft elliptic scaling) rather than structural holes. A hard journal referee may still ask for one more page of chart-by-chart extension estimates in Step 1; they no longer block the dual-certification dependency chain.

| Success criterion | Round 1 | Round 2 |
|-------------------|---------|---------|
| Lemma K from upper \(K\)-bound + H | Partial | **PASS** |
| \(C_1\) from named lemmas, no magic 4 | Partial | **PASS** (factor 2 derived) |
| No \(N\mathfrak{p}\)/gluing | PASS | **PASS** |
| Code Arb + M=100,200,400 | PASS | **PASS** |
| No false Luke lower bound | FAIL | **PASS** |
| Non-circular \(\eta_0\) | FAIL | **PASS** |
| Hyperbolic vs Euclidean | Missing | **PASS** (Lemma metric) |

---

## Chain check

```
Luke / asymptotic UPPER only  ──►  PASS (eq:Luke; no r-free lower bound)
        │
Lemma K pointwise (A)(B)(C)   ──►  PASS  |K_ir|≤K_0≤K_{1/2}=√(π/2y) e^{-y}, C_K=1
        │
Tail under H                  ──►  PASS  r₂≤6d(n); incomplete gamma; Arb
        │
Trace C_tr = (4/3)|F|/|T|     ──►  PASS  single constant
        │
Residual + metric y_min       ──►  PASS  named A_bdry, A_ell, A_res, A_met
        │
Poincaré D_Y/π                ──►  PASS
        │
Theorem D(K)                  ──►  PASS structure; C1=2 A_bdry A_res (1+A_cusp)
```

---

## Numerical (referee re-run)

| Quantity | \(\mathbb{Z}[i]\) | \(\mathbb{Z}[\omega]\) |
|----------|------------------:|----------------------:|
| \(C_K\) sharp | 1 | 1 |
| \(C_1\) | \(\approx 1.324\cdot 10^{6}\) | \(\approx 3.564\cdot 10^{5}\) |
| \(C_2\) | \(\approx 11.73\) | \(\approx 3.35\) |
| Tail \(M=400\) | \(\sim 2.84\cdot 10^{-83}\) | same |
| \(\eta_0\) | \(\sim 1.4\cdot 10^{-13}\) | \(\sim 2.0\cdot 10^{-12}\) |

Improvement vs Round 1: \(C_1\) dropped \(\sim 10^{3}\) by using sharp \(C_K=1\) and removing double-counted \((1+C_K)\).

---

## Minor non-blocking notes

1. **Monotonicity of \(K_\nu\) in \(\nu\)** (Step B) is cited to DLMF §10.37 — acceptable; an optional one-line reference to the integral \(\partial_\nu K_\nu\) sign is enough if a pedantic referee objects.
2. **Step 1 extension operator** is still somewhat schematic (tubular charts); constants are conservative. Fine for architecture; tighten before journal submission.
3. **\(A_{\mathrm{ell}}=1+\lambda\)** is a scaling choice, not a sharp elliptic estimate — label as such (paper does via named constants).
4. Optional: numerical check `|scipy.kve(1j*r, y)|` vs majorant at a few points (diagnostic only).

---

## Recommendation

**Accept into the dual-certification stack as Rung 0–1.**  
Next work: two-cusp block matrix on top of D(K); production of a putative form with certified \(\eta\le\eta_0\) near Then’s \(r_1\).
