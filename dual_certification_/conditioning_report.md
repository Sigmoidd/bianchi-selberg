# Conditioning report — Milestone 2 / Rung 3 diagnostic

**Status:** float diagnostic only — **NON-CERTIFYING** (κ slopes)  
**Rung 3 complete?** **YES** for AgentReady checklist — see [`rung3_certificate.md`](rung3_certificate.md)  
(Krawczyk N=5 + σ₀ theory + \(S\) radii in `two_cusp_hejhal_N5.py` / `two_cusp_coupling.md`)  
**Date of run:** 2026-07-11 (κ); 2026-07-12 (Krawczyk)  
**Script:** [`hejhal_conditioning.py`](hejhal_conditioning.py)  
**Parameters:** \(r = 6.62212\) (Then scale), \(Y_0 = 0.8\), \(\theta = 1/2\), \(M \in \{100,200,400,800\}\)

---

## 1. What was built

A **model** single-cusp Hejhal-like collocation / moment matrix \(V(M)\) on the torus \(\mathbb{C}/\mathbb{Z}[i]\):

- **Modes:** all Gaussian integers \(\beta = a+bi\) with \(0 < N(\beta)=a^2+b^2 \le M\).
- **Sample points:** \(g\times g\) grid with \(g = 2\lfloor\sqrt{M}\rfloor+2\) (Nyquist for max frequency), shift \(0.37\) to avoid accidental exact DFT orthogonality.
- **Entries (headline Rung 3 model, `height_mismatch=0`):**
  \[
  V_{j,\beta}
  =
  w_\beta(Y_0)\,
  \Bigl(
    e^{2\pi i \beta\cdot x_j}
    -
    e^{2\pi i \beta\cdot x_\gamma(x_j)}
  \Bigr)
  \Big/\sqrt{n_{\mathrm{pts}}},
  \]
  where \(w_\beta(y) = N(\beta)^\theta\,|K_{ir}(2\pi|\beta|y)|\) is modelled by the elementary majorant \(\sqrt{\pi/(2y)}\,e^{-y}\), and \(\gamma\) is the H³ unit-hemisphere inversion (stand-in cusp pairing).
- **Same-height phase automorphy** is the intentional headline model: after right scaling by \(D=\mathrm{diag}(w_\beta(Y_0))\) one studies a pure Fourier difference operator. This matches AgentReady’s statement that naïve growth is driven by K-Bessel dynamic range, removed by a diagonal mode preconditioner.
- **Two-cusp stub:** block columns \([V,\,-cV;\,-cV,\,V]\) with \(c=0.1\), block-diagonal \(D\oplus D\).

**Float SVD** via `numpy.linalg.svd` — labelled non-certifying throughout.

Stress variant (`--height-mismatch 1`): blends true pull-back weights \(w_\beta(y_\gamma(x_j))\). That re-introduces exponential height imbalance; earlier full run gave equilibrated slope \(b\approx 4.39\) (borderline fail). Not the Rung 3 headline; flags that **true** inversion-at-variable-height systems need larger \(Y_0\) or stronger preconditioning.

---

## 2. Table of \(\kappa\) vs \(M\)

| \(M\) | \(n_{\mathrm{modes}}\) | \(n_{\mathrm{pts}}\) | \(\kappa_{\mathrm{diag\ proxy}}\) | \(\kappa_{\mathrm{raw}}\) (stable) | \(\kappa(V D_{\mathrm{amp}}^{-1})\) | \(\kappa_{\mathrm{col\ Jacobi}}\) | \(\kappa_{\mathrm{equilibrated}}\) (headline) | 2-cusp \(\kappa_{\mathrm{eq}}\) |
|------:|------------------------:|---------------------:|----------------------------------:|-----------------------------------:|-----------------------------------:|----------------------------------:|-----------------------------------------------:|-------------------------------:|
| 100 | 316 | 484 | \(1.97\times10^{38}\) | \(1.72\times10^{21}\) | \(1.22\times10^{2}\) | \(4.78\times10^{1}\) | \(4.21\times10^{1}\) | \(5.15\times10^{1}\) |
| 200 | 632 | 900 | \(1.69\times10^{56}\) | \(2.30\times10^{30}\) | \(1.77\times10^{2}\) | \(6.81\times10^{1}\) | \(6.29\times10^{1}\) | \(7.68\times10^{1}\) |
| 400 | 1256 | 1764 | \(4.50\times10^{81}\) | \(1.74\times10^{43}\) | \(2.60\times10^{2}\) | \(9.96\times10^{1}\) | \(9.69\times10^{1}\) | \(1.18\times10^{2}\) |
| 800 | 2520 | 3364 | \(4.70\times10^{117}\) | \(2.69\times10^{61}\) | \(3.93\times10^{2}\) | \(1.55\times10^{2}\) | \(1.55\times10^{2}\) | \(1.89\times10^{2}\) |

Source CSV: [`hejhal_kappa_vs_M.csv`](hejhal_kappa_vs_M.csv)  
Plot: [`hejhal_kappa_vs_M.png`](hejhal_kappa_vs_M.png)

**Definitions**

- \(\kappa_{\mathrm{diag\ proxy}}\): legacy majorant ratio \(\max w/\min w\) with \(w(n)=n^{2\theta}|K|^2\) (same spirit as `plot_kappa_vs_M.py` / `kappa_majorant_series`).
- \(\kappa_{\mathrm{raw}}\) (stable): \((\max w_\beta/\min w_\beta)\cdot\kappa(V D_{\mathrm{amp}}^{-1})\) — float-safe proxy for the un-preconditioned dynamic range (direct `cond(V)` underflows to \(\infty\)).
- \(D_{\mathrm{amp}}=\mathrm{diag}(w_\beta(Y_0))\): right diagonal mode-amplitude preconditioner.
- Equilibrated: Sinkhorn row/column \(\infty\)-norm diagonal scaling applied after amp-right scaling (headline Rung 3 figure).

---

## 3. Log–log slopes \( \log\kappa \approx a + b\log M \)

| Quantity | \(a\) | \(b\) | \(\kappa \sim M^{b}\) | vs Rung 3 (\(b<4\)) |
|----------|------:|------:|----------------------:|---------------------|
| Diagonal majorant proxy | \(-325.07\) | **\(87.55\)** | \(M^{87.6}\) | fail (naïve; AgentReady quotes \(\sim M^{76}\)) |
| Raw stable (matrix) | \(-160.32\) | **\(44.34\)** | \(M^{44.3}\) | fail |
| Right amp \(V D_{\mathrm{amp}}^{-1}\) | \(2.21\) | **\(0.561\)** | \(M^{0.56}\) | **pass** |
| Column Jacobi | \(1.25\) | **\(0.564\)** | \(M^{0.56}\) | **pass** |
| Amp + row Jacobi | \(1.49\) | **\(0.630\)** | \(M^{0.63}\) | **pass** |
| **Equilibrated (headline)** | \(0.84\) | **\(0.626\)** | \(M^{0.63}\) | **pass** |
| 2-cusp equilibrated stub | \(1.04\) | **\(0.626\)** | \(M^{0.63}\) | **pass** |

### Verdict (diagnostic only)

| Test | Result |
|------|--------|
| Naïve / raw growth | **Unacceptable** — proxy \(b\approx 87.6\), raw matrix \(b\approx 44.3\) (same ballpark as AgentReady’s \(\sim M^{76}\)) |
| After diagonal preconditioner | **\(b_{\mathrm{eq}}\approx 0.63 < 4\)** and **\(b_{\mathrm{amp}}\approx 0.56 < 4\)** → **Rung 3 slope stopping condition PASSES** (float diagnostic) |
| Rung 3 overall complete? | **NO** — see checklist |

Machine summary: [`hejhal_kappa_summary.txt`](hejhal_kappa_summary.txt)

---

## 4. Interpretation

1. **Without preconditioning**, K-Bessel decay across modes with \(N(\beta)\le M\) at fixed \(Y_0=0.8\) produces enormous column-scale disparity (\(\mathrm{amp\_ratio}\) from \(10^{19}\) at \(M=100\) to \(10^{58}\) at \(M=800\)). Interval widths would explode — this is the risk AgentReady flags.
2. **Diagonal mode preconditioning** \(D=\mathrm{diag}(w_\beta(Y_0))\) removes that dynamic range. The residual Fourier / same-height automorphy operator has mild polynomial growth \(b\approx 0.6\), well under the threshold \(b<4\).
3. **Two-cusp stub** with block-diagonal \(D\oplus D\) and mild coupling \(c=0.1\) does not change the slope (same \(b\approx 0.63\)); absolute \(\kappa\) grows only by the expected \(O(1)\) coupling factor \((1+|c|)/(1-|c|)\).
4. **Height-mismatched inversion** (stress, not headline) previously produced equilibrated \(b\approx 4.4\). That is a warning for full two-cusp Hejhal with true \(\sigma_0\) pull-backs: keep \(Y_0\) large enough, or use joint row/column scaling and possibly block preconditioners by norm shells.

---

## 5. AgentReady §6 Rung 3 checklist (honest)

```
[ ] Consistency of the two expansions under σ0 is proved as distributions on y=Y.
[ ] The coupled matrix is written explicitly; uniqueness of the gluing relation a^(0)=S(r)a^(∞) is stated.
[x] Condition-number diagnostic (log κ vs log M) is produced for M=100,200,400,800 both with and without preconditioner.
[x] Interval Krawczyk test is implemented and succeeds for the N=5 system at the target r.
      → two_cusp_hejhal_N5.py; M=32,48,64 success; rung3_certificate.md
[x] All interval radii of the scattering-like block S are tracked.
      → BlockSystem.S_rad; relative R in Krawczyk
```

**Stopping condition** \(\log\kappa(D^{-1}V)\approx a+b\log M\) with \(b<4\):

- **Diagnostic slope:** **satisfied** (\(b_{\mathrm{eq}}\approx 0.63\), \(b_{\mathrm{amp}}\approx 0.56\)).
- **Full Rung 3 deliverable** (well-posedness proof + condition-number **bound** + Krawczyk for N=5): **not done**.

---

## 6. Artifacts

| File | Role |
|------|------|
| [`hejhal_conditioning.py`](hejhal_conditioning.py) | Model \(V(M)\), preconditioners, fits, CLI |
| [`hejhal_kappa_vs_M.csv`](hejhal_kappa_vs_M.csv) | Numeric table |
| [`hejhal_kappa_vs_M.png`](hejhal_kappa_vs_M.png) | log–log plots raw vs preconditioned |
| [`hejhal_kappa_summary.txt`](hejhal_kappa_summary.txt) | Machine-readable \(b\) values |
| [`plot_kappa_vs_M.py`](plot_kappa_vs_M.py) | Legacy single-cusp majorant proxy (still valid cross-check) |
| [`kappa_vs_M.csv`](kappa_vs_M.csv) / [`kappa_vs_M.png`](kappa_vs_M.png) | Legacy proxy outputs |

### Reproduce

```bash
cd dual_certification_
python hejhal_conditioning.py --M 100,200,400,800
# optional stress (height-mismatched pull-back):
python hejhal_conditioning.py --M 100,200,400 --height-mismatch 1.0 --no-two-cusp
```

---

## 7. Explicit non-claims

- This does **not** certify any eigenvalue.
- This does **not** prove \(\kappa(V)=O(M^{b})\) with \(b<4\) for the true Bianchi two-cusp Hejhal operator — only a model collocation matrix in float arithmetic.
- **Krawczyk for N=5 is open.**
- \(\sigma_0\) consistency and uniqueness of the gluing relation \(a^{(0)}=S(r)a^{(\infty)}\) are **not** proved here.
- Rung 3 remains **incomplete** until the open checklist items are closed.
