# Rung 4 quantitative gap (honest stopping condition)

**Status:** `pipeline_ok = true`, `rung4_certified = false` — **correct** for current residuals.  
**Do not overclaim.** Keep the `hard` map as the stopping condition.

---

## What is GREEN

| Item | Value / evidence |
|------|------------------|
| F5 gluing | `true` — perms T1,R,TiR,S; cusp 1+5; round-trips |
| FEM exclusion (0,1) | `true` — Theorem 2; mesh 8×4×3 24-split; n≈28400; Y=1.25; ν\*=1.02; float pencil min-eig ≈1.44; NP=5 index 6 |
| FEM language | *Only* exclusion of (0,1). Upper bound near first eig = **engineering target** |
| Two-cusp Hejhal scan | `true` — Rung 3 operator |

Typical scan (model residual, not production Maass residual):

```
r ≈ 6.0
sigma_min ≈ 1e-4
delta ≈ 1e-4
tau ≈ 4e-4
eta ≈ 5e-4
```

---

## Why `eta_le_eta0 = false` blocks dual width

Defect theorem requires \(\eta\le\eta_0\) for the full L² proximity claim. Numerically:

### Before sharp geometry (conservative \(D_Y=10\), \(h_{\min}=0.25\), \(C_{\mathrm{tr}}=16\))

```
C1 ≈ 1.02e6
C2 ≈ 11.7
eta0 = min(1/4, (2 C1)^{-2}) ≈ 2.4e-13
eta  ≈ 5e-4
eta/eta0 ≈ 2e9   (~9 orders)
```

Then

```
dλ = C1·η + C2 e^{-2π Y0} ≈ 515
hejhal λ-interval ≈ [37−515, 37+515] = [−478, 552]
overlap with [1,∞) width ≈ 550 ≫ 0.1
```

### After sharp geometry (default now: \(C_K=1\), \(D_Y=3\), \(C_{\mathrm{tr}}=3/h_{\min}\approx 6.67\))

```
C1 ≈ 7.3e4     (~14× smaller)
C2 ≈ 6.7
eta0 ≈ 4.6e-11 (~200× larger)
eta  ≈ 5e-4
eta/eta0 ≈ 1e7   (~7 orders)   ← still blocked, better constants
```

With \(U_{\mathrm{norm}}=2\): \(\eta_0\sim 1.9\cdot 10^{-10}\) (further ×4). Still \(\eta\gg\eta_0\).

---

## Dominant uncertainty (shifted)

| Was | Now |
|-----|-----|
| Conditioning \(\kappa\sim M^{76}\) | **Mitigated** float: \(b_{\mathrm{eq}}\approx 0.63<4\) at tested params |
| Primary blocker | **Mathematical:** \(\eta\le\eta_0\) (defect residual), then counting |

Conditioning risk is **reduced empirically**, not removed as a theorem. Reproducibility across \(Y_0,r\), precision, levels still required.

---

## Immediate fix plan (tracked)

| # | Fix | Effect | Status |
|---|-----|--------|--------|
| 1 | Sharp \(C_K=1\), \(D_Y=3\), \(C_{\mathrm{tr}}=3/h_{\min}\), \(y_{\min}\) metric | \(C_1\downarrow\sim 10^2\), \(\eta_0\uparrow\sim 10^2\) | **Done** in `C1_C2_constants(sharp_geom=True)` |
| 2 | Relax \(\eta_0\) via explicit \(\|U\|_2\) | \(\eta_0 \propto U_{\mathrm{norm}}^2\) | **Done** (`U_norm` param) |
| 3 | Reduce \(\eta\) by 5–9 orders (M↑, Arb pairing, better \(\delta,\tau\)) | Primary remaining work | **Open** — model residual \(\sim 10^{-4}\) |
| 4 | Route A counting cert on \((1,\lambda_\star)\) | Makes width meaningful | **Open** — prototype only |

---

## Hard map (do not change semantics)

```
F5_gluing: true
FEM_exclusion_01: true
two_cusp_hejhal_scan: true
defect_bound_applied: true
eta_le_eta0: false          ← blocks dual GREEN
overlap_nonempty_disjoint_unit: false  (wide interval crosses nonsense)
width_lt_tol: false
counting_certified: false
⇒ pipeline_ok: true
⇒ rung4_certified: false
```

---

## What “GREEN dual” would require (numbers)

Need roughly

\[
C_1\eta + C_2 e^{-2\pi Y_0} < 0.1.
\]

At sharp \(C_1\sim 7\cdot 10^4\), \(Y_0=1.5\) (\(C_2 e^{-2\pi Y_0}\sim 10^{-3}\)):

\[
\eta \lesssim \frac{0.1}{C_1} \sim 10^{-6}
\]

for a **width** target alone; for \(\eta\le\eta_0\sim 10^{-10}\)–\(10^{-11}\) (theorem hypothesis as written) need **~6–7 more orders** on \((\delta,\tau)\) beyond current \(10^{-4}\).

Synthetic demo (`eta=1e-14`, `Y0=1.5`) already shows arithmetic path to \(\mathrm{d}\lambda\sim 10^{-3}<0.1\).

---

## Reproduce

```bash
cd dual_certification_
python -c "from lemma_K import C1_C2_constants, ball_mid; d=C1_C2_constants(r=6,Y0=0.8,sharp_geom=True); print(ball_mid(d['C1']), ball_mid(d['eta0']))"
python rung4_N5_dual.py --M 32 --Y0 1.25 --Y0-defect 1.5
```
