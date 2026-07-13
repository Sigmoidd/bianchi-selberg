# Dual certification — Rung 0–1: Lemma K and Theorem D(K)

Self-contained analytic core for single-cusp Bianchi defect bounds.
**No gluing, no \(N\mathfrak{p}\), no reference cell, no two-cusp coupling in the main theorems.**

Revision 2026-07-12: full elementary \(C_K\) proof, consistent \(C_{\mathrm{tr}}\), metric comparison, non-circular \(\eta_0\), named \(C_1/C_2\). See [`REVISION_NOTES.md`](REVISION_NOTES.md) and [`REFEREE_REPORT.md`](REFEREE_REPORT.md).

## Contents

| File | Role |
|------|------|
| [`theorem_DK.tex`](theorem_DK.tex) | Full paper (amsart): Assumption H, Lemma K, Theorem D(K), explicit \(C_1,C_2\) |
| [`lemma_K.py`](lemma_K.py) | Arb-enclosable \(C_K(r)\), tail majorant \(S_{M,\infty}\), `C1_C2_constants` |
| [`plot_kappa_vs_M.py`](plot_kappa_vs_M.py) | Diagnostic: single-cusp \(\kappa_{\mathrm{proxy}}\) vs \(M\) → CSV/PNG |
| [`constants_DK.md`](constants_DK.md) | Numerical table of geometric constants, \(C_K\), \(C_1\), \(C_2\) |
| [`REVISION_NOTES.md`](REVISION_NOTES.md) | Referee blocking items → equation-level fixes |
| [`DualCertificationroadmap.md`](DualCertificationroadmap.md) | Architecture roadmap (context) |
| [`annotated_bibliography_dual.md`](annotated_bibliography_dual.md) | Literature map |

## What the paper claims

Under **Assumption H\((K,\theta)\)** (Hecke growth \(|a_\beta|\le C_H(\varepsilon)\,N(\beta)^{\theta+\varepsilon}\)):

1. **Lemma K.** Elementary pointwise **upper** bound
   \[
   |K_{ir}(y)|\le\sqrt{\frac{\pi}{2y}}\,e^{-y}\,C_K(r),\quad
   C_K(r)=1
   \]
   (proved from the integral representation via \(K_0\le K_{1/2}\); no \(r\)-independent lower bound), and a fully explicit, Arb-computable bound on the Fourier tail
   \(S_{M,\infty}=\sum_{N(\beta)>M}|a_\beta|^2|K_{ir}(2\pi|\beta|Y_0)|^2\)
   under H and \(r_2(n)\le 6d(n)\).

2. **Theorem D(K).** A \(K\)-finite truncated expansion with automorphy defect
   \(\delta_{\mathrm{aut}}\) and residual \(\tau_{\mathrm{tail}}\)
   (\(\eta=\delta_{\mathrm{aut}}+\tau_{\mathrm{tail}}\) small) is within
   \[
   |\lambda_{\mathrm{true}}-\lambda|\le C_1\eta+C_2 e^{-2\pi Y_0},\qquad
   \|U_{\mathrm{cusp}}-f\|_2\le C_1\sqrt{\eta}
   \]
   of a true cusp form (cuspidal projection), with \(C_1,C_2\) assembled from
   **named** contributions:
   \[
   C_1=2\,A_{\mathrm{bdry}}A_{\mathrm{res}}(1+A_{\mathrm{cusp}}),\qquad
   C_2=A_{\mathrm{cut}},
   \]
   using trace (\(C_{\mathrm{tr}}=\frac43|F|/|T|\)), metric comparison on \(\{y\ge y_{\min}\}\), Payne–Weinberger, Sobolev, and Lemma K. The threshold \(\eta_0\) is defined from a priori constants only.

## What the paper does **not** claim

- Two-cusp coupling or the block matrix \(\begin{pmatrix}V_\infty&-S\\-S^*&V_0\end{pmatrix}\)
- Congruence level \(N\mathfrak{p}\), reference-cell gluing, FEM certificates
- Certified counting \(N(\lambda)\); Weyl law as an exact count
- Uniform gap \(\mu(N)\) (heuristic appendix only)
- That Then’s \(r_1\approx 6.62212\) is a theorem (it is an **engineering target** until a concrete \(\eta\) is small enough)

## Requirements

- Python 3.10+
- **Preferred:** [`python-flint`](https://github.com/flintlib/python-flint) (`from flint import arb`) for certifying Arb balls
- Optional: `matplotlib` for `kappa_vs_M.png`; `scipy` only if flint is absent (float incomplete-gamma fallback)

This repository environment has python-flint available.

## How to run

From this directory:

```bash
# Self-tests + M=100,200,400 majorant benchmark + C1/C2 table
python lemma_K.py --test --bench --constants

# Single evaluation
python lemma_K.py --M 400 --Y0 0.8 --r 6.6 --theta 0.5

# Classical (looser) C_K = exp(π r/2) comparison
python lemma_K.py --bench --classical

# Condition-number proxy diagnostic (writes kappa_vs_M.csv and .png)
python plot_kappa_vs_M.py
```

### API (`lemma_K.py`)

```python
from lemma_K import C_K, lemma_K_tail, C1_C2_constants, HAS_FLINT

ck = C_K(6.6)                    # sharp = 1 (Arb if HAS_FLINT)
ck_cl = C_K(6.6, classical=True) # exp(π r/2)
bound = lemma_K_tail(            # enclosure of analytic majorant for S_{M,∞}
    M=400, Y0=0.8, r=6.6, theta=0.5, C_H=1.0, eps=0.0
)
consts = C1_C2_constants(field="i")   # C1, C2, A_bdry, ...
# consts = C1_C2_constants(field="omega")
```

- If `HAS_FLINT` is true, return values are **certifying** Arb enclosures of the *majorant / constant formulas*.
- If flint is missing, floats are returned and must be treated as **non-certifying diagnostics**.

### Paper build

```bash
pdflatex theorem_DK.tex
pdflatex theorem_DK.tex   # for TOC
```

## Hecke input (where it enters)

Assumption H is stated once (paper §3). It enters **only** in:

- the factors \(C_H(\varepsilon)^2\) and the power \(2\theta+2\varepsilon\) inside Lemma K;
- the cusp piece \(A_{\mathrm{cusp}}\) inside \(C_1\).

Trace, Poincaré, Sobolev, metric comparison, and sharp \(C_K\equiv 1\) are pure analysis / geometry.

## Success criteria (referee checklist)

1. **Lemma K tail** is checkable from the elementary pointwise **upper** \(K_{ir}\) bound + Assumption H alone (`lemma_K.py` + paper §4). No false Luke lower bound.
2. **Theorem D(K) constant \(C_1\)** is Arb-computable from named factors: \(C_{\mathrm{trace}}\), \(C_{\mathrm{Poincaré}}\), \(C_{\mathrm{met}}\), Lemma K (`C1_C2_constants` + `constants_DK.md` + paper §5). Every factor maps to a lemma; no hidden \(N\mathfrak{p}\)/gluing.
3. **\(\eta_0\)** is defined from a priori \(A_{\mathrm{bdry}},C_1\) only (no circular justification).
4. Separation known / engineering / unknown preserved; Then \(r_1\) remains engineering.

### Reference numerics at \((Y,Y_0,r)=(1.25,0.8,6.6)\), sharp \(C_K=1\)

| | \(\mathbb{Z}[i]\) | \(\mathbb{Z}[\omega]\) |
|--|------------------:|----------------------:|
| \(C_1\) | \(\approx 1.324\cdot 10^{6}\) | \(\approx 3.564\cdot 10^{5}\) |
| \(C_2\) | \(\approx 11.728\) | \(\approx 3.346\) |
| Lemma K tail \(M=400\) | \(\sim 2.8\cdot 10^{-83}\) | same |

## Relation to the roadmap

| Rung | Status in this folder |
|------|------------------------|
| 0 — Lemma K | **Done** — `lemma_K_certificate.md` |
| 1 — Theorem D(K) | **Done** — `rung1_certificate.md`, `defect_bound_arb.py` |
| 2 — Dual overlap level 1 | Partial (FEM half-line + Hejhal when η small) |
| 3 — Two-cusp + Krawczyk N=5 | **Done** — `rung3_certificate.md` |
| 4 — N=5 dual pipeline | **Infra done** — `rung4_certificate.md`; §13 cert open |
| 5+ — scaling / field port | Open |

```bash
python residue_F5_tests.py
python two_cusp_hejhal_N5.py --sweep-M 32,48,64
python rung4_N5_dual.py --M 32 --Y0 1.25 --Y0-defect 1.5
python hejhal_conditioning.py --M 100,200,400,800
```
