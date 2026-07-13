# Revision notes — referee blocking items → fixes

**Manuscript:** `theorem_DK.tex`  
**Code:** `lemma_K.py`  
**Tables:** `constants_DK.md`  
**Date:** 2026-07-12  

Each blocking item from `REFEREE_REPORT.md` is closed below with a pointer to the equation / lemma that implements the fix.

---

## 1. Luke / K-Bessel — no false lower bound; upper bound proved

| Referee demand | Fix |
|----------------|-----|
| Do **not** claim \(\sqrt{\pi/(2y)}\,e^{-y}\le\|K_{ir}(y)\|\) independent of \(r\) | Deleted. Paper states explicitly that no such lower bound is claimed (Lemma~\ref{lem:Kpoint} / eq.~\eqref{eq:Luke} discussion). |
| Prove or cite upper bound with \(C_K\) | **Full elementary proof** of the **sharp** bound \(C_K(r)=1\): Steps (A)–(C) from the integral representation \(K_{ir}(y)=\int_0^\infty e^{-y\cosh t}\cos(rt)\,dt\) — (A) \(\lvert\cos\rvert\le 1\Rightarrow\|K_{ir}\|\le K_0\); (B) monotonicity \(K_0\le K_{1/2}\); (C) closed form \(K_{1/2}=\sqrt{\pi/(2y)}\,e^{-y}\). Equations: \eqref{eq:Kint}, \eqref{eq:KA}, \eqref{eq:KB}, \eqref{eq:KC}, \eqref{eq:Kpoint}, \eqref{eq:CK}. |
| Classical \(e^{\pi r/2}\) | Retained as optional looser constant \(C_K^{\mathrm{cl}}\) \eqref{eq:CKclass}; theorems use sharp \(C_K=1\). |
| Optional poly refinement as UPPER only | \eqref{eq:Luke}: \(\|K_{ir}(y)\|\le\sqrt{\pi/(2y)}\,e^{-y}(1+1/(8y))\) for \(y\ge 1\). |
| Code docstrings | `lemma_K.py`: `C_K` default returns `1`; `classical=True` for \(e^{\pi r/2}\); comments state upper-only. |

---

## 2. Lemma K tail

| Referee demand | Fix |
|----------------|-----|
| Structure Hecke \(\lvert a_\beta\rvert^2\le C_H^2 n^{2(\theta+\varepsilon)}\), group by \(n\) | Unchanged algebra: Theorem~\ref{thm:lemmaK}, \eqref{eq:lemmaK}–\eqref{eq:lemmaKcrude}, proof via \eqref{eq:Ksq}. |
| Prefer \(r_2(n)\le 6d(n)\) as theorem default | Paper \eqref{eq:r2maj}; code default `r2_mode="div"`. Crude \(6n\) optional (`mode="crude"`). |
| Incomplete-gamma; state monotonicity threshold | Lemma~\ref{lem:integral}, \(N_{\mathrm{mono}}\) in \eqref{eq:Nmono}; code `monotonicity_threshold`. |
| Arb enclosure when flint present; M=100,200,400 | `lemma_K_tail` returns Arb; tests/bench green (see § benchmarks below). |

---

## 3. Trace inequality (Lemma T) — one consistent constant

| Referee demand | Fix |
|----------------|-----|
| Fix \(\lvert F\rvert/\lvert T\rvert\) vs \(\frac43\lvert F\rvert/\lvert T\rvert\) contradiction | **Single** constant throughout: \(C_{\mathrm{tr}}(T,F)=\frac43\frac{\lvert F\rvert}{\lvert T\rvert}\) \eqref{eq:Ctr}, derived from Young \eqref{eq:Young} in the proof of Lemma~\ref{lem:T}. |
| Same in `constants_DK.md` | §4: \(C_{\mathrm{trace}}^{\mathrm{cons}}=4/h_{\min}=16\). |
| Euclidean trace; hyp weights via \(y_{\min}\) | Stated in Lemma~\ref{lem:T} and Lemma~\ref{lem:metric}. |

---

## 4. Residual estimate — constant tracking

| Referee demand | Fix |
|----------------|-----|
| Expand Steps 1–2 with explicit factors | Step 1: \eqref{eq:step1a}–\eqref{eq:step1} (jump → trace reverse → \(e_{\mathrm{bdry}}\) → hyp norms). Step 2: \eqref{eq:step2a}–\eqref{eq:step2} (\(\Delta_{\mathrm{hyp}}\) pointwise bound, \(H^2\) scaling with \(h_{\min}\), \(A_{\mathrm{ell}}=1+\lambda\), \(A_{\mathrm{res}}\)). |
| Hyperbolic vs Euclidean on \(F_Y\subset\{y\ge y_{\min}\}\) | Lemma~\ref{lem:metric}: \eqref{eq:M-L2}, \eqref{eq:M-grad}, \(C_{\mathrm{met}}\) \eqref{eq:Cmet}, \(A_{\mathrm{met}}\) \eqref{eq:Amet}. |
| No magic undocumented factors | Factor map Remark~\ref{rem:factor-map}. |

---

## 5. Poincaré / Sobolev

| Referee demand | Fix |
|----------------|-----|
| PW/Bebendorf for convex **Euclidean** domains | Lemma~\ref{lem:PW}; metric comparison separated. |
| \(D_Y\le 10\) with containing box | Lemma~\ref{lem:Sobolev}, \eqref{eq:diam}. |
| \(C_{\mathrm{Poincaré}}=D_Y/\pi\), \(C_{\mathrm{Sob}}=1+C_{\mathrm{Poincaré}}\) | \eqref{eq:CPoincSob}, \eqref{eq:CPoinc}–\eqref{eq:CSob}. |

---

## 6. Theorem D(K) — \(C_1,C_2\), \(\eta_0\), continuous spectrum

| Referee demand | Fix |
|----------------|-----|
| \(C_1,C_2\) as SUM/PRODUCT of **named** contributions | \eqref{eq:Amet}–\eqref{eq:Acut}, \eqref{eq:C1}–\eqref{eq:C2}: \(C_1=2\,A_{\mathrm{bdry}}A_{\mathrm{res}}(1+A_{\mathrm{cusp}})\), \(C_2=A_{\mathrm{cut}}\). |
| Magic factor 4 removed | Prefactor is **2**, derived from \(\|U\|\ge 1/2\) in Step 4 \eqref{eq:step4b}–\eqref{eq:step4d}. |
| Revisit \((1+C_K)\) double-counting | **Removed.** \(C_K\) lives only inside \(A_{\mathrm{cusp}}\) \eqref{eq:Acusp}. |
| Non-circular \(\eta_0\) | \(\eta_0=\min\bigl(1/(4(1+A_{\mathrm{bdry}})^2),\,1/(4C_1^2)\bigr)\) \eqref{eq:eta0}, using a priori constants defined first. |
| Spectral pigeonhole proved | Step 4: \eqref{eq:step4c}–\eqref{eq:step4d} (second-moment / support intersection). |
| Continuous spectrum \(\lambda>1\) | Theorem stated for **cuspidal projection** \(U_{\mathrm{cusp}}\); Step 5 tracks BSV/Child soft-cutoff factor as \(\le A_{\mathrm{bdry}}\) into \(C_1\) \eqref{eq:step5}. |
| NOT proved list | Remark~\ref{rem:notproved}: no two-cusp, no \(N\mathfrak{p}\), no counting, no Weyl-as-count, no \(\mu(N)\). |

---

## 7. constants_DK.md + code

| Referee demand | Fix |
|----------------|-----|
| Regenerate tables at \((Y,Y_0,r)=(1.25,0.8,6.6)\) for both \(\lvert T\rvert\) | `constants_DK.md` fully regenerated. |
| Document \(y_{\min}\) | \(1/2\) for \(\mathbb{Z}[i]\), \(\sqrt{2/3}\) for \(\mathbb{Z}[\omega]\). |
| `C1_C2_constants(...)` | Implemented in `lemma_K.py`; returns Arb when flint present. |
| Run `python lemma_K.py --test --bench` | Passed (HAS_FLINT=True). |
| README success criteria | Updated. |

---

## Benchmark snapshot (post-revision, sharp \(C_K=1\))

```text
HAS_FLINT = True
C_K(6.6) sharp = 1
C_K(6.6) classical = [31801.0871994664 ± 8.69e-11]

tail M=100: [7.768951054136e-41 ± 5.95e-54]
tail M=200: [2.631582112047e-58 ± 1.58e-71]
tail M=400: [2.838307385819e-83 ± 6.11e-96]

C1(Z[i])  = [1324375.37123066 ± 1.69e-9]
C2(Z[i])  = [11.7280139288508 ± 4.66e-14]
C1(Z[ω])  = [356428.145921738 ± 6.00e-10]
C2(Z[ω])  = [3.34567592894081 ± 4.44e-15]
```

---

## Separation preserved

- No two-cusp theorem / block matrix in main theorems.  
- No \(N\mathfrak{p}\), no reference-cell gluing.  
- Then \(r_1\) remains an engineering target only.  
- Known / engineering / unknown separation kept in the abstract.  
- Hecke enters only via \(C_H,\theta\) in Lemma K / \(A_{\mathrm{cusp}}\).
