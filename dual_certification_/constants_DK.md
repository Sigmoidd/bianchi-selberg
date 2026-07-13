# Explicit constants for Lemma K and Theorem D(K)

Parameters used throughout this file (unless noted):

| Symbol | Value | Meaning |
|--------|-------|---------|
| \(Y\) | 1.25 | core truncation height |
| \(Y_0\) | 0.8 | collar / Fourier evaluation height |
| \(r\) | 6.6 | spectral parameter (near Then \(r_1\approx 6.62212\)) |
| \(\lambda\) | \(r^2+1 = 44.56\) | Laplace eigenvalue parameter |
| \(\theta\) | \(1/2\) | Hecke exponent (Rankin–Selberg; unconditional) |
| \(C_H\) | 1.0 | Hecke constant (unit convention; rescale with form) |
| \(\varepsilon\) | 0.0 | Hecke \(\varepsilon\) |
| \(D_Y\) | 10 | conservative Euclidean diameter of a body containing \(F_Y\) |
| \(h_{\min}\) | 0.25 | conservative minimal altitude for mesh-free trace majorant |
| \(M\) | 400 | Fourier truncation (for \(A_{\mathrm{cusp}}\)) |

**Backend.** Values tagged **[Arb]** are produced with `python-flint` ball arithmetic (certifying enclosure of the *formula*). Values tagged **[float]** are non-certifying diagnostics.

Regenerate:

```text
python lemma_K.py --test --bench --constants
```

API: `from lemma_K import C1_C2_constants, lemma_K_tail, C_K`.

---

## 1. Geometric data \(|T|\) and \(y_{\min}\)

| Field | \(\mathcal{O}_K\) | \(\|T\|\) (exact) | \(\|T\|\) numeric | \(y_{\min}\) (conservative) |
|-------|-------------------|-------------------|-------------------|------------------------------|
| \(\mathbb{Q}(i)\) | \(\mathbb{Z}[i]\) | \(1/2\) | 0.5 | \(1/2 = 0.5\) |
| \(\mathbb{Q}(\sqrt{-3})\) | \(\mathbb{Z}[\omega]\) | \(\sqrt{3}/6\) | 0.2886751345948129… | \(\sqrt{2/3}\approx 0.816496580927726\) |

**Choice of \(y_{\min}\).** On the Ford/Humbert floor one has \(y\ge\sqrt{1-\max|z|^2}\) over the planar section. The values above are **conservative lower bounds** (strictly below the classical Humbert floor \(1/\sqrt{2}\) for \(\mathbb{Z}[i]\), and a standard conservative floor for the EGM section of \(\mathbb{Z}[\omega]\)). They may be sharpened without changing the structure of the constants.

Hyperbolic cusp tail volume: \(\mathrm{vol}(C_Y)=\|T\|/(2Y^2)\).

At \(Y=1.25\):

| Field | \(\mathrm{vol}(C_Y)\) |
|-------|------------------------|
| \(\mathbb{Z}[i]\) | \(0.5/(2\cdot 1.5625)=0.16\) |
| \(\mathbb{Z}[\omega]\) | \(\sqrt{3}/6/(2\cdot 1.5625)\approx 0.09237\) |

---

## 2. \(C_K(r)\) — K-Bessel prefactor (UPPER bound only)

**Sharp (theorem default), elementary proof** (paper Lemma Kpoint, Steps A–C):

\[
|K_{ir}(y)|\le K_0(y)\le K_{1/2}(y)=\sqrt{\frac{\pi}{2y}}\,e^{-y},
\qquad
C_K(r):=1.
\]

**Classical looser constant** (optional, literature comparison):

\[
C_K^{\mathrm{cl}}(r)=\exp\Bigl(\frac{\pi r}{2}\Bigr)\ge 1.
\]

**Optional poly envelope for \(y\ge 1\) only** (upper refinement):

\[
C_K^{\mathrm{poly,sharp}}=1+\tfrac18,
\qquad
C_K^{\mathrm{poly,cl}}(r)=C_K^{\mathrm{cl}}(r)\Bigl(1+\frac{1+r^2}{4}\Bigr).
\]

No \(r\)-independent **lower** bound is claimed.

At \(r=6.6\):

| Quantity | Value | Tag |
|----------|-------|-----|
| \(C_K(6.6)\) sharp | `1` | **[Arb]** |
| \(C_K^{\mathrm{cl}}(6.6)\) | `[31801.0871994664 ± 8.69e-11]` | **[Arb]** |
| \(C_K^{\mathrm{poly,sharp}}\) | `1.125` | **[Arb]** |
| \(C_K^{\mathrm{cl}}(6.62212)\) (Then target) | \(\approx 3.216\cdot 10^{4}\) | **[float]** |

---

## 3. Lemma K tail majorant \(S_{M,\infty}\)

Formula (theorem default \(r_2(n)\le 6d(n)\le 6n\), sharp \(C_K=1\)):

\[
S_{M,\infty}
\le
\frac{3\,C_H^2\,C_K(r)^2}{2Y_0}
\sum_{n>M}
n^{2\theta+2\varepsilon+1/2}\,
e^{-4\pi\sqrt{n}\,Y_0},
\]

enclosed via upper incomplete gamma after the monotonicity threshold
\(N_{\mathrm{mono}}=\lceil(2\alpha/c)^2\rceil+1\) (\(\alpha>0\)).

At \((Y_0,r,\theta,C_H,\varepsilon)=(0.8,6.6,1/2,1,0)\), **sharp \(C_K=1\)**:

| \(M\) | Tail enclosure | Truncated majorant sum \(n\le 20\mathrm{k}\) (exact \(r_2\) for \(\mathbb{Z}[i]\)) | Tag |
|------:|----------------|-------------------------------------------------------------------:|-----|
| 100 | `[7.768951054136e-41 ± 5.95e-54]` | \(\approx 4.46\cdot 10^{-43}\) | **[Arb]** |
| 200 | `[2.631582112047e-58 ± 1.58e-71]` | \(\approx 5.86\cdot 10^{-61}\) | **[Arb]** |
| 400 | `[2.838307385819e-83 ± 6.11e-96]` | \(\approx 4.17\cdot 10^{-86}\) | **[Arb]** |

Enclosure / truncation ratios \(\sim 10^2\)–\(10^3\) are the price of \(r_2\le 6n\); both sides are negligible for any realistic defect \(\eta\).

*(With classical \(C_K^{\mathrm{cl}}\), multiply enclosures by \(C_K^{\mathrm{cl}}(6.6)^2\approx 1.01\cdot 10^{9}\).)*

---

## 4. Trace constant \(C_{\mathrm{trace}}\) — **one consistent formula**

From Lemma T (paper §5): for a tetrahedron \(T\) and face \(F\), after Young,

\[
C_{\mathrm{tr}}(T,F)=\frac43\frac{|F|}{|T|},
\qquad
C_{\mathrm{trace}}(K,Y)=\max_{T,F\subset\partial F_Y}C_{\mathrm{tr}}(T,F).
\]

This is a **Euclidean** trace inequality; hyperbolic weights are converted by the metric-comparison constants of §6.

**Mesh-free conservative bound** (no FEM dependence).  
For a right tetrahedron of altitude \(h\), \(|F|/|T|=3/h\), so

\[
C_{\mathrm{trace}}^{\mathrm{cons}}
=\frac43\cdot\frac{3}{h_{\min}}
=\frac{4}{h_{\min}}.
\]

Taking \(h_{\min}=0.25\):

\[
C_{\mathrm{trace}}^{\mathrm{cons}}=16.
\]

| Field | \(C_{\mathrm{trace}}\) used | Tag |
|-------|----------------------------:|-----|
| both | 16 | **[float]** formula-exact with \(h_{\min}=0.25\); conservative geometric majorant, not sharp |

---

## 5. Poincaré / Sobolev constants

Payne–Weinberger / Bebendorf: for **convex Euclidean** domains,

\[
\|w-\bar w\|_{L^2(\Omega)}\le\frac{\mathrm{diam}(\Omega)}{\pi}\,\|\nabla_{\mathrm{euc}} w\|_{L^2}.
\]

| Bound | Value | Justification |
|-------|------:|---------------|
| sharp-ish box diam | \(<3\) | Picard box sides \((1,0.5,Y-y_{\min})\); EGM prism comparable |
| **conservative \(D_Y\)** | **10** | safety margin for charts / partition of unity |

\[
C_{\mathrm{Poincaré}}=\frac{D_Y}{\pi}=\frac{10}{\pi},
\qquad
C_{\mathrm{Sob}}=1+C_{\mathrm{Poincaré}}.
\]

| Constant | Value | Tag |
|----------|------:|-----|
| \(D_Y\) | 10 | conservative explicit |
| \(C_{\mathrm{Poincaré}}\) | `[3.18309886183791 ± 2.98e-15]` | **[Arb]** from \(10/\pi\) |
| \(C_{\mathrm{Sob}}\) | `[4.18309886183791 ± 2.98e-15]` | **[Arb]** |

---

## 6. Metric comparison \(C_{\mathrm{met}}\)

On \(F_Y\subset\{y_{\min}\le y\le Y\}\) (paper Lemma metric):

\[
\begin{aligned}
C_{\mathrm{met},0}&=y_{\min}^{-3/2},&
C_{\mathrm{met},0}^{-}&=Y^{3/2},\\
C_{\mathrm{met},1}&=y_{\min}^{-1/2},&
C_{\mathrm{met},1}^{-}&=Y^{1/2},\\
C_{\mathrm{met}}&=\max\bigl(C_{\mathrm{met},0},C_{\mathrm{met},1},C_{\mathrm{met},0}^{-},C_{\mathrm{met},1}^{-}\bigr).
\end{aligned}
\]

| Field | \(y_{\min}\) | \(C_{\mathrm{met}}\) | Tag |
|-------|-------------:|---------------------:|-----|
| \(\mathbb{Z}[i]\) | 0.5 | `[2.82842712474619 ± 2.92e-16]` \(=2^{3/2}\) | **[Arb]** |
| \(\mathbb{Z}[\omega]\) | \(\sqrt{2/3}\) | `[1.39754248593737 ± 1.38e-15]` \(=Y^{3/2}\) | **[Arb]** |

---

## 7. Composite constants \(C_1\), \(C_2\) for Theorem D(K)

Paper formulas (named contributions; **no** magic factor 4; **no** extra \((1+C_K)\)):

\begin{align*}
A_{\mathrm{met}}
&=C_{\mathrm{met}},\\
A_{\mathrm{bdry}}
&=C_{\mathrm{trace}}\,C_{\mathrm{Sob}}\,A_{\mathrm{met}}\,(1+\sqrt{\lambda}),\\
A_{\mathrm{ell}}
&=1+\lambda,\\
A_{\mathrm{res}}
&=A_{\mathrm{ell}}\bigl(1+C_{\mathrm{Poincaré}}\,A_{\mathrm{met}}\bigr),\\
A_{\mathrm{cusp}}
&=\frac{C_H\,C_K(r)}{2\sqrt{Y_0}}
\Bigl(\sum_{n>M}6n\cdot n^{2\theta+2\varepsilon-1/2}e^{-4\pi\sqrt{n}Y_0}\Bigr)^{1/2},\\
A_{\mathrm{cut}}
&=2\,|T|\,Y^{-1}\,(1+C_{\mathrm{Sob}})\,A_{\mathrm{met}},\\
C_1
&=2\,A_{\mathrm{bdry}}\,A_{\mathrm{res}}\,(1+A_{\mathrm{cusp}}),\\
C_2
&=A_{\mathrm{cut}},\\
\eta_0
&=\min\Bigl(\tfrac{1}{4(1+A_{\mathrm{bdry}})^2},\,\tfrac{1}{4C_1^2}\Bigr).
\end{align*}

The prefactor \(2\) in \(C_1\) is the \(\|U\|\ge 1/2\) normalization from the spectral pigeonhole (paper Step 4).  
\(C_K\) enters **only** inside \(A_{\mathrm{cusp}}\) (and is \(1\) in the sharp case).  
Hecke \((C_H,\theta,\varepsilon)\) enters **only** inside \(A_{\mathrm{cusp}}\) / Lemma K.

### Numerical evaluation at \((Y,Y_0,r,M)=(1.25,0.8,6.6,400)\), sharp \(C_K=1\)

\(\lambda=44.56\), \(\sqrt{\lambda}\approx 6.675326\).

| Quantity | \(\mathbb{Z}[i]\) | \(\mathbb{Z}[\omega]\) | Tag |
|----------|-----------------:|----------------------:|-----|
| \(A_{\mathrm{met}}\) | \(\approx 2.828427\) | \(\approx 1.397542\) | **[Arb]** |
| \(A_{\mathrm{bdry}}\) | `[1452.98132393805 ± 4.02e-12]` | `[717.926622082993 ± 3.91e-13]` | **[Arb]** |
| \(A_{\mathrm{ell}}\) | 45.56 | 45.56 | **[Arb]** |
| \(A_{\mathrm{res}}\) | `[455.744113641178 ± 4.09e-13]` | `[248.234384238041 ± 2.02e-13]` | **[Arb]** |
| \(A_{\mathrm{cusp}}\) (\(M=400\)) | `[5.327576734144e-42 ± 3.45e-55]` | same | **[Arb]** |
| **\(C_1\)** | **`[1324375.37123066 ± 1.69e-9]`** | **`[356428.145921738 ± 6.00e-10]`** | **[Arb]** |
| **\(C_2\)** | **`[11.7280139288508 ± 4.66e-14]`** | **`[3.34567592894081 ± 4.44e-15]`** | **[Arb]** |
| \(\eta_0\) | \(\approx 1.43\cdot 10^{-13}\) | \(\approx 1.97\cdot 10^{-12}\) | **[Arb]** |

### Factor map (referee checklist)

| Factor in \(C_1\) / \(C_2\) | Named source |
|----------------------------|--------------|
| \(C_{\mathrm{trace}}=16\) | Lemma T, \(C_{\mathrm{tr}}=\frac43\|F\|/\|T\|\), cons.\ \(h_{\min}\) |
| \(C_{\mathrm{Sob}}=1+D_Y/\pi\) | Payne–Weinberger + mean control |
| \(A_{\mathrm{met}}=C_{\mathrm{met}}\) | hyperbolic vs Euclidean on \(\{y\ge y_{\min}\}\) |
| \(A_{\mathrm{bdry}}\) | Step 1 (automorphy → \(H^1\) corrector) |
| \(A_{\mathrm{ell}},A_{\mathrm{res}}\) | Step 2 (residual / elliptic scaling) |
| \(A_{\mathrm{cusp}}\) | Lemma K + collar (Hecke only here) |
| \(A_{\mathrm{cut}}=C_2\) | Step 3 (collar cutoff) |
| factor \(2\) in \(C_1\) | Step 4 (\(\|U\|\ge 1/2\)) |
| \(\eta_0\) | defined from a priori \(A_{\mathrm{bdry}},C_1\) only |

### Practical consequence

- \(C_2 e^{-2\pi Y_0}\): with \(e^{-2\pi\cdot 0.8}\approx 6.55\cdot 10^{-3}\), one gets \(\approx 0.077\) for \(\mathbb{Z}[i]\) at these parameters. For a tight eigenvalue window, increase \(Y_0\) or sharpen the cutoff error in \(C_2\).
- \(C_1\eta\): to obtain \(|\lambda_{\mathrm{true}}-\lambda|\le 10^{-4}\) one needs roughly \(\eta\lesssim 10^{-10}\) (\(\mathbb{Z}[i]\)). The Lemma K tail at \(M=400\) is already \(\ll 10^{-80}\) and is not the bottleneck; automorphy defect quality is.
- Improvement vs pre-revision: sharp \(C_K=1\) and removal of double-counted \((1+C_K)\) drop \(C_1\) by \(\sim 10^3\) relative to the old \(\sim 10^9\) composite that used \(C_K^{\mathrm{cl}}\) twice.

---

## 8. Summary table (copy-paste)

| Constant | Formula | Value at \((Y,Y_0,r)=(1.25,0.8,6.6)\) |
|----------|---------|----------------------------------------|
| \(\|T\|_{\mathbb{Z}[i]}\) | \(1/2\) | 0.5 |
| \(\|T\|_{\mathbb{Z}[\omega]}\) | \(\sqrt{3}/6\) | 0.288675… |
| \(y_{\min}(\mathbb{Z}[i])\) | conservative | 0.5 |
| \(y_{\min}(\mathbb{Z}[\omega])\) | \(\sqrt{2/3}\) | 0.81650… |
| \(C_K(r)\) sharp | \(1\) | 1 **[Arb]** |
| \(C_K^{\mathrm{cl}}(r)\) | \(e^{\pi r/2}\) | \(\approx 31801.087\) **[Arb]** |
| \(C_{\mathrm{trace}}\) | \(\frac43\max\|F\|/\|T\|\) cons. | 16 |
| \(C_{\mathrm{Poincaré}}\) | \(D_Y/\pi\), \(D_Y=10\) | \(10/\pi\approx 3.1831\) **[Arb]** |
| \(C_{\mathrm{Sob}}\) | \(1+C_{\mathrm{Poincaré}}\) | \(\approx 4.1831\) **[Arb]** |
| \(C_{\mathrm{met}}(\mathbb{Z}[i])\) | \(\max y_{\min}^{\pm O(1)},Y^{O(1)}\) | \(\approx 2.8284\) **[Arb]** |
| \(C_1\) (\(\mathbb{Z}[i]\)) | \(2 A_{\mathrm{bdry}}A_{\mathrm{res}}(1+A_{\mathrm{cusp}})\) | \(\approx 1.324\cdot 10^{6}\) **[Arb]** |
| \(C_1\) (\(\mathbb{Z}[\omega]\)) | same | \(\approx 3.564\cdot 10^{5}\) **[Arb]** |
| \(C_2\) (\(\mathbb{Z}[i]\)) | \(A_{\mathrm{cut}}\) | \(\approx 11.728\) **[Arb]** |
| \(C_2\) (\(\mathbb{Z}[\omega]\)) | same | \(\approx 3.346\) **[Arb]** |
| Lemma K tail \(M=400\) | incomplete-gamma majorant, \(C_K=1\) | \(\sim 2.8\cdot 10^{-83}\) **[Arb]** |

---

## 9. What enters where

| Input | Enters |
|-------|--------|
| Assumption H \((\theta,C_H,\varepsilon)\) | \(A_{\mathrm{cusp}}\), Lemma K tail only |
| Assumption A (analyticity) | Theorem D(K) statement / unique continuation |
| Assumption S (spectral decomp.) | Theorem D(K) spectral pigeonhole |
| Elementary \(K_{ir}\) upper bound | \(C_K(r)=1\) (sharp) |
| Trace inequality (Lemma T) | \(C_{\mathrm{trace}}\to A_{\mathrm{bdry}}\to C_1\) |
| Payne–Weinberger | \(C_{\mathrm{Poincaré}}\to A_{\mathrm{res}},C_{\mathrm{Sob}}\) |
| Diameter bound \(D_Y\le 10\) | \(C_{\mathrm{Poincaré}}\) |
| Metric comparison \(y_{\min},Y\) | \(A_{\mathrm{met}}\to A_{\mathrm{bdry}},A_{\mathrm{res}},C_2\) |
| \(\|T\|\) | \(C_2\) only |
| **Not used** | \(N\mathfrak{p}\), gluing, reference cell, Weyl count, \(\mu(N)\) |

---

## 10. Evaluation at Then engineering target \(r=6.622\) (not certified \(\lambda\))

Regenerate: `python -c "from lemma_K import C1_C2_constants, ball_str; d=C1_C2_constants(r=6.622); print(ball_str(d['C1']), ball_str(d['C2']))"`

At \((Y,Y_0,r,M)=(1.25,0.8,6.622,400)\), sharp \(C_K=1\), \(\theta=1/2\), \(C_H=1\):

| Quantity | \(\mathbb{Z}[i]\) | Tag |
|----------|------------------:|-----|
| \(\lambda=r^2+1\) | \(1+6.622^2=44.850884\) | exact float input |
| **\(C_1\)** | **`[1336608.38547962 ± 1.43e-9]`** | **[Arb]** |
| **\(C_2\)** | **`[11.7280139288508 ± 4.66e-14]`** | **[Arb]** (independent of \(r\)) |
| Lemma K tail \(M=400\) | `[2.838307385819e-83 ± 6.11e-96]` | **[Arb]** (same exponential; \(C_K=1\)) |

Then \(r_1\approx 6.62212\) remains an **engineering target** only. The values above enclose the *formula* for \(C_1,C_2\); they do not certify an eigenvalue.

---

## 11. Sensitivity of \(C_1,C_2\) to \(Y_0\pm 0.1\), \(r\pm 0.1\)

Regenerate: `python defect_bound_arb.py --sensitivity` or `--md-sensitivity`.

**Structural note.**  
- \(C_2=A_{\mathrm{cut}}\) depends on \((K,Y)\) only (via \(\|T\|\), \(Y\), \(C_{\mathrm{Sob}}\), \(A_{\mathrm{met}}\)) — **independent of \(Y_0\) and \(r\)**.  
- \(C_1\) depends on \(r\) through \(A_{\mathrm{bdry}}\propto(1+\sqrt{\lambda})\) and \(A_{\mathrm{res}}\propto(1+\lambda)\); dependence on \(Y_0\) is only through the negligible \(A_{\mathrm{cusp}}\) at \(M=400\).  
- The **eigenvalue-error bound** \(C_1\eta+C_2 e^{-2\pi Y_0}\) *does* depend strongly on \(Y_0\) via the exponential collar term.

### \(\mathbb{Z}[i]\) at \(Y=1.25\), \(M=400\), \(\theta=1/2\), \(C_H=1\) **[Arb]**

| Config | \(Y_0\) | \(r\) | \(C_1\) | \(C_2\) | \(C_2 e^{-2\pi Y_0}\) |
|--------|--------:|------:|--------:|--------:|----------------------:|
| baseline | 0.80 | 6.60 | \(1.324375\cdot 10^{6}\) | 11.7280 | 0.076952 |
| \(Y_0+0.1\) | 0.90 | 6.60 | \(1.324375\cdot 10^{6}\) | 11.7280 | 0.041053 |
| \(Y_0-0.1\) | 0.70 | 6.60 | \(1.324375\cdot 10^{6}\) | 11.7280 | 0.144244 |
| \(r+0.1\) | 0.80 | 6.70 | \(1.380598\cdot 10^{6}\) | 11.7280 | 0.076952 |
| \(r-0.1\) | 0.80 | 6.50 | \(1.269728\cdot 10^{6}\) | 11.7280 | 0.076952 |
| \(Y_0+0.1,\,r+0.1\) | 0.90 | 6.70 | \(1.380598\cdot 10^{6}\) | 11.7280 | 0.041053 |
| \(Y_0-0.1,\,r-0.1\) | 0.70 | 6.50 | \(1.269728\cdot 10^{6}\) | 11.7280 | 0.144244 |

Relative change in \(C_1\) for \(\Delta r=\pm 0.1\): about \(\pm 4.2\%\).  
Relative change in collar term for \(\Delta Y_0=\pm 0.1\): factor \(\approx e^{\pm 0.2\pi}\approx 1.87^{\pm 1}\).

### \(\mathbb{Z}[\omega]\) (same parameters) **[Arb]**

| Config | \(Y_0\) | \(r\) | \(C_1\) | \(C_2\) | \(C_2 e^{-2\pi Y_0}\) |
|--------|--------:|------:|--------:|--------:|----------------------:|
| baseline | 0.80 | 6.60 | \(3.56428\cdot 10^{5}\) | 3.34568 | 0.021952 |
| \(Y_0+0.1\) | 0.90 | 6.60 | \(3.56428\cdot 10^{5}\) | 3.34568 | 0.011711 |
| \(Y_0-0.1\) | 0.70 | 6.60 | \(3.56428\cdot 10^{5}\) | 3.34568 | 0.041149 |
| \(r+0.1\) | 0.80 | 6.70 | \(3.71559\cdot 10^{5}\) | 3.34568 | 0.021952 |
| \(r-0.1\) | 0.80 | 6.50 | \(3.41721\cdot 10^{5}\) | 3.34568 | 0.021952 |
| \(Y_0+0.1,\,r+0.1\) | 0.90 | 6.70 | \(3.71559\cdot 10^{5}\) | 3.34568 | 0.011711 |
| \(Y_0-0.1,\,r-0.1\) | 0.70 | 6.50 | \(3.41721\cdot 10^{5}\) | 3.34568 | 0.041149 |

---

## 12. Hejhal plug-in API (Rung 1 deliverable)

```text
from defect_bound_arb import defect_to_lambda_error

out = defect_to_lambda_error(delta, tau, Y=1.25, Y0=0.8, r=6.6,
                             field='i', M=400, theta=0.5, C_H=1.0)
# out['lambda_error_bound']  = Arb enclosure of C1*(δ+τ) + C2*exp(-2π Y0)
# out['L2_error_bound']      = Arb enclosure of C1*sqrt(δ+τ)  (needs η≤η0)
# out['assumptions']         = [H, A, S]
```

Demo: `python defect_bound_arb.py --demo`  
Sensitivity: `python defect_bound_arb.py --sensitivity`

**Language.** The API returns enclosures of the *Theorem D(K) defect bound under Assumptions H, A, S*. It does **not** by itself certify an eigenvalue: a later Hejhal run must supply certified \((\delta,\tau)\) with \(\eta\le\eta_0\).
