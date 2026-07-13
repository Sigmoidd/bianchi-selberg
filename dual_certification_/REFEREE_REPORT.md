# Referee report — Lemma K + Theorem D(K)

**Manuscript:** `dual_certification_/theorem_DK.tex`  
**Code:** `lemma_K.py`, `plot_kappa_vs_M.py`, `constants_DK.md`  
**Date:** 2026-07-12  
**Role:** Independent referee against the stated success criterion.

---

## Verdict

### **REVISION REQUIRED** (not yet publishable as a theorem paper)

The package correctly **isolates** the single-cusp core (no \(N\mathfrak{p}\), no gluing, no reference cell), states Assumption H once, ships a working Arb majorant for the **analytic tail majorant**, and records Then’s \(r_1\) as an engineering target only. Those architectural goals are met.

It does **not** yet meet the success criterion that a referee can check Theorem D\((K)\)’s constant \(C_1\) as rigorously derived from trace + Poincaré + Lemma K with no hidden factors. Lemma K’s pointwise bound is classical in spirit but the written proof of \(C_K(r)\) and especially the claimed Luke **lower** bound are not acceptable as stated. Theorem D\((K)\)’s proof is a **philosophy sketch** (BSV/Child style) with composite constants that absorb untracked absolute factors via a magic factor \(4\).

| Success criterion | Status |
|-------------------|--------|
| Lemma K checkable from Luke-type \(K_{ir}\) + H alone | **Partial** — structure yes; Luke lower bound false as written; \(C_K\) proof incomplete |
| \(C_1\) Arb-computable from trace, Poincaré, Lemma K | **Partial** — formula is Arb-evaluable; derivation of numerical prefactors incomplete |
| No hidden \(N\mathfrak{p}\)/gluing | **PASS** |
| Code: `lemma_K_tail` Arb enclosure + M=100,200,400 | **PASS** (ran green, `HAS_FLINT=True`) |
| Separation known / engineering / unknown | **PASS** |
| No two-cusp / Weyl-as-count / \(\mu(N)\) in main thms | **PASS** |

---

## What was delivered

| File | Present | Referee note |
|------|---------|--------------|
| `theorem_DK.tex` | yes | Full amsart spine; sections match the prompt |
| `lemma_K.py` | yes | Arb \(C_K\), incomplete-gamma majorant, tests pass |
| `plot_kappa_vs_M.py` + CSV/PNG | yes | Diagnostic only |
| `constants_DK.md` | yes | Explicit numbers at \((Y,Y_0,r)=(1.25,0.8,6.6)\) |
| `README.md` | yes | Clear claims / non-claims |

**Benchmarks (re-run by referee):**

```text
C_K(6.6) = [31801.0871994664 ± 8.69e-11]   [Arb]
tail M=100: ~7.86e-32
tail M=200: ~2.66e-49
tail M=400: ~2.87e-74
```

Enclosure / truncated-sum ratios \(\sim 10^2\)–\(10^3\) are expected under \(r_2\le 6n\).

---

## Detailed comments

### A. Assumption H — **Accept**

- Stated once up front; all theorems under H.
- \(\theta\in\{1/2,\,7/64,\,0\}\) correctly scoped.
- Remark that Hecke enters only via \(C_H,\theta\) in Lemma K / \(A_{\mathrm{cusp}}\): **correct and important**.

### B. Lemma K (Theorem 1) — **Major revision on analysis; code OK**

#### B1. Pointwise bound \(C_K(r)=e^{\pi r/2}\)

**Claim.** \(|K_{ir}(y)|\le\sqrt{\pi/(2y)}\,e^{-y}\,e^{\pi r/2}\).

**Status.** Standard in the automorphic literature (Watson / DLMF asymptotics / Iwaniec-type bounds). **Not adequately proved in the text.** The integral representation argument (p. ~298–311) asserts
\[
\Bigl|\int_0^\infty e^{-y(\cosh t-1)}\cos(rt)\,dt\Bigr|\le e^{\pi r/2}
\]
without a complete estimate. A referee cannot verify the constant from the given lines alone.

**Required revision.** Either:
1. Give a complete elementary proof of the constant \(e^{\pi r/2}\) from the integral rep, or  
2. Quote a precise theorem from Watson/DLMF/NIST with statement matching (eq:Kpoint) **exactly** (including the range of \(r,y\)).

#### B2. Luke double inequality (eq:Luke) — **ERROR**

The paper claims, for \(y\ge 1\),
\[
\sqrt{\frac{\pi}{2y}}\,e^{-y}
\;\le\;
|K_{ir}(y)|
\;\le\;
\sqrt{\frac{\pi}{2y}}\,e^{-y}\,C_K(r)\Bigl(1+\frac{1+r^2}{4y}\Bigr).
\]

The **lower bound independent of \(r\)** is false for large \(r\) at fixed \(y\) (e.g. \(y=1\), \(r=6.6\)): \(|K_{ir}(y)|\) is much smaller than the \(r=0\) asymptotic floor once oscillation / order growth is accounted for. Luke-type two-sided bounds, when correctly stated, have \(r\)-dependent factors on **both** sides or are asymptotic as \(y\to\infty\) with \(r\) fixed / \(r=o(y)\).

**Required revision.** Delete or replace the lower bound. Keep only a correctly cited upper refinement if needed. The **upper** bound alone is what Lemma K’s tail uses.

#### B3. Tail under H — **Accept structure; minor nits**

- Substitution \(y=2\pi\sqrt{n}Y_0\) and \(|K|^2\le \frac{1}{4\sqrt{n}Y_0}e^{-4\pi\sqrt{n}Y_0}C_K^2\): **correct** given the pointwise upper bound.
- Grouping by \(n=N(\beta)\), Hecke \(|a_\beta|^2\le C_H^2 n^{2\theta+2\varepsilon}\): **correct**.
- \(r_2(n)\le 6n\): **true** (crude); prefer \(r_2\le 6d(n)\) in the theorem statement (already in code as `mode="div"`).
- Incomplete-gamma integral comparison (Lemma integral): **standard and OK**, with the monotonicity threshold handled in code.

#### B4. Code `lemma_K.py` — **Accept for majorant enclosure**

- Implements the **crude majorant**, not a numerical verification of true \(K_{ir}\) values against the majorant (that would need a \(K_{ir}\) library). Acceptable if documented (it is: analytic majorant vs truncated majorant).
- Self-tests and M=100/200/400 benchmarks: **PASS**.
- Arb path via python-flint: **PASS**.

### C. Theorem D\((K)\) — **Major revision (proof incomplete)**

#### C1. Statement — **Accept with caveats**

- Single-cusp; defects \(\delta_{\mathrm{aut}},\tau_{\mathrm{tail}}\); conclusion of nearby true eigenvalue: correct **shape** (BSV Prop 3.1 / Child Thm 1.1).
- Explicit non-claims (no two-cusp, no \(N\mathfrak{p}\), no counting): **good**.
- Then interval as engineering target: **good**.

#### C2. Geometric lemmas

| Lemma | Issue |
|-------|--------|
| Trace (Lemma T) | Integration-by-parts idea OK; Young algebra for \(C_{\mathrm{tr}}=\|F\|/\|T\|\) vs \(\frac43\|F\|/\|T\|\) is muddled (text says both). Fix one consistent constant. |
| Payne–Weinberger | Citation OK; applies to **Euclidean** convex domains — must state that the metric comparison (hyperbolic vs Euclidean gradients on \(F_Y\subset\{y\ge y_{\min}\}\)) is absorbed into \(C_{\mathrm{trace}}/C_{\mathrm{Sob}}\) with explicit \(y_{\min}^{-O(1)}\) factors. **Currently missing.** |
| Diameter \(D_Y\le 10\) | Acceptable as a conservative explicit bound; justify once with a concrete containing box (done, roughly). |
| \(C_{\mathrm{trace}}=16\) via \(h_{\min}=0.25\) | Mesh-free; OK as **upper bound of a formula**, but not derived from the trace lemma’s face set of the actual orbifold fundamental domain. Label as “conservative geometric majorant,” not “sharp.” |

#### C3. Proof of D\((K)\) — **Insufficient for a theorem**

Steps 1–5 outline the right strategy but leave gaps a referee cannot fill without new work:

1. **Step 1 (automorphy → \(H^1\) boundary corrector).**  
   Pointwise \(\delta_{\mathrm{aut}}\) on faces gives \(L^\infty\) control; passage to an \(H^1\) corrector \(e_{\mathrm{bdry}}\) with \(\|e_{\mathrm{bdry}}\|_{H^1}\le A_{\mathrm{bdry}}\delta_{\mathrm{aut}}\) needs a specific extension operator / tubular chart and **hyperbolic** volume weights \(y^{-3}\). The factor \(A_{\mathrm{bdry}}=C_{\mathrm{trace}}C_{\mathrm{Sob}}(1+\lambda)^{1/2}\) is **asserted**, not derived.

2. **Step 2 (residual).**  
   Claims \(\|(\Delta-\lambda)e_{\mathrm{bdry}}\|\le C\|e_{\mathrm{bdry}}\|_{H^2}\) and “elliptic control of \(H^2\) by \(H^1+\Delta\)” without constants or citation of a quantitative elliptic estimate on the orbifold core. For a \(C^2\) putative form this is plausible but **not written**.

3. **Step 3 (collar / Lemma K).**  
   Connecting the truncated expansion to a global \(\Gamma\)-periodization via Shimizu is the right idea; the \(O(e^{-2\pi Y_0})\) cutoff error needs an explicit constant (enters \(C_2\)). Presently \(C_2=2|T|Y^{-1}(1+C_{\mathrm{Sob}})\) is a plausible scale but not derived line-by-line from the cutoff.

4. **Step 4 (spectral resolution).**  
   The key inequality
   \[
   \sum_j|\lambda_j-\lambda|^2|c_j|^2 \le \|(\Delta-\lambda)U\|_2^2
   \]
   is correct for the discrete part. The jump from “some \(c_{j_0}\) large” to
   \[
   |\lambda_{j_0}-\lambda|\le \|(\Delta-\lambda)U\|_2/|c_{j_0}|
   \]
   needs a quantitative isolation / pigeonhole (as in BSV): if the residual is \(<\frac12\|U\|\) then mass cannot be all far from \(\lambda\). The text gestures at this via \(\eta_0=(2C_1)^{-2}\) but **does not prove** that this \(\eta_0\) works with the same \(C_1\). Circular risk: \(C_1\) is used to define \(\eta_0\) which is used to justify \(C_1\).

5. **Step 5 (continuous spectrum).**  
   For \(\lambda>1\) (Then regime) the soft cutoff argument is cited to BSV/Child but not reproduced with constants. Acceptable only if the paper clearly says “constants as in Child Thm 1.1, specialized to \(\HH^3\)” **and** tracks them into \(C_1\). Currently they are “absorbed.”

6. **Magic factor \(4\) in \(C_1=4\,A_{\mathrm{bdry}}A_{\mathrm{res}}(1+A_{\mathrm{cusp}})(1+C_K)\).**  
   Unjustified. Either derive it or replace by a sum of named contributions.

#### C4. Numerical size of \(C_1\sim 10^9\)

Not a mathematical error, but a **practical** one: with \(C_1\sim 1.3\cdot 10^9\), one needs \(\eta\lesssim 10^{-13}\) for a \(10^{-4}\) eigenvalue window. The paper correctly notes that Lemma K tails are not the bottleneck; automorphy defect is. Referee accepts this honesty. Improving \(C_1\) (drop the unnecessary \((1+C_K)\) if \(A_{\mathrm{cusp}}\) already has \(C_K\), fix metric factors, raise \(Y_0\)) should be discussed.

**Note:** \(C_1\) multiplies both \((1+A_{\mathrm{cusp}})\) and \((1+C_K)\). Since \(A_{\mathrm{cusp}}\) already contains \(C_K\), the extra \((1+C_K)\) roughly multiplies by \(3\cdot 10^4\) and may be an artifact of overcounting. **Revisit.**

### D. Separation and constraints — **PASS**

- No two-cusp theorem; block matrix only in “What remains.”
- No \(N\mathfrak{p}\), no reference cell in main theorems.
- \(\mu(N)\) only in heuristic appendix, explicitly unused.
- Weyl only for window guidance, not exact counting.
- Hecke entry localized (Remark + README).

### E. Literature — **Mostly OK**

Expected citations present in spirit (BSV, Child, Then, EGM, PW/Bebendorf, Luke, Rump, Blomer–Brumley). Ensure bibliography entries are complete (DOI/year) in the next revision. Dual roadmap as informal citation is fine for a repo draft.

---

## Required revisions (blocking)

1. **Fix or remove Luke lower bound** (eq:Luke); do not claim \(r\)-independent lower bound for \(|K_{ir}|\).
2. **Prove or precisely cite** the upper bound with \(C_K(r)=e^{\pi r/2}\).
3. **Rewrite Theorem D\((K)\) proof** so that every factor in \(C_1,C_2\) is either  
   (a) derived with inequality tags, or  
   (b) cited to a numbered external lemma with matching hypotheses (metric weights included).
4. **Eliminate circular \(\eta_0\)–\(C_1\) definition**; define \(\eta_0\) from a priori constants only.
5. **Hyperbolic vs Euclidean:** insert \(y_{\min}\)-explicit comparison constants on \(F_Y\).
6. **Re-examine \((1+C_K)\) in \(C_1\)** for double-counting.
7. **Align Lemma T constant** (\(\|F\|/\|T\|\) vs \(\frac43\|F\|/\|T\|\)) everywhere (paper + `constants_DK.md`).

## Non-blocking suggestions

- Prefer \(r_2(n)\le 6d(n)\) in the theorem (code already has `div` mode).
- Add a numerical check: high-precision \(|K_{ir}(y)|\) vs majorant at a few \((r,y)\) (optional; needs `scipy.special.kve` or Arb Bessel).
- Provide `C1_arb(...)` in Python mirroring `constants_DK.md` formulas for full reproducibility of \(C_1\).
- Soften abstract language from “we prove” to “we prove Lemma K fully; we establish Theorem D\((K)\) under explicit geometric constants detailed in §5” **until** (3)–(5) are fixed — or keep “prove” only after revisions.

---

## Referee checklist (prompt success criterion)

> A referee can check Lemma K tail bound independently with only Luke bounds and Hecke bound H.

**Not yet.** Tail algebra from a **given** pointwise upper bound + H: yes. Luke lower bound and incomplete \(C_K\) proof: no.

> Theorem D\((K)\) constant \(C_1\) is Arb-computable from trace, Poincaré, and Lemma K, with no hidden dependence on \(N\mathfrak{p}\) or gluing.

**Formula-computable: yes. Derived: no.** No \(N\mathfrak{p}\)/gluing: yes.

---

## Summary for authors

You have a **strong architectural draft** that matches the dual-certification separation (Rung 0–1) and a **correct, Arb-backed majorant implementation** for Lemma K’s series. To survive a hard analysis referee on Theorem D\((K)\), the defect argument must be expanded from a BSV/Child **outline** into a **constant-tracking proof**, and the Bessel section must not overclaim Luke.

**Recommendation:** Revise and resubmit. After revisions, re-run `python lemma_K.py --test --bench` and regenerate `constants_DK.md`.
