# Audit: load-bearing claims for ℤ[ω] cuspidal elliptic (CE)

**Status:** *not* accepted as paper-ready without independent writeup.  
This note answers a collaborator audit of four claims. Tone: separate
**verified software** from **mathematics that must appear as propositions**.

Primary sources in-repo:

| Claim | Code | Math note |
|-------|------|-----------|
| CE formula / kernel | `cuspidal_ce.py` header, `bianchi_omega*.py` | Friedman arXiv:math/0612807 Lemmas 4.3.2, 4.4.4 (OCR: `bianchiselberg-refs/friedman_thesis.txt`) |
| Six classes + \|C\|=6 + constants | `cuspidal_ce.exact_eisenstein_CE_check`, `run_CE("omega")` | `old_RIGOR_GAPS.md` §2 |
| NCE classification | `elliptic_inventory.py` | EGM book Ch.4 §4.3 (cited in RIGOR_GAPS §1; **not** substituted by 1982 OCR) |

---

## Q1. Where is the Eisenstein kernel \(\dfrac{\sinh x}{\cosh x+\tfrac12}\) derived?

**Answer: Friedman’s general CE integral, specialized to order-3.**

Friedman Lemma 4.3.2 (trace formula CE contribution) has the shape

\[
\int_0^\infty g(x)\,\frac{\sinh x}{\cosh x - 1 + \tfrac12\,|1-\epsilon_i^2|^2}\,dx
\]

(see thesis text around (4.3.3); same structure in the global formula with
\(\alpha_{\alpha k}\)).

For a **cuspidal elliptic of order 3** in \(\mathrm{PSL}(2,\mathbb Z[\omega])\),
the rotation multiplier satisfies \(\epsilon^2\) a primitive cube root of unity
in the appropriate sense of the setup, and the code (and EGM/Friedman data
for this field) use

\[
|1-\epsilon^2|^2 = 3.
\]

Hence

\[
\cosh x - 1 + \tfrac12\cdot 3 = \cosh x + \tfrac12,
\]

so the kernel is exactly \(\sinh x/(\cosh x+\tfrac12)\). This is **not** an
accidental replacement of \(\tanh(x/2)\); it is the same formula that for
order-2 Picard data (\(|1-\epsilon^2|^2=4\)) becomes

\[
\cosh x - 1 + 2 = \cosh x + 1,
\]

which is the usual \(\tanh(x/2)\) integrand after a hyperbolic identity.

**Paper requirement:** Quote Lemma 4.3.2 (or EGM equivalent) and record the
one-line specialization \(|1-\epsilon^2|^2=3\Rightarrow c_{\mathrm{kern}}=\tfrac12\).
Do **not** cite only `ckern=1/2` in `bianchi_omega.py`.

**Software status:** Hardcoded consistently in `bianchi_omega.py` /
`bianchi_omega_arb.py` (`ckern=1/2` for omega). Integration itself is Arb-
quadrature; the kernel choice is mathematical input.

---

## Q2. Is “exactly six CE classes” proved, or only enumerated?

**Answer: hybrid — algebraic upper bound + finite check + Friedman identity.
Not “the code found six.” Not a free-standing classification theorem in
our prose yet.**

`old_RIGOR_GAPS.md` §2 (reproduced in spirit):

1. **Upper bound (algebra).** Class number one ⇒ every CE is conjugate to
   \(g_{u,b}=\bigl(\begin{smallmatrix}u&b\\0&u^{-1}\end{smallmatrix}\bigr)\)
   with \(u\in\mathbb Z[\omega]^\times\), \(u^2\neq 1\). Mod \(\{\pm I\}\) the
   nontrivial multipliers are \(u=\omega,\omega^2\). Translations act on \(b\)
   by shifts of \((u^{-1}-u)\); that ideal has norm 3, so **at most**
   \(2\cdot 3=6\) residue candidates.

2. **Finite exact arithmetic.** For each candidate, reduce denominators and
   compute centralizer order (Q3). Code: `exact_eisenstein_CE_check`.

3. **Identity (Friedman).** With \(k_\infty=l_\infty=1\),
   \([\Gamma_\infty:\Gamma'_\infty]=3\), \(|1-\epsilon^2|^2=3\), Lemma 4.4.4
   (as cited in RIGOR_GAPS) gives
   \[
   2\sum_i\frac1{3\,|C(g_i)|}+\frac13=1.
   \]
   If every actual class has \(|C|=6\), each contributes weight \(1/18\), so
   the sum has **exactly six** terms.

4. **Conclusion.** At most six candidates, each realisable weight matches,
   identity saturates ⇒ those six are distinct and exhaustive.

**Where “forces” is load-bearing:**

- The identity must be Friedman’s **after** the CE classes enter with those
  weights — not an identity that *assumes* you already listed six classes
  and then “checks” them. RIGOR_GAPS cites Lemma 4.4.4 with fixed
  \(k,l,[\Gamma:\Gamma']\); **paper must quote that lemma’s hypotheses**.
- Completeness of the upper bound (only \(u=\omega,\omega^2\), only those
  residue modules) is representation theory / class number one, not software.

**Software status:** `exact_eisenstein_CE_check` proves the **integer**
arithmetic (six candidates, \(|C|=6\), multiset of \(|c|^2\), weights sum to
\(1/3\)). Exhaustiveness is **mathematical**, documented in RIGOR_GAPS, not
re-proved by a new paper-quality writeup in this directory.

**Risk flag (reviewer is right):** If Lemma 4.4.4 is misread as “assuming the
list,” the argument is circular. Verify the lemma statement in Friedman
before using the word **forces** in a publication.

---

## Q3. Why is \(|C(g)|=6\) for every candidate?

**Answer: finite exact computation from an algebraic criterion — not an
abstract “always 6” theorem independent of representatives.**

Criterion (RIGOR_GAPS §2; code `centralizer_order`): an element centralizing
\(g_{u,b}\) that does not swap the two fixed points is upper triangular with
diagonal unit \(p\in O_K^\times\), and

\[
b\frac{p-p^{-1}}{u-u^{-1}}\in O_K.
\]

There are **six** units. Loop over them; count how many satisfy the
divisibility test (`OK.divides`, exact integer arithmetic). For all six
\((u,b)\) residue candidates the count is 6.

**Strengths:** deterministic, no float, reproduces on Picard (\(|C|=4\)) in
`run_CE("i")`.

**Weaknesses for a paper:**

- Does not replace a conceptual reason (e.g. full unit group of the
  centralizer torus equals \(O_K^\times\)).
- Relies on “cannot interchange fixed points” (order-3 vs inverse) — needs
  one sentence in the writeup.
- SL vs PSL conventions for \(|C|\) must match Friedman’s normalization
  (code comments say centralizer in SL).

**Software status:** Asserted in `exact_eisenstein_CE_check` (`assert C==6`).

---

## Q4. Are \(\frac29\log 3\) and \(\frac13\) symbolic or numerical?

**Answer: symbolic consequences of the multiset of \(|c|\) and \(|C|\), once
that multiset is accepted.**

From the six candidates, exact reduction gives (code counters)

\[
|c_i|^2 \in \{1,1,3,3,3,3\},\qquad |C(g_i)|=6,\qquad |1-\epsilon^2|^2=3.
\]

Weight per class \(w_i=1/(6\cdot 3)=1/18\):

\[
\sum_i w_i = 6\cdot\frac1{18}=\frac13,
\]

\[
\sum_i 2\log|c_i|\,w_i
= 2\cdot\frac1{18}\Bigl(2\log 1 + 4\log\sqrt3\Bigr)
= \frac4{18}\cdot\tfrac12\log 3
= \frac29\log 3.
\]

So **if** the multiset is proved (Q2–Q3), the CE coefficients in
`bianchi_omega*.py` (`CEg0=(2/9)log3`, `CEint=1/3`) are **exact**, not
Arb-fitted. Arb is used only to integrate the continuous kernel against \(g\).

**Software status:** `Fraction` checks `weights==1/3` and
`c2=={1:2,3:4}`; float `run_CE` matches Picard target \((5/16)\log2\), \(1/4\).

---

## Regression plan (sensible; keep it)

1. Picard CE: `run_CE("i")` → \((5/16)\log2\), \(1/4\), identity \(=1\).
2. Exact Eisenstein finite algebra: `exact_eisenstein_CE_check`.
3. Friedman identity numeric/symbolic saturation.
4. Only then: Arb CE integral in `bianchi_omega_arb.py`.

This matches the Gaussian discipline. **Do not** declare the ℤ[ω] engine
“done” until Q1–Q4 appear as propositions with citations.

---

## One-sentence status for collaborators

| Claim | Verdict |
|-------|---------|
| Kernel \(\sinh/(\cosh+\frac12)\) | **OK if** specialized from Friedman 4.3.2 + \(|1-\epsilon^2|^2=3\) written out |
| Exactly six CE classes | **Conditional OK** — upper bound + \(|C|=6\) + Lemma 4.4.4; **prove non-circularity of the identity** |
| \(|C|=6\) always (for these) | **Computational exact check**, not abstract; acceptable if criterion is stated |
| \(\frac29\log3\), \(\frac13\) | **Symbolic** from multiset, not float inference |

**Cautiously optimistic on architecture; not publication-complete on CE
classification prose.** Same standard as the ℤ[i] congruence ladder: every
load-bearing constant needs a named lemma, not only a green assert.

---

## Before treating CE as frozen for a paper

- [ ] Quote Friedman 4.3.2 and derive \(c_{\mathrm{kern}}=\frac12\).
- [ ] Quote Friedman 4.4.4 (or EGM) and write the six-class proposition
      without the verb “forces” unless the lemma’s logic is unpacked.
- [ ] State the centralizer criterion and that the unit loop is exhaustive.
- [ ] Derive \(\frac29\log3\) and \(\frac13\) from the multiset in the text.
- [ ] Keep Picard regression as a permanent gate in CI / smoke tests.

*Preliminary sandbox: `eisenstein_numbers/`. Production CE code: repo root
`cuspidal_ce.py`, `bianchi_omega_arb.py`. Gaps note: `old_RIGOR_GAPS.md`.*
