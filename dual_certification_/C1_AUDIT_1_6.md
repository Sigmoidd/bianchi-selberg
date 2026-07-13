# C1 audit §1.6 — \((1+\sqrt\lambda)\) in \(A_{\mathrm{bdry}}\) vs \(A_{\mathrm{ell}}=1+\lambda\)

**Status:** diagnostic only. Production default keeps `bdry_sqrt_lam=True` (theorem as written).  
**API:** `C1_C2_constants(..., bdry_sqrt_lam=False)` drops the factor for factor measurement.  
**Hard map:** unchanged.

---

## Formulae

Paper / code (default):

\[
A_{\mathrm{bdry}}
= C_{\mathrm{tr}}\,C_{\mathrm{Sob}}\,A_{\mathrm{met}}\,(1+\sqrt\lambda),
\qquad
A_{\mathrm{ell}}=1+\lambda,
\qquad
A_{\mathrm{res}}=A_{\mathrm{ell}}(1+C_P A_{\mathrm{met}}),
\qquad
C_1 = 2\,A_{\mathrm{bdry}}\,A_{\mathrm{res}}\,(1+A_{\mathrm{cusp}}).
\]

So \(\lambda\)-dependence in \(C_1\) is at least

\[
C_1 \propto (1+\sqrt\lambda)\,(1+\lambda)
\approx \lambda^{3/2}
\quad(\lambda\gg 1).
\]

---

## Where each factor comes from (sketch)

| Factor | Role in defect argument |
|--------|-------------------------|
| \(C_{\mathrm{tr}}\) | Trace \(H^1(T)\to L^2(F)\) on tet faces |
| \(C_{\mathrm{Sob}}\) | \(H^1\to L^\infty\) / Poincaré chart on core |
| \(A_{\mathrm{met}}\) | Hyperbolic vs Euclidean metric comparison at \(y\ge y_{\min}\) |
| \((1+\sqrt\lambda)\) in \(A_{\mathrm{bdry}}\) | Gradient / elliptic energy control when lifting face jumps to an \(H^1\) corrector (Child / BSV style) |
| \(A_{\mathrm{ell}}=1+\lambda\) | Residual of \((\Delta+\lambda)\) on the periodized trial after corrector |

**Possible double count:** if the corrector estimate already uses \(\|(\Delta+\lambda)u\|\) (or an energy identity that absorbs \(\sqrt\lambda\) into the residual channel), then putting both \((1+\sqrt\lambda)\) and \(1+\lambda\) is conservative but not sharp.

**What would justify dropping \((1+\sqrt\lambda)\):** a rewritten proof that the boundary corrector satisfies

\[
\|e_{\mathrm{bdry}}\|_{H^1}
\le C_{\mathrm{tr}}C_{\mathrm{Sob}}A_{\mathrm{met}}\,\delta_{\mathrm{aut}}
\]

with **no** \(\sqrt\lambda\) (or with \(\sqrt\lambda\) only inside a residual already counted in \(A_{\mathrm{res}}\)). That rewrite is **not** done here.

---

## Measured factor (sharp geom, r=6, Y=1.25)

| Mode | C1 | η₀ | factor vs default |
|------|---:|---:|------------------:|
| `bdry_sqrt_lam=True` (default) | 1.4027e4 | 1.2706e-9 | 1× |
| `bdry_sqrt_lam=False` (audit) | 1.9805e3 | 6.3739e-8 | **7.083×** on C1 |

At \(r=6\): \(1+\sqrt\lambda = 1+\sqrt{37}\approx 7.082763\) — matches the C1 ratio to 3 digits.

---

## Decision

| Choice | When |
|--------|------|
| Keep default `True` | Until theorem_DK.tex proof drops \((1+\sqrt\lambda)\) with a complete estimate |
| Use `False` | Sensitivity / optimistic C1 only; **do not** claim certified smaller η₀ |

**Orders impact (η=5e-4):** default gap ~5.6 orders; audit-only C1 would give ~4.8 orders — still residual-dominated, not geometry-dominated.

Reproduce:

```bash
python -c "from lemma_K import C1_C2_constants, ball_mid
for flag in (True, False):
 d=C1_C2_constants(field='i', sharp_geom=True, r=6.0, bdry_sqrt_lam=flag)
 print(flag, ball_mid(d['C1']), ball_mid(d['eta0']))"
```
