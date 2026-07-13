# C1 §1.6 — **ACCEPTED**: drop \((1+\sqrt\lambda)\) from \(A_{\mathrm{bdry}}\)

**Status:** **proved by bookkeeping rewrite** in `theorem_DK.tex` (Remark `rem:no-double-sqrt`).  
**Production default:** `bdry_sqrt_lam=False` in `C1_C2_constants`.  
**Hard map:** unchanged (residual still blocks certification).

---

## Verdict

| Question | Answer |
|----------|--------|
| Is \((1+\sqrt\lambda)\) needed in Step 1 \(H^1\) corrector? | **No** — Step 1 is λ-free: \(\|e\|_{H^1}\le C_{\mathrm{tr}}C_{\mathrm{Sob}}A_{\mathrm{met}}\delta\) |
| Is λ-dependence needed for \((\Delta-\lambda)e\)? | **Yes, once** — lives in \(A_{\mathrm{ell}}=1+\lambda\subset A_{\mathrm{res}}\) (Step 2) |
| Keeping both factors | **Double count** — product \(\propto(1+\sqrt\lambda)(1+\lambda)\approx\lambda^{3/2}\) |
| Drop from \(A_{\mathrm{bdry}}\) only | **Safe** — no estimate weakened; factor gain \(1+\sqrt\lambda\) |

---

## Formulae (current production)

\[
A_{\mathrm{bdry}}
= C_{\mathrm{tr}}\,C_{\mathrm{Sob}}\,A_{\mathrm{met}},
\qquad
A_{\mathrm{ell}}=1+\lambda,
\qquad
A_{\mathrm{res}}=A_{\mathrm{ell}}(1+C_P A_{\mathrm{met}}),
\qquad
C_1 = 2\,A_{\mathrm{bdry}}\,A_{\mathrm{res}}\,(1+A_{\mathrm{cusp}}).
\]

Legacy double-count (regression only): `bdry_sqrt_lam=True` multiplies \(A_{\mathrm{bdry}}\) by \((1+\sqrt\lambda)\).

---

## Measured factor (sharp geom, r=6, Y=1.25)

| Mode | C1 | η₀ | factor |
|------|---:|---:|------:|
| legacy `bdry_sqrt_lam=True` | 1.4027e4 | 1.2706e-9 | 1× |
| **production** `False` | **1.9805e3** | **6.3739e-8** | **7.083×** |

\(1+\sqrt{37}\approx 7.082763\) matches C1 ratio.

**Orders @ η=5e-4:** legacy sharp gap ~5.6 → **~3.9** with §1.6.

Reproduce:

```bash
python c1_audit_16.py
python -c "from lemma_K import C1_C2_constants, ball_mid
for flag in (False, True):
 d=C1_C2_constants(field='i', sharp_geom=True, r=6.0, bdry_sqrt_lam=flag)
 print('bdry_sqrt_lam', flag, 'C1', ball_mid(d['C1']), 'eta0', ball_mid(d['eta0']))"
```
