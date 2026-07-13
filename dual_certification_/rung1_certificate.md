# Rung 1 certificate — Theorem D(K) single-cusp defect bound

**Date:** 2026-07-11  
**Scope:** DualCertification_AgentReady.md §6 Rung 1 only  
**Artifacts:** `theorem_DK.tex`, `constants_DK.md`, `defect_bound_arb.py`, `lemma_K.py` (Rung 0, unchanged)  
**Backend:** python-flint Arb (`HAS_FLINT = True`)

---

## Language discipline

- Then \(r_1\approx 6.62212\) / \(\lambda_1\approx 44.85\) is an **engineering target**, not a certified eigenvalue of this rung.
- Outputs of `defect_to_lambda_error` are **enclosures of the defect bound** under Assumptions H, A, S — not a dual-certified \(\lambda_1\).
- Assumption H is an explicit input \((\theta,C_H,\varepsilon)\); never hidden.
- No \(r\)-independent lower bound on \(|K_{ir}|\) is claimed (sharp upper \(C_K=1\) only).
- No Hejhal solver, two-cusp coupling, or counting implemented on this rung.

---

## Stopping conditions

| Condition | Result | Evidence |
|-----------|--------|----------|
| \(C_1\sim 10^5\)–\(10^6\) at \((Y,Y_0,r)=(1.25,0.8,6.6)\) | **PASS** | \(C_1(\mathbb{Z}[i])=[1324375.37123066\pm 1.69\cdot 10^{-9}]\) **[Arb]** \(\sim 1.32\cdot 10^6\) |
| \(C_2\) fully explicit | **PASS** | \(C_2=A_{\mathrm{cut}}=2\|T\|Y^{-1}(1+C_{\mathrm{Sob}})A_{\mathrm{met}}\); \(\mathbb{Z}[i]\): \([11.7280139288508\pm 4.66\cdot 10^{-14}]\) |
| Lemma K tail \(M=400<10^{-70}\) | **PASS** | `[2.838307385819e-83 ± 6.11e-96]` **[Arb]** |
| No \(r\)-independent lower bound on \(\|K_{ir}\|\) | **PASS** | `lemma_K.C_K` default \(=1\) upper only; paper Lemma Kpoint |

Rung 0 dependency preserved: sharp \(C_K=1\), `lemma_K_tail` API, `C1_C2_constants` formulas unchanged.

---

## Agent checklist

```
[x] PASS  Theorem statement lists Assumptions H, A, S and geometric parameters (Y, Y0, r, K[, M]).
          Evidence: theorem_DK.tex Ass. H (ass:H), A (ass:A), S (ass:S);
          Theorem D(K) (thm:DK) opens with H+A+S and parameters Y≥Y0≥1/2, r≥0, M≥1, field K.

[x] PASS  C(K,Y,r) is an explicit Arb-computable function of known quantities only.
          Evidence: eq:Amet–eq:C2; lemma_K.C1_C2_constants; defect_bound_arb.defect_to_lambda_error.

[x] PASS  Trace inequality on the boundary torus (Lemma T) with explicit constant.
          Evidence: lem:T, C_tr=(4/3)|F|/|T|; C_trace^cons=4/h_min=16.

[x] PASS  Weighted elliptic / Sobolev on {y≥y_min} explicit.
          Evidence: lem:PW, lem:Sobolev, lem:metric; A_met, A_ell, A_res named.

[x] PASS  Automorphy defect δ and Fourier tail τ appear only through C·(δ+τ).
          Evidence: eq:lamerr = C1·(δ_aut+τ_tail)+C2 e^{-2π Y0}; API uses eta=delta+tau.

[x] PASS  Numerical evaluation of C at (Y=1.25, Y0=0.8, r=6.622) in Arb recorded in constants_DK.md.
          Evidence: constants_DK.md §10; C1(Z[i])≈1.3366e6, C2≈11.728 [Arb].

[x] PASS  Sensitivity of C to small changes in Y0 and r is tabulated.
          Evidence: constants_DK.md §11; python defect_bound_arb.py --sensitivity.

[x] PASS  Final defect theorem written so Hejhal can plug measured (δ,τ) and get |λ̃−λ| enclosure.
          Evidence: thm:DK plug-in paragraph; defect_bound_arb.defect_to_lambda_error;
          demo: python defect_bound_arb.py --demo.
```

**Checklist score: 8 / 8 PASS**

---

## Deliverables

| File | Role | Status |
|------|------|--------|
| `defect_bound_arb.py` | Main API `defect_to_lambda_error` + `--demo` + sensitivity | **done** |
| `constants_DK.md` | Arb tables + §10 (r=6.622) + §11 sensitivity + §12 API | **updated** |
| `theorem_DK.tex` | Assumptions A, S; Thm D lists H,A,S + Hejhal plug-in | **updated** |
| `rung1_certificate.md` | This checklist | **done** |
| `lemma_K.py` | Rung 0 (dependency; not modified for Rung 1 logic) | **intact** |

---

## Snapshot (reproduce)

```text
cd dual_certification_
python lemma_K.py --test --constants
python defect_bound_arb.py --demo
python defect_bound_arb.py --sensitivity
```

Representative Arb values at \((Y,Y_0,r,M)=(1.25,0.8,6.6,400)\), \(\theta=1/2\), \(C_H=1\), sharp \(C_K=1\):

| Quantity | \(\mathbb{Z}[i]\) | \(\mathbb{Z}[\omega]\) |
|----------|------------------:|----------------------:|
| \(C_1\) | \(1.324375\cdot 10^{6}\) | \(3.56428\cdot 10^{5}\) |
| \(C_2\) | \(11.7280\) | \(3.34568\) |
| \(\eta_0\) | \(\sim 1.43\cdot 10^{-13}\) | \(\sim 1.97\cdot 10^{-12}\) |
| tail \(M=400\) | \(\sim 2.84\cdot 10^{-83}\) | same |
| \(C_2 e^{-2\pi\cdot 0.8}\) | \(\approx 0.0770\) | \(\approx 0.0220\) |

**Toy Hejhal call** (`--demo` Case A: \(\delta=10^{-14}\), \(\tau=10^{-15}\)):

\[
|\lambda_{\mathrm{true}}-\lambda|
\le
C_1\eta+C_2 e^{-2\pi Y_0}
\approx
1.46\cdot 10^{-8}+0.0770
\approx 0.077
\]

The collar term \(C_2 e^{-2\pi Y_0}\) dominates until \(Y_0\) is raised or \(C_2\) is sharpened; automorphy quality controls the \(C_1\eta\) piece. This is the engineering bottleneck for a tight window near Then’s target — not the Lemma K tail.

---

## Explicit non-claims (Rung 1)

- No certified first eigenvalue.
- No interval Hejhal implementation / Krawczyk solve.
- No two-cusp block matrix.
- No counting certificate \(N(\lambda)=0\).
- No FEM dual overlap.

Those are Rungs 2+.

---

## Verdict

**Rung 1: PASS** (all stopping conditions and checklist items).
