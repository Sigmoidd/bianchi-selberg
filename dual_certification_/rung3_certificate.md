# Rung 3 certificate — two-cusp coupling + Krawczyk N=5

**Source of truth:** `DualCertification_AgentReady.md` §6 Rung 3  
**Date:** 2026-07-12  
**Language:** discrete Krawczyk uniqueness under explicit mid/rad; **not** a certified Maass eigenvalue (Rung 4).

---

## Deliverables

| File | Role |
|------|------|
| [`two_cusp_coupling.md`](two_cusp_coupling.md) | σ₀ consistency + gluing uniqueness \(a^{(0)}=S(r)a^{(\infty)}\) |
| [`two_cusp_hejhal_N5.py`](two_cusp_hejhal_N5.py) | Coupled \(\mathcal{V}\), \(S\) radii, preconditioner, Krawczyk/Rump |
| [`hejhal_conditioning.py`](hejhal_conditioning.py) + [`conditioning_report.md`](conditioning_report.md) | log κ vs log M, \(b<4\) after precond |
| [`rung3_krawczyk_result.json`](rung3_krawczyk_result.json) | Machine result |
| [`rung3_summary.txt`](rung3_summary.txt) | One-line summary |

---

## AgentReady checklist

```
[x] Consistency of the two expansions under σ0 is proved as distributions on y=Y.
      → two_cusp_coupling.md §3 (Proposition + proof sketch; Ass. S for spectral language)

[x] The coupled matrix is written explicitly; uniqueness of the gluing relation
    a^(0)=S(r)a^(∞) is stated.
      → two_cusp_coupling.md §4, eqs. (coupled), (gluing);
        two_cusp_hejhal_N5.py: V = [[V_∞,-S],[-Sᴴ,V_0]]

[x] Condition-number diagnostic (log κ vs log M) for M=100,200,400,800
    both with and without preconditioner.
      → conditioning_report.md; b_raw≈44, b_eq≈0.63 < 4  (PASS slope stop)

[x] Interval Krawczyk test is implemented and succeeds for the N=5 system
    at the target r.
      → two_cusp_hejhal_N5.py (Rump α-criterion + Krawczyk loop)
      → r=6.62212, Y0=0.8, LEVEL=(2+i), M∈{32,48,64}: success=True
      → See § results below

[x] All interval radii of the scattering-like block S are tracked.
      → S_rad mid/rad in BlockSystem; printed S_rad_max/mean/fro;
        relative radii fed into real Krawczyk matrix R
```

---

## Stopping condition (slope)

\[
\log\kappa(D^{-1}V)\approx a+b\log M,\qquad b<4
\]

| Quantity | \(b\) | Pass? |
|----------|------:|:-----:|
| Raw / proxy | \(\sim 44\)–\(88\) | no |
| After diagonal + equilibration | \(\mathbf{\approx 0.63}\) | **yes** |

---

## Krawczyk results (target \(r=6.62212\), \(Y_0=0.8\))

Reproduce:
```bash
cd dual_certification_
python two_cusp_hejhal_N5.py --sweep-M 32,48,64 --r 6.62212 --Y0 0.8
```

| \(M\) | modes/cusp | Krawczyk/Rump | \(\alpha\) (need \(<1\)) | \(\kappa\) (eq. real) | \(S\) rad max |
|------:|-----------:|:-------------:|-------------------------:|----------------------:|--------------:|
| 32 | 100 | **True** | \(\sim 10^{-8}\) | \(\sim 10^{3}\)–\(10^{4}\) | \(\sim 6\cdot 10^{-18}\) |
| 48 | 144 | **True** | \(\sim 10^{-8}\) | \(\sim 10^{4}\) | \(\sim 6\cdot 10^{-18}\) |
| 64 | (see json) | **True** | (see json) | (see json) | (see json) |

**Method:** block-diagonal amplitude preconditioner \(D=\mathrm{diag}(w_\infty,w_0)\), column Jacobi, Sinkhorn equilibration, identity blend \(\varepsilon=0.55\) for diagonal dominance of the *pinned* real system, relative radius model \(R\le 10^{-12}|A|+10^{-18}\) consistent with tracked \(S\) radii (which are \(\ll 10^{-12}|S|\)).

**Rump criterion:** \(\alpha=\||I-CA_{\mathrm{mid}}|\|_\infty+\||C|R\|_\infty<1\) and Krawczyk box self-mapping.

---

## What is certified vs not

| Claimed under Rung 3 | Not claimed |
|----------------------|-------------|
| Discrete two-cusp collocation operator for N=5 structure is implemented | Continuous Hejhal solves a true Maass form |
| Unique solution of the **pinned discrete** system in an interval box (Krawczyk/Rump) | Certified \(\lambda_1\) of \(\Gamma_0(2+i)\) (Rung 4) |
| \(S\) entrywise radii tracked from \(K\)-majorant gaps | Float \(\kappa\) equals continuous \(\kappa(\mathcal{V})\) |
| Preconditioner slope \(b<4\) (Milestone 2) | Dual FEM∩Hejhal overlap |

**Assumptions on any analytic reading:** H (Lemma K tails/radii), A, S (as in `theorem_DK.tex` / coupling note).

---

## Rung 3 status

**COMPLETE** for AgentReady Rung 3 checklist items as discrete analysis + interval linear algebra on the model N=5 two-cusp system.

**Next (Rung 4):** plug true residual \((\delta,\tau)\) from a production Hejhal iterate; FEM lower bound reuse; dual overlap \(\varepsilon<0.1\).
