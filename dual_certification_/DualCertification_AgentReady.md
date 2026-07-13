# Dual Certification Roadmap for \(\mathrm{Spec}_\Delta(M_{\mathfrak{p}})\)

**Version 5 — Agent-Executable Research Architecture**  
**Date:** 2026-07-12  
**Status:** Research program roadmap + complete certification checklist  
**Audience:** Human mathematician + autonomous research/coding agent

---

## 0. Executive Summary for an Agent

You are to execute a dual-certification pipeline that produces a rigorously certified first Laplace eigenvalue on Bianchi congruence orbifolds \(M_{\mathfrak{p}} = \Gamma_0(\mathfrak{p})\backslash\mathbb{H}^3\).

**Core principle.**  
A result is **certified** only when every numerical claim is accompanied by an explicit Arb enclosure (or Rump certificate), every analytic hypothesis is declared, every mesh/precision independence check passes, and the dual-verification overlap condition holds. Agreement between methods is persuasive; the mathematical proof still resides in each method satisfying its own hypotheses.

**Primary output.**  
A certified interval \([\lambda_1 - \varepsilon, \lambda_1 + \varepsilon]\) with \(\varepsilon < 10^{-4}\) (level 1) or \(\varepsilon < 0.1\) (first congruence levels), together with a counting certificate that it is the first eigenvalue, subject to the explicit hypotheses listed below.

**The larger picture is the dual loop, not any single rung.**

---

## 1. Contribution Statement (do not weaken)

The contribution is a certification architecture that combines independent variational (Crouzeix–Raviart + Rump + Arb) and automorphic (interval Hejhal + defect bound) methods into a single rigorous workflow for arithmetic hyperbolic orbifolds arising from Bianchi groups and their congruence subgroups.

A single method cannot close the bottom of \(\mathrm{Spec}_\Delta(M_{\mathfrak{p}})\). Certified exclusion of \((0,1)\) plus existence of an eigenvalue near \(45\) still leaves open holes at \(3,7,12,18,\dots\). The architecture supplies the missing piece — certified counting — while making every logical dependency and every unproved analytic hypothesis explicit.

Independent certification leads to agreement; agreement leads to high confidence; the proof still comes from each method satisfying its own hypotheses. Agreement is not itself a mathematical proof.

---

## 2. Separation of Concerns

**Known Mathematics**  
Friedman Thm 3.8.1 (spectral decomposition), Shimizu horoball lemma, Lax–Phillips cusp reduction, Payne–Weinberger Poincaré constant, Crouzeix–Raviart interpolation, Rump positive-definiteness verification (BIT 46, 2006), Selberg trace formula as identity.

**Engineering already realized**  
Arb enclosures of all geometric constants \(h_T,\gamma,\tau,S_M,S_Q,S_\Sigma,V_\Sigma\); \(n\sim 10^4\)–\(3\cdot 10^4\) dofs; residue rings \(\mathbb{Z}[i]/\mathfrak{p}\cong\mathbb{F}_N\); reference-cell gluing of \(N\mathfrak{p}+1\) copies; Rump shifts \(\sim 2\cdot 10^{-8}\).

**Unknown Theorems (explicit research objectives)**  
1. Explicit defect bounds for Bianchi groups (single-cusp and two-cusp) — Theorem D(\(K\)).  
2. Uniform \(K_{ir}\) tails with Hecke-controlled coefficients — Lemma K.  
3. Conditioning of the coupled two-cusp Hejhal system.  
4. Certified eigenvalue counting via two independent routes.  
5. Scaling of interval widths under reference-cell tiling.  
6. Uniform spectral-gap lower bound \(\mu(N)\) — heuristic only; never claimed.

---

## 3. Explicit Analytic Hypotheses (must appear in every certificate)

**Assumption H (Hecke-type growth).**  
There exist constants \(C_\varepsilon = C_\varepsilon(\Gamma,\chi)>0\) such that the Fourier coefficients of any Maass form of eigenvalue \(\lambda=1+r^2\) satisfy
\[
|a_\beta| \le C_\varepsilon\, N(\beta)^{1/2+\varepsilon}
\]
for every \(\varepsilon>0\). All subsequent constants (especially the defect constant \(C(K,Y,r)\)) depend on \(C_\varepsilon\) and must track that dependence. Never hide this assumption.

**Assumption A (analyticity / unique continuation).**  
Eigenfunctions are real-analytic; unique continuation holds across the compact core (Morrey–Nirenberg / Aronszajn).

**Assumption S (spectral decomposition).**  
Friedman Thm 3.8.1 applies to the congruence subgroups under consideration (trivial character).

Any certificate that fails to list the hypotheses it relies on is incomplete.

---

## 4. Dependency Graph

```text
Paper I–III (existing FEM exclusion)
      │
      ▼
Certified exclusion of (0,1)
already done for N𝔭 = 1,5,9,13 and Eisenstein–Picard

      │
      ├──────────────────────┐
      ▼                      ▼
Route A                  Route B
Dirichlet–Neumann        Lehmann–Goerisch +
bracketing               interval Hejhal
counting                 (defect + Krawczyk)
      │                      │
      └──────────┬───────────┘
                 ▼
          Dual certification
          FEM interval [a,b] ∩ Hejhal interval [c,d] ≠ ∅
          + counting certificate that N(λ) = 0 on (1,λ)
                 │
                 ▼
          Rigorous first eigenvalue
          (subject to counting, conditioning,
           implementation, and analytic hypotheses)
                 │
                 ▼
          Congruence levels N=5,9,13,\ldots
          (where present trace-formula B > 1)
```

Route A and Route B are **not** interchangeable; they rest on different hypotheses. The architecture is robust precisely because both are pursued.

---

## 5. Dual Certification Loop (the persuasive core)

```
FEM certifies          λ₁ ∈ [a,b]
                       lower bound from Theorem G1 / CR
                       upper bound from conforming Rayleigh quotient
                       — engineering target for window near 45

Hejhal independently   λ₁ ∈ [c,d]
                       defect bound of Child/BSV type
                       with explicit C(K,Y,r) from Theorem D(K)

If [a,b] ∩ [c,d] ≠ ∅ and diam([a,b] ∪ [c,d]) < tolerance
then the two intervals mutually verify each other.
```

**Language rule (mandatory).**  
- Use “engineering target” for any FEM upper bound near 45 until a certified computation exists.  
- Use “addresses the remaining logical gap” for counting until a counting theorem is proved.  
- Never write “only our pipeline can”; write “at levels where the present trace-formula exclusion criterion no longer closes (B ≳ 1) the dual-certification pipeline becomes the proposed route to a rigorous first eigenvalue.”

---

## 6. Challenge Ladder with Deliverables, Dependencies, Stopping Conditions & Full Certification Checklists

### Rung 0 — Lemma K (K-Bessel + Hecke growth)  
**Status:** free win, highest priority  
**Deliverable:** `lemma_K.py` + `lemma_K_certificate.md`  
**Dependency:** none  
**Stopping condition (must all hold):**  
- Arb enclosure of the tail majorant \(\varepsilon(M,Y_0,r,\theta) < 10^{-30}\) at \(M=100\), \(Y_0=0.8\), \(r=6.62212\).  
- Ratio (enclosure / truncated sum up to \(n\le 20000\)) \(< 10^3\).  
- Self-test against Luke double inequality passes for 20 random \((r,y)\) pairs.  
- Exact \(r_2\) for \(\mathbb{Z}[i]\) used for the benchmark.

**Agent certification checklist for Rung 0**  
```
[ ] Assumption H is stated and C_ε (or θ) is an explicit input parameter.
[ ] tail_majorant(M, Y0, r, theta) returns an Arb ball.
[ ] For M=100,200,400 the enclosure is < 10^{-30} (or tighter).
[ ] Comparison against direct truncated sum (n ≤ 20000) using Then’s first Fourier coefficients shows the majorant is not more than 10^3 times larger.
[ ] Luke-type double inequality is recovered as a special case (or validated numerically).
[ ] All special-function evaluations (K_ir, gamma, etc.) are performed in Arb.
[ ] Floating-point backend is used only for diagnostics; final certificate is pure Arb.
[ ] File lemma_K_certificate.md records every constant, every precision used, and the final enclosure.
```

---

### Rung 1 — Single-cusp defect bound, Theorem D(K)  
**Deliverable:** `theorem_DK.tex` + `constants_DK.md` + `defect_bound_arb.py`  
**Dependency:** Rung 0  
**Stopping condition:**  
- \(C_1 \sim 10^5\)–\(10^6\) at the concrete point \((Y,Y_0,r)=(1.25,0.8,6.6)\).  
- \(C_2\) fully explicit.  
- Lemma K tail at \(M=400\) is \(< 10^{-70}\) and therefore negligible compared with any realistic automorphy defect \(\delta\).  
- No \(r\)-independent lower bound is claimed.

**Agent certification checklist for Rung 1**  
```
[ ] Theorem statement lists Assumptions H, A, S and the geometric parameters (Y, Y0, r, K).
[ ] The constant C(K,Y,r) is an explicit expression (or Arb-computable function) of known quantities only.
[ ] Trace inequality on the boundary torus (Lemma T or equivalent) is proved with explicit constant.
[ ] Weighted elliptic regularity / Sobolev estimates on {y ≥ y_min} are explicit.
[ ] Automorphy defect δ and Fourier tail τ appear only through the linear combination C·(δ+τ).
[ ] Numerical evaluation of C at the target point (Y=1.25, Y0=0.8, r=6.622) is performed in Arb and recorded in constants_DK.md.
[ ] Sensitivity of C to small changes in Y0 and r is tabulated.
[ ] The final defect theorem is written so that any later Hejhal run can simply plug in its measured (δ,τ) and obtain a certified |λ̃−λ| enclosure.
```

---

### Rung 2 — Level-1 dual certification (engineering target)  
**Deliverable:** dual interval + counting certificate for Picard level 1  
**Dependency:** Rungs 0 and 1  
**Stopping condition:**  
- FEM interval \([a,b]\) and Hejhal interval \([c,d]\) satisfy \([a,b]\cap[c,d]\ne\emptyset\) and \(\mathrm{diam}([a,b]\cup[c,d]) < 10^{-4}\).  
- Counting (Route A or B) certifies \(N(44.8)=0\).  
- Final certified first eigenvalue is the intersection (or a rigorous hull) of the two intervals.

**Agent certification checklist for Rung 2**  
```
[ ] FEM lower bound uses the already-certified Theorem G1 / CR machinery (exclusion of (0,1) is free).
[ ] FEM upper bound (conforming Rayleigh) is either
      (i) a true certified upper bound, or
      (ii) explicitly labelled “engineering target / non-certified diagnostic”.
[ ] Interval Hejhal run uses M=400, Y0=0.8, 128-bit (or higher) Arb, and the defect bound from Theorem D(K).
[ ] The measured automorphy defect δ and tail τ are recorded; the defect theorem is applied.
[ ] Krawczyk (or equivalent) interval solver for the Hejhal linear system returns a unique solution enclosure.
[ ] Overlap condition is verified by rigorous interval arithmetic (not float comparison).
[ ] Counting certificate (Route A or B) is present and its own checklist (below) is passed.
[ ] All floating-point intermediate results are guarded by outward rounding or replaced by Arb.
[ ] Reproducibility suite (Milestone 1) has been run on this certificate.
```

---

### Rung 3 — Two-cusp coupling analysis  
**Deliverable:** well-posedness proof + condition-number bound + interval Krawczyk for N=5  
**Dependency:** Rung 1  
**Stopping condition:**  
After a diagonal (or block-diagonal) preconditioner one has
\[
\log\kappa(D^{-1}V) \approx a + b\log M \qquad\text{with } b < 4
\]
(current naïve growth is \(b\sim 76\) at \(Y_0=0.8\)).

**Agent certification checklist for Rung 3**  
```
[ ] Consistency of the two expansions under σ0 is proved as distributions on y=Y.
[ ] The coupled matrix is written explicitly; uniqueness of the gluing relation a^(0)=S(r)a^(∞) is stated.
[ ] Condition-number diagnostic (log κ vs log M) is produced for M=100,200,400,800 both with and without preconditioner.
[ ] Interval Krawczyk test is implemented and succeeds for the N=5 system at the target r.
[ ] All interval radii of the scattering-like block S are tracked.
```

---

### Rung 4 — N=5 dual certification (first congruence level where B ≳ 1)  
**Deliverable:** certified interval for \(\Gamma_0(2+i)\)  
**Dependency:** Rungs 2 and 3  
**Stopping condition:** certified interval \([\lambda_1-\varepsilon,\lambda_1+\varepsilon]\) with \(\varepsilon<0.1\), disjoint from (0,1), and counting certifies it is the first eigenvalue.

**Agent certification checklist for Rung 4**  
```
[ ] FEM lower bound already certified (μ≈4.4) is re-used or re-verified.
[ ] Two-cusp Hejhal uses the coupling analysis of Rung 3.
[ ] Overlap + counting succeed exactly as in Rung 2, with the looser tolerance ε<0.1.
[ ] Residue-ring arithmetic for 𝔽5 is verified by an independent unit-test suite.
[ ] Gluing of the 6 copies is exact (combinatorial identity, not numerical).
```

---

### Rung 5 — Scaling theorem / safety of the jump  
**Deliverable:** numerical evidence + (if possible) theorem that \(W(N)=O(1)\) or \(O(\log N)\)  
**Stopping condition:** plot of max entry-wise radius of the glued interval matrix vs N is essentially flat (supports O(1)) or grows at most logarithmically; κ diagnostic remains stable.

**Agent certification checklist for Rung 5**  
```
[ ] congruence_prototype.py (or equivalent) run for N=5,9,13,25,37.
[ ] Max radius, mean radius, and 99-percentile radius of Q and M are recorded.
[ ] κ(V) (or the float pencil) is plotted against N.
[ ] Reference-cell principle (isometric copies + bounded-degree gluing) is used to explain the observed growth.
[ ] Conjecture A (O(1)) is stated as primary; O(log N) as fallback; neither is claimed as proved unless a theorem is written.
```

---

### Rung 6 — Field independence  
**Deliverable:** port of the entire pipeline to \(\mathbb{Z}[\omega]\) and to inert primes (N=9)  
**Stopping condition:** same dual-certification success criteria as Rung 4, with the correct area |T∞|=√3/6 and rotation order 3.

---

## 7. Two Routes for Certified Counting — Full Checklists

### Route A — Dirichlet–Neumann Bracketing  
**When to use:** preferred first attempt; reuses existing FEM meshes.

**Agent checklist**  
```
[ ] Compact core K is truncated at a certified height Y with product collar structure.
[ ] Mesh on K is product-like near the artificial boundary (or the error of non-product structure is enclosed).
[ ] Dirichlet and Neumann eigenvalue problems are both solved with guaranteed bounds (CR or Morley + post-processing).
[ ] The bracketing inequality λ_k^{N,K_h} ≤ λ_k ≤ λ_k^{D,K_h} is proved with explicit constants that account for the truncation error.
[ ] The artificial-boundary error is itself enclosed by an Arb ball (or by a comparison with a larger domain).
[ ] Final counting function enclosure N(λ) is an integer interval; for the target λ it must be [0,0] on (1,λ₁).
```

### Route B — Lehmann–Goerisch + approximate eigenfunctions  
**When to use:** when high-quality trial functions from Hejhal become available.

**Agent checklist**  
```
[ ] A certified lower bound ρ on λ_{m+1} is available (from exclusion, ρ=1 is free).
[ ] Approximate eigenfunctions (from Hejhal) and flux reconstructions (hypercircle or residual equilibration) are supplied.
[ ] Residuals of the trial functions are small enough that the Lehmann–Goerisch hypotheses hold; the residual norms are Arb-enclosed.
[ ] The resulting lower bounds on λ_{m-n+1},\ldots,λ_m are rigorously larger than any previously excluded region.
[ ] Upper bounds come from the Rayleigh quotients of the same trial functions (or from conforming FEM).
[ ] Two-sided counting follows.
```

Both routes must be implemented; success of either is sufficient for the dual loop.

---

## 8. Required Benchmarks (Scientific Infrastructure)

### Milestone 1 — Reproducibility (elevated to mandatory)  
Take the published Picard value \(r_1\approx 6.62212\).  
Re-run the entire dual-certification pipeline under:

- completely independent mesh generators,  
- different quadrature rules,  
- different node numberings / dof orderings,  
- different compilers and floating-point libraries,  
- Arb precisions 64 / 128 / 256.

**Acceptance criterion:** the final certified interval is identical up to the declared enclosure radius.  
Any difference larger than the radius is a hard failure.

**Agent checklist**  
```
[ ] At least three independent meshes are generated.
[ ] At least two different quadrature schemes are used.
[ ] Node renumbering (random or reverse Cuthill–McKee) is tested.
[ ] Results at 64-bit, 128-bit and 256-bit Arb are compared.
[ ] A machine-readable log of all variations and the resulting intervals is produced.
[ ] The certificate file states “reproducibility suite passed”.
```

### Milestone 2 — Conditioning diagnostic  
Plot \(\log\kappa(V(r))\) vs \(\log M\) for \(M=100,200,400,800\) at fixed \(Y_0\).  
Current naïve growth \(\sim M^{76}\) is unacceptable. After diagonal preconditioner one must obtain \(b<4\).

### Milestone 3 — Width-versus-N diagnostic  
Free with existing prototypes; plot max radius of glued interval matrices.

### Milestone 4 — Cross-validation against Then  
Once dual certification succeeds, compare the rigorous enclosure with Then’s high-precision non-rigorous value. Agreement to many digits is strong corroboration (not a proof).

---

## 9. Heuristic Motivation for \(\mu(N)\) (Appendix — explicitly not claimed)

This section is labelled heuristic and forms no part of any proof or certificate.

As \(N\) grows the discrete operator is essentially \(N+1\) copies of the level-1 operator plus a gluing permutation. The first new eigenvalue is controlled by the spectral gap of the Schreier graph of \(\mathrm{PSL}_2(\mathbb{F}_N)/B\). These graphs are expanders (Lubotzky–Phillips–Sarnak, Bourgain–Gamburd); diameter \(O(\log N)\) yields, via Cheeger–Buser, a gap \(\gtrsim c/\log N\). This motivates the logarithmic form of the uniform-gap conjecture. It is never used in a certificate.

---

## 10. Risk Register (must be monitored by the agent)

| Risk | Severity | Mitigation / Agent Action |
|------|----------|---------------------------|
| Defect constant \(C(K,Y,r)\) blows up as \(r\to 0\) | High | Increase \(Y_0\); track C in Arb; abort if C > 10^8 |
| Hejhal conditioning exponential in \(M\) | High | Preconditioner mandatory; abort if b ≥ 4 after preconditioning |
| FEM upper bound never reaches ~45 | Medium | Accept coarser window; rely more on Hejhal + counting |
| Interval width grows with \(N\) | Medium | Width-vs-N plot; fall back to O(log N) claim |
| Route A counting error too large | Medium | Switch to Route B |
| Spectral clustering \(\lambda_2-\lambda_1 \approx 10^{-8}\) | Medium | Higher precision; block Lehmann–Goerisch; check Then spectrum gaps |
| Hidden Hecke dependence | Low | Assumption H must appear in every paper and every certificate |
| Reproducibility failure | Critical | Hard stop; do not publish any certificate that fails Milestone 1 |

---

## 11. Immediate 30-Day Experiment Queue (ordered)

| Pri | Experiment | Unlocks | Effort | Deliverable |
|-----|------------|---------|--------|-------------|
| 1 | Lemma K + `constants_DK.md` Arb tables | All later defect bounds | 3–5 d | `lemma_K.py`, certificate |
| 2 | Reproducibility suite (independent meshes) | Scientific credibility | 1 w | reproducibility log |
| 3 | Route A counting prototype on product mesh \(K_Y\) | Closes logical gap without Hejhal | 1–2 w | counting certificate for level 1 |
| 4 | Single-cusp interval Hejhal + diagonal preconditioner | Existence half of dual loop | 2 w | Hejhal interval + defect application |
| 5 | Width-vs-N and κ-vs-M plots (with preconditioning) | Safety of scaling | 2 d | diagnostic plots + decision on Conjecture A |

---

## 12. Language Discipline (agent must obey)

- “engineering target” until a certified upper bound near 45 exists.  
- “addresses the remaining logical gap” until a counting theorem is proved.  
- “Then \(N=5\) is the first congruence level at which the present trace-formula exclusion criterion no longer closes (B ≳ 1), making the dual-certification pipeline the proposed route to a rigorous first eigenvalue.”  
- Dual certification produces high confidence **subject to** the listed analytic, counting, conditioning and implementation hypotheses. The proof is the conjunction of those hypotheses being satisfied, not the mere numerical agreement of intervals.

---

## 13. Final Operational Definition of “Certified”

A numerical claim about an eigenvalue or a counting function is **certified** if and only if:

1. Every floating-point intermediate is either avoided or guarded by rigorous outward rounding / Arb.  
2. Every matrix positive-semidefiniteness claim is accompanied by a Rump certificate (or equivalent verified method).  
3. Every special-function evaluation that enters the final interval is performed in Arb.  
4. The dual-overlap condition (when claimed) is verified by interval arithmetic.  
5. The counting certificate (Route A or B) is present and its own checklist has passed.  
6. Assumption H (and any other analytic hypotheses) is explicitly listed.  
7. The reproducibility suite (Milestone 1) has been executed and has passed.  
8. All source files, mesh files, precision settings and random seeds are archived so that a third party can re-run the certificate.

If any item fails, the claim is **not certified**. The agent must refuse to emit a “certified” label.

---

*This document is the single source of truth for an autonomous agent executing the dual-certification research program. Every unproved analytic statement is deliberately labelled. The goal is not premature theorem claims, but precise, independent, attackable theorems together with complete, machine-checkable certification checklists. The larger picture remains the dual loop.*
