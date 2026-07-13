# Dual Certification Roadmap for $\mathrm{Spec}_\Delta(M_{\mathfrak{p}})$

**Version 4 — Referee Proof Architecture Document**  
**Date:** 2026-07-12  
**Status:** Research program roadmap, not a completed theorem

The larger picture is the dual loop, not any single rung.

---

## Contribution Statement

The contribution is a certification architecture that combines independent variational and automorphic methods into a single rigorous workflow for arithmetic hyperbolic orbifolds arising from Bianchi groups and their congruence subgroups.

A single method cannot close the bottom of $\mathrm{Spec}_\Delta(M_{\mathfrak{p}})$. Certified exclusion of $(0,1)$ plus existence of an eigenvalue near $45$ still leaves open holes at $3,7,12,18,\dots$. The architecture supplies the missing piece, certified counting, while making every logical dependency and every unproved analytic hypothesis explicit.

```
Arithmetic group
    ↓ exact congruence tiling over O_K/p
    ↓ certified FEM (Lax Phillips + Crouzeix Raviart + Rump + Arb)
    ├── certified exclusion on (0,1)
    ├── certified upper bounds (engineering target)
    └── Route A / Route B counting
    ↓ interval Hejhal with two cusp coupling
    ↓ certified existence in [lambda - eps, lambda + eps]
    ↓ mutual verification (overlap)
    ↓ rigorous first eigenvalue, subject to each method satisfying its hypotheses
```

Independent certification leads to agreement, agreement leads to trust, but the proof still comes from each method satisfying its own hypotheses. Agreement between independent methods is very persuasive, it is not itself a mathematical proof.

---

## Separation of Concerns

**Known Mathematics**  
Spectral decomposition for cofinite Kleinian groups, Friedman Thm 3.8.1, Shimizu horoball lemma, Lax Phillips cusp reduction, Payne Weinberger Poincare constant, Crouzeix Raviart interpolation estimates, Rump positive definiteness verification, Selberg trace formula as an identity.

**Engineering, already realized in repository**  
Turning known inequalities into Arb enclosures with explicit constants $h_T,\gamma,\tau,S_M,S_Q,S_\Sigma,V_\Sigma$, managing $n\sim10^4$ to $3\cdot10^4$ degrees of freedom, residue ring arithmetic $\mathbb{Z}[i]/\mathfrak{p}\cong\mathbb{F}_N$ and $\mathbb{Z}[i]/(3)\cong\mathbb{F}_9$, gluing of $N\mathfrak{p}+1$ isometric copies of a reference cell, Rump certificates with shifts of order $10^{-8}$.

**Unknown Theorems, explicit research objectives**  
1. Explicit defect bounds for Bianchi groups, single cusp and two cusp.
2. Uniform $K_{ir}$ tails with Hecke controlled coefficients, Lemma K.
3. Conditioning of coupled two cusp Hejhal system.
4. Certified eigenvalue counting, two independent routes.
5. Scaling of interval widths under reference cell tiling.
6. Uniform spectral gap lower bound $\mu(N)$, heuristic only.

---

## Dependency Graph

This graph makes dependencies transparent and shows counting is not an independent oracle, it typically reuses FEM meshes or approximate eigenfunctions, and Hejhal is the only component currently absent from the repository.

```text
Paper I-III
      │
      ▼
Certified exclusion of (0,1)
already done for Np = 1,5,9,13 and Eisenstein

      │
      ├─────────────────┐
      ▼                 ▼
Route A             Route B
Dirichlet Neumann   Lehmann Goerisch + interval Hejhal
bracketing          defect bound + Krawczyk
counting            → existence [lambda - eps, lambda + eps]
      │                 │
      └──────┬──────────┘
             ▼
      Dual certification
      FEM interval [a,b] ∩ Hejhal interval [c,d] ≠ ∅
      diameter < tolerance, mutual verification

             ▼
      Rigorous first eigenvalue
      subject to counting, conditioning, implementation hypotheses

             ▼
      Congruence levels N=5,9,13,...
      where present trace formula exclusion no longer closes
```

Route A and B are not interchangeable, they have different hypotheses. That makes the roadmap more robust than pretending there is one magic counting theorem.

---

## Dual Certification Loop, the persuasive core

```
FEM certifies          lambda1 ∈ [a,b]
                       lower bound from Theorem G1 / CR,
                       upper bound from conforming Rayleigh quotient
                       — engineering target for window near 45

Hejhal independently   lambda1 ∈ [c,d]
                       defect bound of Child / BSV type
                       with explicit C(K,Y,r) from Theorem D(K)

If [a,b] ∩ [c,d] ≠ ∅ and diam([a,b] ∪ [c,d]) < tolerance
then the two intervals mutually verify each other.
```

Two mathematically distinct methods, one variational, one automorphic, certify the same eigenvalue. Overlap is far stronger evidence than either method alone. Non overlap immediately localizes the defect, mesh size vs Fourier tail vs conditioning.

**Important caveat on language.**  
The statement FEM upper bound $44 < \lambda_1 < 46$ is an engineering target, not an existing capability. The present CR machinery was engineered and certified for exclusion of $(0,1)$. Reaching a certified upper bound near $45$ requires $hp$ refinement or substantially larger mesh and is therefore listed as Rung 2 work, not as already available data. Using engineering target language matters, a referee will not accuse overselling of future work clearly identified as future work.

---

## Two Independent Routes for Certified Counting

Counting addresses the remaining logical gap identified by the reviewer, certifying no eigenvalues in $(0,1)$ plus an eigenvalue in $[44.8,44.9]$ does not exclude $3,7,12,18,25,39$ without a counting argument.

### Route A — Dirichlet Neumann Bracketing
* Natural for truncated domains that respect cusp collar product structure.
* Produces $\lambda_k^{N,K_h} \le \lambda_k \le \lambda_k^{D,K_h}$ by imposing Neumann or Dirichlet on artificial boundary of compact core $K$.
* Requires mesh product like in collar and guaranteed FEM bounds, already available from CR engine.
* Advantage, purely variational, reuses reference cell infrastructure.
* Risk, artificial boundary error must be enclosed.

### Route B — Lehmann Goerisch plus approximate eigenfunctions
* Uses certified lower bound $\rho$ for $\lambda_{m+1}$ from exclusion, $\rho=1$, to get lower bounds for $\lambda_{m-n+1}\dots\lambda_m$ via Lehmann Goerisch theorem.
* Uses Hejhal approximate eigenfunctions as trial functions for upper bounds.
* Advantage, can reuse Hejhal output, gives two sided $N(\lambda)$.
* Risk, needs small residual for trial functions.

Route A and B solve related but non identical problems. Both must be pursued, neither subsumes the other.

---

## Challenge Ladder with Deliverables, Dependencies, Stopping Conditions

**Rung 0 — Lemma K, K Bessel plus Hecke growth, free win**  
Deliverable: `lemma_K.py` with Arb enclosure `tail_majorant(M,Y0,r,theta)` and benchmark $M=100,200,400$ vs truncated sum $n\le20000$ using exact $r_2$ for $\mathbb{Z}[i]$. Stopping condition: tail enclosure $<10^{-30}$ at $M=100$, ratio enclosure to trunc $<10^3$, self test passes. Dependency: none. Status: implementation runs, float backend, needs python flint Arb wiring.

**Rung 1 — Single cusp defect bound, Theorem D(K)**  
Deliverable: `theorem_DK.tex` with Lemma K upper bound from integral representation, no $r$ independent lower bound claimed, $C_{\mathrm{tr}}$ consistent, metric comparison on $\{y\ge y_{\min}\}$ explicit, $\eta_0$ defined from a priori constants only. Stopping condition: $C_1\sim10^5$ to $10^6$ at $(Y,Y_0,r)=(1.25,0.8,6.6)$, $C_2$ explicit, Lemma K tail $<10^{-70}$ at $M=400$ negligible vs automorphy defect. Dependency: Rung 0. Status: draft complete, `constants_DK.md` missing.

**Rung 2 — Level 1 dual certification, engineering target**  
Deliverable: FEM exclusion $(0,1)$ already done, plus certified upper bound in window around $45$ from conforming method, engineering target, plus interval Hejhal at Then $r_1\approx6.62212$, $M=400$, $Y_0=0.8$, 128 bit Arb, plus Route A or B counting $N(44.8)=0$. Stopping condition: $[a,b]\cap[c,d]\neq\emptyset$ and diameter $<10^{-4}$, counting certifies no eigenvalue in $(1,44.8)$. Dependency: Rungs 0 and 1. Status: exclusion done, upper bound and Hejhal and counting are future work.

**Rung 3 — Two cusp coupling analysis**  
Deliverable: well posedness, consistency of $\sigma_0$ gluing, condition number bound for block matrix $\begin{pmatrix}V_\infty & -S\\ -S^* & V_0\end{pmatrix}$, implementation of interval Krawczyk for $N=5$. Stopping condition: $\log\kappa(D^{-1}V)\approx a+b\log M$ with $b<4$ after diagonal preconditioner, vs current naive $b\sim76$ at $Y_0=0.8$ showing steep growth. Dependency: Rung 1. Status: prompt created, prototype needed.

**Rung 4 — $N=5$ dual certification, first congruence level where present trace formula exclusion no longer closes**  
Wording softened from only your pipeline can. Deliverable: $\mathfrak{p}=(2+i)$, index 6, FEM lower bound already certified $\mu\approx4.4$, add Hejhal existence and counting, overlap gives mutual verification at first nontrivial congruence level where $B\approx1.9>1$. Stopping condition: certified interval $[\lambda_1-\varepsilon,\lambda_1+\varepsilon]$ with $\varepsilon<0.1$ disjoint from $(0,1)$ and counting certifies it is first. Dependency: Rungs 2 and 3.

**Rung 5 — Scaling theorem, safety of jump**  
Deliverable: prove or rigorously verify on sequence $N=5,9,13,25$ Conjecture A, width $W(N)=O(\log N)$ or $O(N^\alpha)$ with explicit $\alpha<1$, and $\kappa(V)=O(N^\alpha)$. Numerical evidence from `congruence_prototype.py gluing` and `kappa_vs_M` diagnostics required before claiming jump safe. Stopping condition: plot $\|Q_r\|$ vs $N$ flat, $b<4$ stable.

**Rung 6 — Field independence**  
Deliverable: port to $\mathbb{Z}[\omega]$ and inert primes $N=9$, area $|T|=\sqrt3/6$, rotation order 3, residue field $\mathbb{F}_9$, tests reference cell principle independent of ring.

**Explicitly left as conjectures, not claims, with heading explicitly not claimed**  
* Uniform $\mu(N)\ge c/\log N$
* Explicit Weyl remainder with constants for cusped Bianchi orbifolds
* Simplicity of $\lambda_1$

The heading explicitly not claimed does a tremendous amount of work, it tells the referee you know difference between motivation and proof. Never remove that heading.

---

## Required Benchmarks, scientific infrastructure

**Milestone 1 — Reproducibility, elevated to explicit milestone**  
Take published Picard eigenvalue $r_1\approx6.62212$. Re run entire dual certification pipeline with completely independent mesh generators, different quadrature rules, different node numberings, different compilers and floating point libraries, different Arb precisions 64, 128, 256. Final certified interval must be identical up to declared enclosure radius. This is software analogue of independent laboratory replication and is required for any claim that code constitutes scientific infrastructure.

**Milestone 2 — Conditioning diagnostic**  
For Hejhal system at fixed $Y_0$, plot $\log\kappa(V(r))$ vs $\log M$ for $M=100,200,400,800$. Exponential growth flags interval failure, polynomial growth acceptable. Current diagnostic without preconditioning shows $M^{76}$, steep growth verdict, increase $Y_0$ or reduce $M$ or add diagonal preconditioner.

**Milestone 3 — Width versus $N$ diagnostic**  
Free with existing prototypes, already described under scaling conjectures.

**Milestone 4 — Cross validation against Then spectrum**  
Once dual certification at level 1 succeeds, compare rigorous enclosure with Then high precision non rigorous value. Agreement to many digits is strong though not logical confirmation.

---

## Heuristic Motivation for $\mu(N)$, Appendix, explicitly not claimed

This section is explicitly labelled heuristic and does not form part of any proof.

As $N$ grows the discrete operator on reference cell tiling is essentially $N+1$ copies of level 1 operator plus gluing permutation. First new eigenvalue controlled by spectral gap of Schreier graph of $\mathrm{PSL}_2(\mathbb{F}_N)/B$. These graphs are expanders, Lubotzky Phillips Sarnak, Bourgain Gamburd, diameter $O(\log N)$, Cheeger Buser comparison yields gap $\gtrsim c/\log N$. This supplies plausible heuristic that $\mu(N)\ge c/\log N$ rather than power law decay. It is not a proof, nor used in any certification. It merely motivates why logarithmic form is natural to test.

---

## Risk Register

| Risk | Severity | Mitigation |
|------|----------|------------|
| Defect constant $C(K,Y,r)$ blows up as $r\to0$ | High | Large $Y$, explicit Arb tracking, increase $Y_0$ |
| Hejhal conditioning exponential in $M$ | High | Diagnostic plot, increase $Y_0$, diagonal preconditioner, reduce $M$ |
| FEM upper bound never reaches 45 | Medium | Accept coarser window, rely more on Hejhal plus counting |
| Interval width grows with $N$ | Medium | Conjecture A plus free diagnostic, width vs $N$ plot |
| Route A counting error too large | Medium | Fall back to Route B |
| Spectral clustering, $\lambda_2-\lambda_1\approx10^{-8}$ makes defect, overlap, counting harder | Medium | Add to risk register, use Lehmann Goerisch with block, require higher precision, note Then spectrum shows gaps $>1$ at low end for Picard, but must be checked for congruence levels |
| Hidden Hecke growth dependence | Low | Explicit Assumption H in every paper, dependence only via $C_H(\varepsilon)$ and $\theta$ |

Spectral clustering belongs in register even if believed not to occur for low Picard eigenvalues, because defect bounds and interval overlap and counting all become harder when eigenvalues cluster.

---

## Immediate High Impact Experiments, next 30 days

| Priority | Experiment | Unlocks | Effort |
|----------|------------|---------|--------|
| 1 | Lemma K plus `constants_DK.md` Arb tables | All later defect bounds | 3 to 5 d |
| 2 | Reproducibility suite, independent meshes | Scientific credibility, explicit milestone | 1 w |
| 3 | Route A counting prototype on $K_Y$ product mesh, addresses remaining logical gap | Closes logical hole without Hejhal | 1 to 2 w |
| 4 | Single cusp interval Hejhal prototype with diagonal preconditioner | Existence half of dual loop | 2 w |
| 5 | Width vs $N$ and kappa vs $M$ plots with preconditioning | Safety of scaling | 2 d |

---

## Final Remark on Language

* Use engineering target for FEM upper bound near 45 until it exists.
* Use addresses the remaining logical gap until counting theorem exists, not closes logical hole.
* Use Then $N=5$ becomes first congruence level where present trace formula exclusion criterion no longer closes the spectral gap, making dual certification pipeline the proposed route to rigorous certification, not only your pipeline can.
* Dual certification leads to high confidence subject to assumptions, counting, conditioning, implementation, proof still comes from each method satisfying its hypotheses.

---

*This document is intended for dialectic iteration. Every box that still contains unproved analytic statement is deliberately labelled as such. Goal is not to claim theorem prematurely, but to make remaining theorems precise, independent, and attackable. The larger picture is the dual loop, not any single rung.*
