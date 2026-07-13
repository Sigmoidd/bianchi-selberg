# Dual Certification Roadmap for \(\mathrm{Spec}_\Delta(M_{\mathfrak{p}})\)

**Version 3 — Architecture Document**  
**Date:** 2026-07-11  
**Status:** Research program roadmap (not a completed theorem)

---

## Contribution Statement

The contribution is a **certification architecture** that combines independent variational (FEM) and automorphic (interval Hejhal) methods into a single rigorous workflow for arithmetic hyperbolic orbifolds arising from Bianchi groups and their congruence subgroups.

A single method cannot close the bottom of \(\mathrm{Spec}_\Delta(M_{\mathfrak{p}})\). Certified exclusion of \((0,1)\) plus existence of an eigenvalue near \(45\) still leaves open holes. The architecture supplies the missing piece — certified counting — while making every logical dependency and every unproved analytic hypothesis explicit.

```
Arithmetic group
    ↓ exact congruence tiling over 𝒪_K/𝔭
    ↓ certified FEM (Lax–Phillips + Crouzeix–Raviart + Rump + Arb)
    ├── certified exclusion on (0,1)
    ├── certified upper bounds (engineering target)
    └── Route A / Route B counting
    ↓ interval Hejhal with two-cusp coupling
    ↓ certified existence in [λ−ε, λ+ε]
    ↓ mutual verification (overlap)
    ↓ rigorous first eigenvalue
```

---

## Separation of Concerns

**Known Mathematics**  
Spectral decomposition for cofinite Kleinian groups (Friedman, Thm 3.8.1), Shimizu horoball lemma, Lax–Phillips cusp reduction, Payne–Weinberger Poincaré constant, Crouzeix–Raviart interpolation estimates, Rump’s positive-definiteness verification, Selberg trace formula as an identity.

**Engineering (already realized in the repository)**  
Turning known inequalities into Arb enclosures with explicit constants \(h_T,\gamma,\tau,S_M,S_Q,S_\Sigma,V_\Sigma\); managing \(n\sim 10^4\)–\(3\cdot 10^4\) degrees of freedom; residue-ring arithmetic \(\mathbb{Z}[i]/\mathfrak{p}\cong\mathbb{F}_N\) and \(\mathbb{Z}[i]/(3)\cong\mathbb{F}_9\); gluing of \(N\mathfrak{p}+1\) isometric copies of a reference cell; Rump certificates with shifts of order \(10^{-8}\).

**Unknown Theorems (explicit research objectives)**  
1. Explicit defect bounds for Bianchi groups (single-cusp and two-cusp).  
2. Uniform \(K_{ir}\) tails with Hecke-controlled coefficients.  
3. Conditioning of the coupled two-cusp Hejhal system.  
4. Certified eigenvalue counting (two independent routes).  
5. Scaling of interval widths under reference-cell tiling.  
6. Uniform spectral-gap lower bound \(\mu(N)\).

---

## Logical Architecture & Dependency Graph

The pipeline is **not** a linear chain. FEM is the foundational engine; Hejhal supplies independent existence; counting is a refinement that may reuse FEM data.

```
                    ┌─────────────────────────────┐
                    │   FEM (Crouzeix–Raviart +   │
                    │   Rump + Arb reference cell)│
                    └─────────────┬───────────────┘
                                  │
          ┌───────────────────────┼───────────────────────┐
          │                       │                       │
          ▼                       ▼                       ▼
   Exclusion of (0,1)      Upper bounds on λ₁      Route A / Route B
   (already certified      (engineering target:    Certified counting
    for N𝔭 = 1,5,9,13)     reach ~45)              N(λ)
          │                       │                       │
          └───────────────────────┼───────────────────────┘
                                  │
                                  ▼
                    ┌─────────────────────────────┐
                    │  Interval Hejhal            │
                    │  (defect bound + Krawczyk)  │
                    │  → existence [λ−ε, λ+ε]     │
                    └─────────────┬───────────────┘
                                  │
                                  ▼
                    ┌─────────────────────────────┐
                    │  Dual verification          │
                    │  FEM interval ∩ Hejhal      │
                    │  interval ≠ ∅               │
                    │  + counting ⇒ first         │
                    │  eigenvalue rigorous        │
                    └─────────────────────────────┘
```

This graph makes two facts transparent:

- Counting is **not** an independent oracle; it typically reuses FEM meshes, flux reconstructions, or approximate eigenfunctions.
- Hejhal is the only component that is currently completely absent from the repository.

---

## Dual Certification Loop (the persuasive core)

```
FEM certifies          λ₁ ∈ [a,b]
                       (lower bound from Theorem G1 / CR;
                        upper bound from conforming Rayleigh quotient
                        — engineering target)

Hejhal independently   λ₁ ∈ [c,d]
                       (defect bound of Child/BSV type
                        with explicit C(K,Y,r))

If [a,b] ∩ [c,d] ≠ ∅ and diam([a,b]∪[c,d]) < tolerance
then the two intervals mutually verify each other.
```

Two mathematically distinct methods — one variational, one automorphic — certify the same eigenvalue. Overlap is far stronger evidence than either method alone. Non-overlap immediately localizes the defect (mesh size vs. Fourier tail vs. conditioning).

**Important caveat.**  
The statement “FEM upper bound \(44 < \lambda_1 < 46\)” is currently an **engineering target**, not an existing capability. The present CR machinery was engineered and certified for exclusion of \((0,1)\). Reaching a certified upper bound near \(45\) requires \(hp\)-refinement or a substantially larger mesh and is therefore listed as Rung 2 work, not as already-available data.

---

## Two Independent Routes for Certified Counting

These methods solve related but non-identical problems and rest on different hypotheses. They must **not** be treated as interchangeable.

### Route A — Dirichlet–Neumann Bracketing
- Natural for truncated domains that respect the cusp-collar product structure.
- Produces two-sided bounds
  \[
  \lambda_k^{N,K_h} \le \lambda_k \le \lambda_k^{D,K_h}
  \]
  by imposing Neumann or Dirichlet conditions on the artificial boundary of a compact core \(K\).
- Requires a mesh that is product-like in the collar and guaranteed FEM eigenvalue bounds on that mesh (already available from the CR engine).
- Advantage: purely variational; reuses existing reference-cell infrastructure.
- Risk: artificial boundary conditions introduce a controllable but non-zero error that must be enclosed.

### Route B — Lehmann–Goerisch Refinement
- Two-stage method: given a lower bound \(\rho\) on \(\lambda_{m+1}\) (supplied by the already-certified exclusion \(\rho=1\)), one obtains rigorous lower bounds on the preceding \(n\) eigenvalues \(\lambda_{m-n+1},\dots,\lambda_m\).
- Requires approximate eigenfunctions and flux reconstructions (hypercircle or residual equilibration).
- Advantage: can bootstrap from a single exclusion all the way to several lower eigenvalues without a full spectrum.
- Risk: quality of the approximate eigenfunctions and the flux reconstruction must themselves be certified.

**Architectural rule.**  
Both routes are first-class citizens of the roadmap. The architecture survives even if one route fails for technical reasons. In practice Route A will be attempted first because it reuses the existing FEM codebase more directly; Route B is the natural refinement once approximate eigenfunctions become available from Hejhal.

---

## Explicit Analytic Assumptions

Every certification path that uses Fourier tails or defect bounds must declare its coefficient growth hypothesis.

**Assumption (Hecke-type growth).**  
There exist constants \(C_\varepsilon>0\) (depending only on the group and the character) such that the Fourier coefficients of any Maass form of eigenvalue \(\lambda=1+r^2\) satisfy
\[
|a_\beta| \le C_\varepsilon\, N(\beta)^{1/2+\varepsilon}
\]
for every \(\varepsilon>0\).

This is the standard bound available from the theory of automorphic forms on \(\mathrm{GL}_2\) over imaginary quadratic fields (or from the Ramanujan conjecture at finite places, which is known in many cases). It is **not** free; it must be cited and the dependence of all subsequent constants on \(C_\varepsilon\) must be tracked.

**Lemma K (foundational, conditional).**  
Under the Hecke-type growth assumption, for \(|\beta|\ge M\) and \(y\ge Y_0\),
\[
\sum_{|\beta|\ge M} |a_\beta|^2 |K_{ir}(2\pi|\beta|Y_0)|^2 \le \varepsilon
\]
with an explicit Arb-enclosable \(\varepsilon = \varepsilon(r,Y_0,M,C_\varepsilon)\).

This lemma is required both for the Hejhal defect bound and for the cusp Fourier tail control inside Lax–Phillips analysis. It is the first item on the challenge ladder precisely because every later step inherits its constants.

---

## Scaling Conjectures

**Conjecture A (primary).**  
Let \(W(N)\) be the supremum of the interval radii appearing in the glued stiffness and mass matrices for the reference-cell tiling at level \(N=N\mathfrak{p}\). Then
\[
W(N) = O(1)
\]
as \(N\to\infty\), provided the local mesh size satisfies \(h_T\le c\, Y/\sqrt{N}\) and the per-row radius absorption of Paper I is used.

**Evidence.**  
The geometry consists of \(N+1\) isometric copies of a fixed reference cell. The stiffness matrix is block-diagonal plus a gluing permutation of bounded degree. Local interval radii therefore depend only on the degree of a vertex and the local mesh geometry, both of which are independent of \(N\).

**Fallback Conjecture.**  
If the primary conjecture fails (e.g., because of subtle accumulation of rounding or residual coupling), then
\[
W(N) = O(\log N).
\]

**Engineering check.**  
Run the existing gluing prototype for \(N=5,9,13,25,37\) and plot the maximum entry-wise radius of the assembled interval matrix against \(N\). Flatness supports Conjecture A; logarithmic growth supports the fallback. This experiment is essentially free with current code.

---

## Challenge Ladder (refined)

**Rung 0 — K-Bessel enclosure (foundational, conditional on Hecke growth)**  
Lemma K with explicit constants, Arb implementation, validation against Luke-type inequalities. No \(N\)-dependence. **Must be finished first.**

**Rung 1 — Single-cusp defect bound**  
Theorem D(\(K\)): if a putative form has automorphy defect \(\le\delta\) on side-pairings and cusp Fourier tail \(\le\tau\) in \(L^2(K_Y)\), then there exists a true cusp form with
\[
|\tilde\lambda - \lambda_{\mathrm{true}}| \le C(K,Y,r)(\delta+\tau)
\]
and \(C\) given by an Arb-enclosable expression. This is the analytic core and will be a full paper.

**Rung 2 — Level-1 dual certification (engineering target)**  
- FEM: certified exclusion of \((0,1)\) (already done) + certified upper bound in a window around \(45\) (requires mesh refinement; currently an engineering target).  
- Interval Hejhal at the published Then value \(r_1\approx 6.62212\) with \(M=400\), \(Y_0=0.8\), 128-bit Arb.  
- Route A or Route B counting establishing \(N(44.8)=0\).  
- Overlap of the two intervals yields mutual verification of the first eigenvalue.

**Rung 3 — Two-cusp coupling analysis**  
Well-posedness, consistency of the \(\sigma_0\)-gluing, and a condition-number bound for the block matrix
\[
\begin{pmatrix} V_\infty & -S \\ -S^* & V_0 \end{pmatrix}.
\]
Implement interval Krawczyk for \(N=5\).

**Rung 4 — \(N=5\) dual certification**  
\(\mathfrak{p}=(2+i)\), index 6. FEM lower bound already certified; add Hejhal existence and counting. Overlap gives mutual verification at the first non-trivial congruence level.

**Rung 5 — Scaling theorem**  
Prove (or rigorously verify on a sequence of levels) Conjecture A or its fallback. Numerical evidence from the free experiment above is required before claiming the jump from \(N=5\) to \(N=13,25,\dots\) is safe.

**Rung 6 — Field independence**  
Port the entire pipeline to \(\mathbb{Z}[\omega]\) and to inert primes (\(N=9\)). Tests that the reference-cell principle is independent of the ring of integers.

**Explicitly left as conjectures (not claims)**  
- Uniform \(\mu(N)\ge c/\log N\).  
- Explicit Weyl remainder with constants for cusped Bianchi orbifolds.  
- Simplicity of \(\lambda_1\).

---

## Required Benchmarks (scientific infrastructure)

1. **Reproducibility benchmark (mandatory)**  
   Take the published Picard eigenvalue \(r_1\approx 6.62212\).  
   Re-run the entire dual-certification pipeline with  
   - completely independent mesh generators,  
   - different quadrature rules,  
   - different node numberings,  
   - different compilers and floating-point libraries,  
   - different Arb precisions (64 / 128 / 256).  
   The final certified interval must be identical (up to the declared enclosure radius).  
   This is the software analogue of independent laboratory replication and is required for any claim that the code constitutes scientific infrastructure.

2. **Conditioning diagnostic**  
   For the Hejhal system at fixed \(Y_0\), plot \(\log\kappa(V(r))\) versus \(\log M\) for \(M=100,200,400,800\).  
   Exponential growth immediately flags that interval arithmetic will fail; polynomial growth is acceptable.

3. **Width-versus-\(N\) diagnostic**  
   Already described under Scaling Conjectures. Free with existing prototypes.

4. **Cross-validation against Then’s spectrum**  
   Once dual certification at level 1 succeeds, compare the rigorous enclosure with Then’s high-precision non-rigorous value. Agreement to many digits is strong (though not logical) confirmation.

---

## Heuristic Motivation for \(\mu(N)\) (Appendix)

*(This section is explicitly labelled heuristic and does not form part of any proof.)*

As \(N\) grows the discrete operator on the reference-cell tiling is essentially \(N+1\) copies of the level-1 operator plus a gluing permutation. The first new eigenvalue is therefore controlled by the spectral gap of the Schreier graph of \(\mathrm{PSL}_2(\mathbb{F}_N)/B\). These graphs are expanders (Lubotzky–Phillips–Sarnak, Bourgain–Gamburd). Expanders of degree \(d\) have diameter \(O(\log N)\), which by the usual Cheeger–Buser comparison yields a spectral gap of size \(\gtrsim c/\log N\).

This supplies a plausible heuristic that
\[
\mu(N) \ge \frac{c}{\log N}
\]
rather than a power-law decay. It is **not** a proof, nor is it used in any certification. It merely motivates why the logarithmic form of the uniform-gap conjecture is the natural one to test.

---

## Immediate High-Impact Experiments (next 30 days)

| Priority | Experiment                              | Unlocks                          | Effort |
|----------|-----------------------------------------|----------------------------------|--------|
| 1        | Lemma K (K-Bessel + Hecke growth)       | All later defect bounds          | 3–5 d  |
| 2        | Level-1 FEM upper-bound target (~45)    | Dual loop at level 1             | 1–2 w  |
| 3        | Single-cusp interval Hejhal prototype   | Existence half of dual loop      | 2 w    |
| 4        | Width-vs-\(N\) and \(\kappa\)-vs-\(M\) plots | Safety of scaling               | 2 d    |
| 5        | Reproducibility suite (independent meshes) | Scientific credibility          | 1 w    |

---

## Risk Register

| Risk                                      | Severity | Mitigation                                      |
|-------------------------------------------|----------|-------------------------------------------------|
| Defect constant \(C(K,Y,r)\) blows up as \(r\to0\) | High    | Large \(Y\), explicit Arb tracking              |
| Hejhal conditioning exponential in \(M\)  | High     | Diagnostic plot; increase \(Y_0\); reduce \(M\) |
| FEM upper bound never reaches 45          | Medium   | Accept coarser window; rely more on Hejhal      |
| Interval width grows with \(N\)           | Medium   | Conjecture A + free diagnostic                  |
| Route A counting error too large          | Medium   | Fall back to Route B                            |
| Hidden Hecke-growth dependence            | Low      | Explicit Assumption statement in every paper    |

---

## Final Remark on Applications

Once the dual-certification pipeline is operational, applications to CMB topology, spectral methods for machine learning, audio processing, and Grassmannian optimization become legitimate **downstream** beneficiaries. They are not the primary justification for the architecture. The primary justification is the production of rigorously certified Laplace spectra on arithmetic hyperbolic 3-orbifolds — a foundational capability that currently does not exist.

---

*This document is intended for dialectic iteration. Every box that still contains an unproved analytic statement is deliberately labelled as such. The goal is not to claim a theorem prematurely, but to make the remaining theorems precise, independent, and attackable.*
