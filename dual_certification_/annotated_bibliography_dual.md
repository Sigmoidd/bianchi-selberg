# Annotated Bibliography v2: Dual Certification for $\mathrm{Spec}_\Delta(M_{\mathfrak{p}})$

**Thesis:** A single method cannot close the bottom of $\mathrm{Spec}_\Delta(M_{\mathfrak{p}})$. Exclusion on $(0,1)$ by verified FEM plus existence above $1$ by interval Hejhal still leaves holes at $3,7,12,18,\dots$. The missing theorem is certified counting. The contribution is the bridge:

```
Arithmetic group
↓ exact congruence tiling over O_K/p
↓ certified FEM (Lax-Phillips + Crouzeix-Raviart + Rump + Arb)
↓ certified exclusion on (0,1)
↓ interval Hejhal with two-cusp coupling
↓ certified existence in [λ-ε,λ+ε]
↓ certified counting (Dirichlet-Neumann bracketing / Lehmann-Goerisch)
↓ rigorous first eigenvalue
```

This version is organized as Known Mathematics / Engineering / Unknown Theorems, as a referee expects.

---

## Separation of Concerns

**Known Mathematics:** Spectral decomposition for cofinite Kleinian groups (Friedman Thm 3.8.1), Shimizu horoball lemma, Lax-Phillips cusp reduction, Payne-Weinberger Poincaré constant, Crouzeix-Raviart interpolation, Rump PSD verification, Selberg trace formula as identity.

**Engineering:** Turning known inequalities into Arb enclosures with explicit $h_T$, $\gamma$, $\tau$, $S_M$, $S_Q$, $S_\Sigma$, $V_\Sigma$, managing $n\sim 10^4$ to $3\cdot10^4$ dofs, keeping interval widths from exploding, implementing $\mathbb{Z}[i]/\mathfrak{p}\cong\mathbb{F}_{N}$ and $\mathbb{Z}[i]/(3)\cong\mathbb{F}_9$ residue arithmetic, and gluing $N\mathfrak{p}+1$ copies.

**Unknown Theorems:** Explicit defect bounds for Bianchi groups, uniform $K_{ir}$ tails, conditioning of coupled two-cusp Hejhal system, certified eigenvalue counting without asymptotics, scaling of interval width with $N\mathfrak{p}$, uniform $\mu(N)$ lower bound.

---

## 1. The biggest logical gap: exclusion + existence ≠ first eigenvalue

**Reviewer attack:** Certifying no eigenvalues in $(0,1)$ and an eigenvalue in $[44.8,44.9]$ does not exclude $3,7,12,18,25,39$.

**What is needed:** Certified counting $N(\lambda)=\#\{j:\lambda_j\le\lambda\}$ with explicit two-sided bounds, not Weyl asymptotics.

**Literature:**

* **Carstensen et al, Direct guaranteed lower eigenvalue bounds with optimal a priori convergence.** An extra-stabilised Morley FEM directly computes guaranteed lower eigenvalue bounds with optimal a priori convergence rates for the bi-Laplace Dirichlet eigenvalues, the smallness assumption $\min\{\lambda_h,\lambda\}h_{\max}^4\le184.95$ makes $\lambda_h\le\lambda$ a lower bound[^1]. Same philosophy gives fully computable two-sided bounds on Laplace eigenvalues via Crouzeix-Raviart plus postprocessing, efficiency proven for graded meshes[^2].

* **Liu — Guaranteed lower bounds via hypercircle.** For Steklov and Laplacian, Liu's approach gives guaranteed lower bounds requiring a priori error estimation for nonhomogeneous Neumann problems, solved by constructing the hypercircle. This is the Dirichlet-Neumann bracketing route: $\lambda_k^{N,K_h}\le\lambda_k\le\lambda_k^{D,K_h}$ on a mesh that respects the cusp collar product structure.

* **Lehmann-Goerisch theorem, Two-Stage approach.** The following theorem gives an application to obtaining lower bounds for $n$ eigenvalues $\lambda_{m-n+1},\dots,\lambda_m$ by using a lower bound $\rho$ of $\lambda_{m+1}$[^3]. This is how you bootstrap from your exclusion $(0,1)$ as $\rho$ to bound $\lambda_2,\lambda_3,\dots$ without computing them.

* **Weyl remainder with explicit constants.** Recent work gives quantitative estimates on the remainder of Weyl's law with explicit constants without using Neumann eigenvalues[^4]. The eigenvalue counting function $N(\lambda)$ satisfies Weyl's asymptotic $N(\lambda)\sim C\,\mathrm{vol}\,\lambda^{3/2}$[^5][^6], but explicit remainder $R(\lambda)=N(\lambda)-C\,\mathrm{vol}\,\lambda^{3/2}$ with interval enclosure is what you would need if you tried to use Weyl for counting. Do not use asymptotics as exact counting; use it to choose search windows.

**Unknown to turn into theorem:** Certified $N(44.8)=0$ for $M_1=\mathrm{PSL}_2(\mathbb{Z}[i])\backslash\mathbb{H}^3$, i.e., a verified upper bound on counting function that matches your lower bound from exclusion. This closes the hole.

---

## 2. FEM and Hejhal certify different objects — the overlap problem

**Reviewer attack:** FEM naturally certifies low spectrum $\lambda\le15$, Hejhal naturally succeeds around $\lambda\approx45$ where $K_{ir}$ oscillates enough to resolve. Giant hole in between.

**Literature:**

* **Booker, Strömbergsson, Venkatesh 2006.** Proves first ten eigenvalues on $\mathrm{PSL}_2(\mathbb{Z})\backslash\mathbb{H}$ correct to at least 100 digits, but also proves certification task is polynomial time given sufficient accuracy input. That sufficiency is the bridge: if you feed Hejhal $100$ digits, BSV certifies; if you feed $10$ digits, it does not. For Bianchi, low eigenvalues have large $y$ decay $K_{ir}(y)\sim e^{-y}$, so Hejhal system becomes ill conditioned as $r\to0$.

* **Then — Large sets of consecutive Maass forms.** Uses Hejhal's identity to compute eigenfunctions corresponding to given eigenvalues and linearisation. Shows consecutive search is possible by sweeping $r$ and tracking sign changes of determinant, which is a counting mechanism if made interval rigorous.

**Unknown to turn into theorem:** Either extend FEM upward via $hp$-refinement to reach $r\sim6.6$ with certified upper bound $<45$, or extend Hejhal downward to $r\sim\sqrt{14}$ ($\lambda=15$) by increasing $M(Y_0)$ and $Y$ and controlling condition number. The dual certification loop the reviewer suggests is: FEM certifies $\lambda_1\in[44,46]$, Hejhal independently certifies $\lambda_1\in[44.8,44.9]$, intervals overlap $\Rightarrow$ mutual verification, far stronger than either alone.

---

## 3. Defect bounds — the hardest analysis

**Reviewer attack:** Child proves small defect $\Rightarrow$ nearby eigenvalue for $\mathrm{SL}_2(\mathbb{R})$. Bianchi needs new Sobolev, cusp, and boundary estimates.

**Literature:**

* Child's method: If $\tilde f$ is putative and $\|(\Delta-\lambda)\tilde f\|_2$ small, spectral resolution shows closeness, and it suffices to obtain strong bounds along boundary of fundamental domain. For Bianchi, boundary is a 2-torus $T^2=\mathbb{C}/\mathcal{O}_K$ times $S^1$ from rotation.

* Need explicit: trace inequality $\|w\|_{L^2(F)}^2\le C_F(\|w\|_{L^2(T)}^2+h_T\|w\|\|\nabla w\|)$, your Lemma T, plus interior elliptic regularity in $\mathbb{H}^3$ with weight $y^{-3}dx\,dy$.

**Unknown to turn into theorem:** Theorem D($K$): Let $\tilde f$ satisfy automorphy defect $\le\delta$ on side pairings and cusp Fourier tail $\le\tau$ in $L^2(K_Y)$. Then exists true cusp form with $|\tilde\lambda-\lambda_{\mathrm{true}}|\le C(K,Y,r)(\delta+\tau)$. Explicit $C$ with Arb enclosure. This is one full paper.

---

## 4. Two-cusp coupling

**Reviewer attack:** Two expansions coupled by $\sigma_0=\left(\begin{smallmatrix}0&-1\\1&0\end{smallmatrix}\right)$ are not duplication.

**Structure:** At $\infty$, $\Lambda_\infty=\mathcal{O}_K$, $|T_\infty|=\mathrm{covol}$. At $0$, $\Lambda_0=\mathfrak{p}$, $|T_0|=N\mathfrak{p}\cdot|T_\infty|$, dual lattice shortest vector $1/\sqrt{N\mathfrak{p}}$. Mode bound $4\pi^2Y^2/N\mathfrak{p}>1$ ensures zero mode dominates in collar, your C4 check.

**Unknowns:**
* Consistency: $f_0 = f_\infty|_{\sigma_0}$ as distributions on $y=Y$.
* Uniqueness: If both expansions satisfy same eigenvalue, they glue to automorphic function iff Fourier coefficients satisfy linear relation $a^{(0)} = S(r) a^{(\infty)}$ where $S(r)$ is scattering-like but at finite $Y$, not $Y\to\infty$.
* Conditioning: Coupled matrix $\begin{pmatrix}V_\infty & -S\\ -S^* & V_0\end{pmatrix}$ may have condition number $\kappa\sim\exp(c\,M)$.
* Interval stability: Need Krawczyk test for this block matrix, not just Rump PSD.

**Literature:** Strömberg's reduction algorithm is the only existing code that handles $\Gamma_0(N)\subset\mathrm{SL}_2(\mathbb{Z})$ with two cusps; his Hilbert modular generalization shows the same coupling appears for $\mathbb{Q}(\sqrt{d})$ real. No interval version exists.

---

## 5. Uniform $K$-Bessel bounds — foundational

**Reviewer attack:** Without explicit interval bounds certification collapses.

**Literature:** Bounds for modified Bessel function $K_\nu(x)$ have attracted attention, Luke proved double inequality for $x>0$[^7], Gaunt and others refined, and recent work gives simple inequalities for integrals involving $I_\nu,K_\nu$ and monotonicity results for $K_\nu$ and new lower bound involving gamma functions[^8]. For $K_{ir}(y)$ with $r\in\mathbb{R}$, need uniform $y\ge Y_0$ tail: $|K_{ir}(y)|\le \sqrt{\pi/(2y)}e^{-y}\cdot C(r)$ with $C(r)$ explicit.

**Unknown to turn into theorem:** Lemma K($r,Y_0,M$): For $|\beta|\ge M$, $\sum_{|\beta|\ge M} |a_\beta|^2 |K_{ir}(2\pi|\beta|Y_0)|^2 \le \varepsilon$ with $\varepsilon$ Arb-enclosed given a priori polynomial bound on $a_\beta$ from Hecke bound $|a_\beta|\le C_\varepsilon N(\beta)^{1/2+\varepsilon}$. This is needed for both Hejhal defect and cusp Fourier tail in Lax-Phillips.

---

## 6. Weyl law — validation, not certification

**Distinction:** Weyl gives $N(\lambda)= \frac{\mathrm{vol}}{6\pi^2}\lambda^{3/2}+a_2\lambda\log\lambda + a_3\lambda +a_4+o(\lambda)$. Your `verify_matthies.py` reconstructs $a_2$ and full $a_3$ bracket from engine constants, which validates elliptic sector completeness. It does **not** count individual eigenvalues. Reviewer will object if you slip from asymptotics to exact count.

**What you need instead:** Explicit remainder with constants, as in recent quantitative Weyl law with explicit constants without Neumann eigenvalues[^4], adapted to Bianchi orbifolds with cusps. This is hard. Better to use Dirichlet-Neumann bracketing on truncated core $K$ to get $N_{D}( \lambda)\le N(\lambda)\le N_{N}(\lambda)$ with both sides computable via FEM with guaranteed bounds[^2].

---

## 7. Interval linear algebra — conditioning and width explosion

**Reviewer attack:** Naively $[A]x=[b]$ intervals explode.

**Literature:**
* Krawczyk's test established by Rump, Theorem 3.1 for interval $X$, if interior condition holds then unique solution exists[^9][^10].
* Rump presented computationally simple and fast sufficient criterion implying positive definiteness of symmetric or Hermitian interval matrix[^11][^12].
* Rump algorithm works up to condition number about $10^{15}$[^13], beyond that fails.
* Determinant and inverse via naive Laplace expansion is factorial in dimension, motivating specialised interval linear algebra[^14].

**What you have:** You already use Rump for PSD certificate, $c_e>d_e$ and Rump shift $\approx2\cdot10^{-8}$ in Paper I Table 1. For Hejhal you need Krawczyk for non symmetric $V(r)$, not just PSD.

**Unknown:** Prove $\kappa(V(r)) = O(M^\alpha)$ for Bianchi Hejhal with $Y_0\asymp1$, not $\exp(M)$. If $\kappa$ grows exponentially, interval width explodes and certification fails regardless of precision. Benchmark on $M=100,200,400$ and fit log $\kappa$ vs $\log M$.

---

## 8. Scaling — do not jump $1\to5\to13$ without a theorem

**Reviewer suggestion:** Insert intermediate theorem on growth of interval width.

**Proposal:** Theorem S: Let $W(N)=\sup$ interval radius for gluing data at $N\mathfrak{p}=N$. Then $W(N)=O(\log N)$ or $W(N)=O(N^\alpha)$ with explicit $\alpha<1$, provided $Y$ fixed and $h_T\le c\,Y/\sqrt{N}$.

**Why plausible:** Reference cell principle says geometry is $N+1$ isometric copies, so stiffness matrix $Q$ is block diagonal plus gluing permutation. Per row radius absorption from Paper I makes radius depend on degree of vertex, not on $N$. So width should be $O(1)$ in $N$, with only linear algebra cost $O(N^3)$ growing.

**Engineering check:** Run `congruence_prototype.py gluing` for $N=5,9,13,25,37$ and plot $\|Q_r\|$ vs $N$. If flat, you have evidence for theorem.

---

## 9. Uniform $\mu(N)$ conjecture — motivate logarithmic

**Current statement:** $\mu(N)\ge c/\log N$ appears out of nowhere.

**Heuristic to add:** $\mu$ is Rayleigh quotient gap $\inf_{v\perp\ker Q} v^T Q v / v^T M v$ minus $s(2+s)$ continuum threshold. As $N$ grows, $Q$ is $N+1$ copies, so first new eigenvalue should scale like first non zero Laplacian eigenvalue of base graph where vertices are copies and edges are gluing permutations. That graph is a Schreier graph of $\mathrm{PSL}_2(\mathbb{F}_{N})/B$, which is an expander with Cheeger constant $\ge c$, so gap $\ge c'/\log N$ by diameter bound $O(\log N)$. This motivates logarithmic, not $1/N$.

**Literature to cite:** Lubotzky, Phillips, Sarnak Ramanujan graphs, Bourgain Gamburd expansion for $\mathrm{SL}_2(\mathbb{F}_p)$. Do not claim proof, present as heuristic that makes conjecture plausible.

---

## 10. What you are underselling — the bridge as architectural contribution

Papers about each box exist: Lax-Phillips, Carstensen guaranteed lower bounds, Rump verification, Hejhal numerics, Then Bianchi numerics, BSV certification. Few papers connect all boxes into one certification pipeline for Bianchi orbifolds.

**Dual certification loop, the most persuasive form:**

```
FEM certifies λ₁ ∈ [a,b]   (lower bound from Theorem G1, upper bound from conforming Rayleigh quotient)
Hejhal independently certifies λ₁ ∈ [c,d]   (defect bound from Child/BSV)
[a,b] ∩ [c,d] ≠ ∅  and diameter < tolerance  ⇒  mutual verification
```

If intervals overlap, two mathematically distinct methods, one variational, one automorphic, certify same eigenvalue. If they do not overlap, gap is a concrete research problem, not vague uncertainty.

**Immediate next month experiment:** Take Paper I level 1, $r_1=6.62212$. Run FEM with refined mesh to get certified upper bound $\lambda_h^{\mathrm{conf}}$ and lower bound $\lambda_h^{\mathrm{CR}}$. Run Hejhal with $M=400$, $Y_0=0.8$, Arb $128$ bits to get interval $[r_1-10^{-6},r_1+10^{-6}]$. Do they overlap? If yes, you have dual certification at level 1. If not, measure distance and identify whether defect is FEM $h_T$ or Hejhal tail.

---

## Revised Challenge Ladder for Dual Certification

**Rung 0 — K-Bessel enclosure (foundational)**
Lemma K with explicit constants, Arb implementation, test against Luke double inequality. No $N$ dependence.

**Rung 1 — Single cusp defect bound**
Theorem D($K$) for $K=\mathbb{Q}(i)$, one cusp, explicit $C(K,Y,r)$. Paper length. This is the analysis core.

**Rung 2 — Level 1 dual certification**
FEM exclusion $(0,1)$ done, plus FEM upper bound $44<\lambda_1<46$, plus Hejhal existence $[44.8,44.9]$, plus counting $N(44.8)=0$ via Dirichlet-Neumann bracketing. Result: first eigenvalue rigorously $44.8\ldots$.

**Rung 3 — Two-cusp coupling analysis**
Prove well posedness and condition number bound for coupled system $V_\infty\oplus V_0$ with $\sigma_0$ consistency. Implement interval Krawczyk for $N=5$.

**Rung 4 — $N=5$ dual certification**
$\mathfrak{p}=(2+i)$, index $6$, volume $1.83$, $B\approx1.9>1$ so trace formula cannot close. Show FEM lower bound $\mu\approx4.4$ still holds, plus Hejhal existence for first eigenvalue above $1$. Overlap gives mutual verification.

**Rung 5 — Scaling theorem**
Prove $W(N)=O(\log N)$ and $\kappa(V)=O(N^\alpha)$ for reference cell tiling, with numerical evidence from $N=5,9,13,25$. This justifies jumping from $5$ to $13$ and beyond without surprise blow up.

**Rung 6 — Eisenstein and $N=9,13$**
Port Rung 2-4 to $\mathbb{Z}[\omega]$, $|T_\infty|=\sqrt3/6$, rotation order $3$, and to inert $N=9$ where residue field is $\mathbb{F}_9$. This tests field independence of reference cell principle.

**Unknowns left as conjectures, not claims:** Uniform $\mu(N)\ge c/\log N$, explicit Weyl remainder with constants for Bianchi orbifolds, simplicity of $\lambda_1$.

---

## Minimal bibliography to include

- BSV 2006 effective computation, 100 digits, polynomial time certification
- Child 2022 arbitrary level and character, bound on difference authentic vs purported, level 5 quadratic character first example,, first with non trivial character
- Then large consecutive sets, Hejhal identity, linearisation, Maass form specified by eigenvalue and finite coefficients
- Then arithmetic quantum chaos, no small eigenvalues $0<\lambda<1$ for Picard, four symmetry classes D G C H
- Strömberg $\Gamma_0(N)$ computation and reduction algorithm as key ingredient for Hejhal,
- Carstensen guaranteed lower bounds, Morley FEM gives $\lambda_h\le\lambda$ under smallness[^1], fully computable two-sided bounds via CR plus postprocessing[^2], Lehmann-Goerisch for $n$ eigenvalues using $\rho$ of $\lambda_{m+1}$[^3]
- Weyl remainder explicit constants without Neumann eigenvalues[^4], counting function $N(\lambda)$ asymptotic[^5][^6]
- K-Bessel bounds attract attention, Luke double inequality[^7], inequalities for integrals of $I_\nu,K_\nu$ and monotonicity[^8]
- Interval arithmetic: Krawczyk test established by Rump[^9][^10], Rump sufficient criterion for positive definiteness[^11][^12], works up to condition $10^{15}$[^13], naive determinant factorial motivating specialised algebra[^14]

This version makes the referee see exactly where analysis ends and engineering begins, and where the dual certification loop closes the logical gap.

[^1]: Direct guaranteed lower eigenvalue bounds with optimal a priori convergence rates for the bi-Laplacian — https://arxiv.org/pdf/2105.01505
[^2]: Guaranteed lower bounds for eigenvalues - UBC Library Open Collections — https://open.library.ubc.ca/cIRcle/collections/48630/items/1.0377220
[^3]: A Two-Stage Finite Element Approach for High-precision Guaranteed Lower Eigenvalue Bounds — https://arxiv.org/pdf/2512.23182
[^4]: [2507.04307v3] Pólya's conjecture up to $ε$-loss and quantitative estimates for the remainder of Weyl's law — https://arxiv.org/abs/2507.04307v3
[^5]: Singular value asymptotics on compact smooth Riemaniann manifolds — https://arxiv.org/pdf/2512.02365
[^6]: Deep estimates for higher eigenvalues of the poly-Laplacian — https://arxiv.org/pdf/2508.04069
[^7]: On approximating the modified Bessel function of the second kind | Journal of Inequalities and Applications | Springer Nature Link — http://link.springer.com/article/10.1186/s13660-017-1317-z
[^8]: [1211.7325v2] Inequalities for modified Bessel functions and their integrals — https://arxiv.org/abs/1211.7325v2
[^9]: VERIFIED COMPUTATIONS FOR HYPERBOLIC 3-MANIFOLDS — http://arxiv.org/pdf/1310.3410v1
[^10]: rump.dvi — https://www.tuhh.de/ti3/paper/rump/Ru10.pdf
[^11]: Verified stability analysis using the Lyapunov matrix equation. - Free Online Library — https://www.thefreelibrary.com/Verified+stability+analysis+using+the+Lyapunov+matrix+equation-a0365070764
[^12]: METHODS FOR VERIFIED SOLUTIONS TO CONTINUOUS-TIME ALGEBRAIC RICCATI EQUATIONS TAYYEBE HAQIRI AND FEDERICO POLONI — https://arxiv.org/pdf/1509.02015v1
[^13]: reliable-computing-15-pp-193-207.dvi — https://jurnal.jumanji.workers.dev/host-http-interval.louisiana.edu/reliable-computing-journal/volume-15/no-3/reliable-computing-15-pp-193-206.pdf
[^14]: Verified bounds for the determinant of real or complex point or interval matrices1 — http://tuhh.de/ti3/paper/rump/Ru20.pdf
