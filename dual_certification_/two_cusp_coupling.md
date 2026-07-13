# Two-cusp coupling under \(\sigma_0\) (Rung 3 theory)

**Status:** analytic setup for AgentReady Rung 3  
**Group:** \(\Gamma=\mathrm{PSL}(2,\mathbb{Z}[i])\), congruence \(\Gamma_0(\mathfrak{p})\) for \(\mathfrak{p}=(2+i)\), \(N\mathfrak{p}=5\), index \(6\), two cusps \(\infty\) and \(0\).

---

## 1. Cusps and conjugating isometry

The two cusps of \(\Gamma_0(\mathfrak{p})\backslash\mathbb{H}^3\) are represented by \(\infty\) and \(0\). The matrix
\[
\sigma_0
=
\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
\in\mathrm{PSL}(2,\mathbb{C})
\]
interchanges them: \(\sigma_0(\infty)=0\), \(\sigma_0(0)=\infty\). On \(\mathbb{H}^3\) (upper half-space model \((z,y)\), \(z\in\mathbb{C}\), \(y>0\)) one has the standard inversion formula
\[
\sigma_0\cdot(z,y)
=
\Bigl(
  \frac{-\overline{z}}{|z|^2+y^2},\;
  \frac{y}{|z|^2+y^2}
\Bigr).
\]
(Cf. EGM; same \(\sigma_0\) as in `independent_exclusion/CONGRUENCE.md` Fact B.)

**Cusp stabilizers (parabolic).**  
- \(\Gamma_\infty\cap\Gamma_0(\mathfrak{p})\) acts by translations by the full lattice \(\Lambda_\infty=\mathbb{Z}[i]\) (folded by units in the orbifold picture).  
- The conjugate stabilizer at \(0\) is \(\sigma_0\Gamma_\infty\sigma_0^{-1}\cap\Gamma_0(\mathfrak{p})\), which acts by translations by the ideal lattice \(\Lambda_0=\mathfrak{p}=\!(2+i)\mathbb{Z}[i]\) (up to units). Euclidean areas of the folded cross-sections satisfy \(|T_\infty|=1/2\) and \(|T_0|=N\mathfrak{p}\cdot|T_\infty|=5/2\) (level-1 fold conventions as in the Gaussian congruence notes).

---

## 2. Fourier expansions at the two cusps

A cuspidal Maass form of eigenvalue \(\lambda=1+r^2\) admits Whittaker–Fourier expansions
\begin{align}
f(z,y)
&=
\sum_{0\neq\beta\in\Lambda_\infty}
a_\beta^{(\infty)}\,
y\,K_{ir}(2\pi|\beta|y)\,
e^{2\pi i\,\mathrm{Re}(\beta\overline{z})}
&&(y\gg 1),\label{eq:exp-inf}\\
(f\circ\sigma_0)(z,y)
&=
\sum_{0\neq\mu\in\Lambda_0^\vee}
a_\mu^{(0)}\,
y\,K_{ir}(2\pi|\mu|y)\,
e^{2\pi i\,\mathrm{Re}(\mu\overline{z})}
&&(y\gg 1),\label{eq:exp-0}
\end{align}
with \(\Lambda_0^\vee\) the dual lattice to \(\Lambda_0\) under the pairing \(\mathrm{Re}(\mu\overline{z})\) (for \(\mathfrak{p}=(2+i)\) one may identify modes with a \(\mathfrak{p}\)-scaled copy of \(\mathbb{Z}[i]\)).

**Truncation.** Fix \(M\ge 1\) and \(Y_0\ge 1/2\). Write \(a^{(\infty)}\) (resp. \(a^{(0)}\)) for the finite vectors of coefficients with \(0<N(\beta)\le M\) (resp. \(0<N(\mu)\le M\)).

---

## 3. Consistency as distributions on \(\{y=Y\}\)

Fix a truncation height \(Y\ge Y_0\). On the collar cross-sections \(T_\infty\times\{Y\}\) and \(T_0\times\{Y\}\) the two expansions define tempered distributions (finite Fourier sums after truncation, plus a tail controlled by Lemma K under Assumption H).

**Pull-back identity.** For a true eigenform, \(f=f\circ\gamma\) for all \(\gamma\in\Gamma_0(\mathfrak{p})\). In particular, taking \(\gamma\) ranging over a set of representatives that realize the cusp change through \(\sigma_0\) (modulo \(\Gamma_0\)), one obtains on the interface
\begin{equation}
\label{eq:consistency}
  f\big|_{T_\infty\times\{Y\}}
  \;=\;
  \bigl(f\circ\sigma_0\bigr)\big|_{\sigma_0^{-1}(T_\infty\times\{Y\})}
\end{equation}
in the sense of distributions on the flat tori (after the height-preserving identification of Fact B in the congruence notes: chimneys at cusp \(0\) assemble at the same height \(Y\)).

**Proposition (consistency under \(\sigma_0\)).**  
Let \(f\) be a true \(L^2\) cusp form on \(\Gamma_0(\mathfrak{p})\backslash\mathbb{H}^3\). Then the pair of coefficient sequences \((a^{(\infty)},a^{(0)})\) is consistent with \eqref{eq:consistency}: the distributional Fourier coefficients of the left- and right-hand sides of \eqref{eq:consistency} agree mode by mode. In particular, writing \(\mathcal{F}_\infty\) and \(\mathcal{F}_0\) for Fourier transforms on the two tori at height \(Y\),
\[
  \mathcal{F}_\infty\bigl(f|_{y=Y}\bigr)
  =
  S(r)\,\mathcal{F}_0\bigl((f\circ\sigma_0)|_{y=Y}\bigr)
\]
for an explicit linear operator \(S(r)\) determined by the Whittaker transform of the \(\sigma_0\)-pullback (the **scattering / coupling block**).

*Proof sketch.* Periodicity under \(\Lambda_\infty\) and \(\Lambda_0\) gives Fourier expansions \eqref{eq:exp-inf}–\eqref{eq:exp-0}. The isometry identity \(f=f\circ\sigma_0\circ\sigma_0\) and the change of variables under \(\sigma_0\) conjugate the two expansions; projecting onto characters of each torus yields a linear relation between coefficient vectors. The radial Whittaker factors transform among themselves (same \(K_{ir}\)), so the relation is \(r\)-dependent but \(y\)-independent once both sides are evaluated at the common height \(Y\) after the height-preserving chimney maps of Fact B. \(\square\)

(The computational matrix \(S\) below is a collocation discretization of this operator at height \(Y_0\), with entrywise radii tracked from the \(K\)-Bessel majorant of Lemma K.)

---

## 4. Coupled collocation matrix and gluing uniqueness

### 4.1 Continuous gluing relation

After truncating tails (Lemma K) and discarding continuous-spectrum packets (cuspidal subspace), the interface matching becomes
\begin{equation}
\label{eq:gluing}
  a^{(0)}
  =
  S(r)\,a^{(\infty)},
\end{equation}
or, allowing residual from truncation and automorphy defect,
\[
  a^{(0)}-S(r)a^{(\infty)}=\rho,\qquad \|\rho\|\le C_{\mathrm{tail}}+C_{\mathrm{aut}}.
\]

**Proposition (uniqueness of the gluing map).**  
For each fixed \(r\ge 0\), the linear map \(S(r)\) realizing \eqref{eq:gluing} on the space of coefficient sequences of rapid decay is **unique**. Indeed, if \(S_1,S_2\) both reproduce the Fourier coefficients of every true cusp form’s \(\sigma_0\)-pullback, then \(S_1-S_2\) annihilates a dense set of coefficient sequences (Hecke eigenforms span a dense subspace of the cuspidal \(L^2\) spectral packet at parameter \(r\) when the packet is nonempty; when empty the relation is vacuous). On the finite-dimensional truncated spaces used computationally, uniqueness of the minimal-norm \(S\) interpolating the collocation equations is the standard uniqueness for a full-rank least-squares map (checked numerically by \(\sigma_{\min}(S)>0\) after preconditioning).

### 4.2 Discrete two-cusp block system

Let \(V_\infty(r)\) (resp. \(V_0(r)\)) be the single-cusp collocation / automorphy matrix at cusp \(\infty\) (resp. \(0\)) for the truncated mode set of size \(n\), and let \(S(r)\) be the \(n\times n\) coupling block. The **coupled Hejhal system** is the \(2n\times 2n\) block matrix
\begin{equation}
\label{eq:coupled}
  \mathcal{V}(r)
  :=
  \begin{pmatrix}
    V_\infty(r) & -S(r) \\
    -S(r)^* & V_0(r)
  \end{pmatrix},
  \qquad
  \mathcal{V}(r)
  \begin{pmatrix} a^{(\infty)} \\ a^{(0)} \end{pmatrix}
  =
  0
\end{equation}
in the homogeneous setting, or \(=\mathbf{b}\) with a pinning / residual right-hand side in the inhomogeneous setting used by Krawczyk.

**Block-diagonal preconditioner.**
\[
  D
  =
  \mathrm{diag}\bigl(D_\infty,D_0\bigr),
  \qquad
  (D_\infty)_{\beta\beta}=w_\beta^{(\infty)}(Y_0),\quad
  (D_0)_{\mu\mu}=w_\mu^{(0)}(Y_0),
\]
with \(w\) the mode amplitude \(|K_{ir}(2\pi|\cdot|Y_0)|\,N(\cdot)^{\theta}\) (Lemma K majorant). Study \(\kappa(D^{-1}\mathcal{V})\) or column-normalized \(\mathcal{V}D^{-1}\).

### 4.3 N=5 specialization

For \(\mathfrak{p}=(2+i)\):
- \(R=\mathbb{Z}[i]/\mathfrak{p}\cong\mathbb{F}_5\), \(i\mapsto 3\), \(|\mathbb{P}^1(R)|=6\);
- residue gluing perms are those of `congruence_prototype.set_level('(2+i)')`;
- computational mode sets: \(\beta\in\mathbb{Z}[i]\) with \(0<N(\beta)\le M\) (cusp \(\infty\)); \(\mu\in\mathbb{Z}[i]\) with \(0<N(\mu)\le M\) as a dual/scaled model of \(\Lambda_0^\vee\) (same cardinality after truncation — exact dual indexing is a constant lattice automorphism and does not change \(\kappa\) asymptotics).

---

## 5. What this file claims / does not claim

| Claimed | Not claimed |
|---------|-------------|
| Consistency of true eigenforms under \(\sigma_0\) as distributions on \(\{y=Y\}\) | Full spectral completeness without Friedman/EGM (Ass. S) |
| Uniqueness of the continuous gluing map \(S(r)\) | That a nontrivial kernel of \(\mathcal{V}(r)\) exists at a given numerical \(r\) without a residual solve |
| Explicit block form \(\mathcal{V}(r)\) for computation | Certified dual eigenvalue (Rung 4) |
| Preconditioning slope \(b<4\) (float diagnostic) | That float \(\kappa\) equals the true continuous operator’s \(\kappa\) |

Assumption **H** enters only through Lemma K tails feeding radii of \(S\) and residual bounds. Assumptions **A, S** underwrite analyticity and the spectral decomposition used to speak of cusp forms.

---

## 6. Code map

| Artifact | Role |
|----------|------|
| `two_cusp_hejhal_N5.py` | Build \(\mathcal{V}\), \(S\) mid/rad, preconditioner, Krawczyk for N=5 |
| `hejhal_conditioning.py` | Milestone 2 slope diagnostics (already \(b_{\mathrm{eq}}<4\)) |
| `rung3_certificate.md` | AgentReady checklist |
