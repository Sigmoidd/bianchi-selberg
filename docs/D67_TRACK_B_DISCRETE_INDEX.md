# Direct transfer of the negative index

This criterion replaces the weighted-mean stabilization as a sufficient
condition. It leaves a finite matrix positivity check, but needs no enclosed
weighted-mean vector. The spectral exclusion remains open until that matrix
check has an independent rigorous replay.

Use the exact global face-moment map P and the local coefficient and mass
upper envelopes from the threshold theory. For eta > 0 define the finite form

\[
C_h(w)=Q_h(w)+\tfrac12 b_h(w)-\tfrac34t_h(w)^2
 -(1+\eta)M_h(w).
\]

For v in the glued core H¹ space, let w=Iv and e=v-Iv. Constant-coefficient
energy splitting, the CR error bound and Young's mass inequality give

\[
C(v)\ge C_h(Iv)+c_0Q_h(e),\qquad
c_0=1-(1+\eta^{-1})\gamma^2.
\]

The top trace means are exactly reproduced; the positive trace term uses
Jensen on the master top-face fragments. No bulk mean occurs in this inequality.

If c0 ≥ 0, every strictly C-negative subspace maps injectively under I to a
strictly C_h-negative subspace. Indeed, for v in that subspace,
C_h(Iv) ≤ C(v)-c0 Q_h(e) < 0. If Iv=0 this is impossible. The same inequality
applies to each nonzero linear combination in the subspace. Therefore

\[
\operatorname{ind}_-(C)\le\operatorname{ind}_-(C_h).
\]

The continuum threshold theorem proves that index_-(C) ≤ 1 excludes every
eigenvalue in (0,1). It is thus sufficient to verify index_-(C_h) ≤ 1.

One particularly convenient sufficient matrix certificate is

\[
C_h+\alpha z z^T\succeq0
\]

for any explicitly specified rational finite vector z and rational alpha > 0.
If C_h had a two-dimensional negative subspace, that subspace would have a
nonzero intersection with ker(zᵀ); on this vector the added form is zero,
contradicting positivity. The vector z need not approximate a geometric mean
and no error bound for such an approximation is required.

A candidate that needs only rational assembled mass and top-area data is
z=M_h 1+t_h/4. This choice is a diagnostic until the full stabilized matrix
has been rigorously verified. Other rational vectors, including a rounded
computed negative mode, are equally valid choices. The scalar sigma is no
longer used in this criterion.

With eta=1/10 the scalar requirement is gamma² ≤ 1/11. The candidate tuned
three-quarter envelope has a floating diagnostic gamma²≈0.07092, suggesting
c0≈0.21991. The exact scalar ledger must replace this floating value, and the
2,042,568-variable finite matrix must still be certified.
