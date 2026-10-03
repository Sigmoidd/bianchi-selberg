# Rational replacement of the weighted core mean

This is a conditional refinement of the mapped-CR threshold criterion.
The scalar and finite-matrix conditions still require exact replay and verified
positivity. No spectral conclusion is asserted here.

Let w=Iv be the CR interpolant on the verified global moment space. Choose
any rational piecewise constant density m0 on the reference tetrahedra, and
write

\[
L_h(w)=\int m w+\tfrac14 t(w),\qquad
\widetilde L_h(w)=\int m_0 w+\tfrac14t(w).
\]

With the same mass upper envelope m_T used for M_h, a scalar enclosure

\[
\Delta^2\ge\sum_T\int_T\frac{(m-m_0)^2}{m_T}
\]

gives |L_h(w)-Ltilde_h(w)|² ≤ Δ² M_h(w) by Cauchy–Schwarz.
For any rational 0 < ν < 1, reverse Young then gives

\[
L_h(w)^2\ge (1-\nu)\widetilde L_h(w)^2
 -(\nu^{-1}-1)\Delta^2 M_h(w).
\]

Consequently the earlier finite condition is implied by positivity of

\[
Q_h+\tfrac12 b_h-\tfrac34t_h^2
 +\rho(1-\theta)(1-\nu)\widetilde L_h^2
 -\bigl[1+\eta+\rho(1-\theta)(\nu^{-1}-1)\Delta^2\bigr]M_h.
\]

The interpolation-error scalar condition c_e ≥ 0 is unchanged. The additional
mean-enclosure error is absorbed in the finite mass coefficient, rather than
being silently omitted from the rank-one mean penalty.

For a constant m0 on a tetrahedron, the mean vector has four local entries
m0 |T|/4 because each CR basis function has average 1/4. All entries can
therefore be rational. An interval [m_lo,m_hi] on a tetrahedron gives the safe
choice m0=(m_lo+m_hi)/2 and contribution
|T|(m_hi-m_lo)²/(4m_T) to Δ². Subdivision for integration improves this scalar
bound without changing the finite-element space: integrate the affine CR
basis exactly on each integration subtetrahedron, preserving the same local
coefficient vector. If m0 varies across those integration subtetrahedra,
the mean entries must use those exact subcell basis integrals.

This route avoids requiring exact closed forms for integrals of the rational
pullback density. It replaces them with explicit rational mean entries and an
explicit error penalty. A floating quadrature vector by itself remains
insufficient for a threshold certificate.
