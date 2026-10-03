# Adjustable vertical energy enclosure

This note refines the coefficient enclosure in
`D67_TRACK_B_THRESHOLD_THEORY.md`. It does not certify spectral exclusion.

On a reference tetrahedron write the pointwise energy as

\[
a(x)\|\xi_h-h(x)\xi_r\|_{G^{-1}}^2+b(x)\xi_r^2,
\qquad a(x)\ge a_{\min}>0,\quad b(x)\ge b_{\min}>0.
\]

Let the exact componentwise interval enclosure for h have midpoint h0,
radii ra,rb, and

\[
R^2=(68r_a^2+4r_ar_b+4r_b^2)/67.
\]

Then ||h-h0||² in the G⁻¹ metric is at most R². Fix a rational
0 < f < 1, and set s=(1-f)b_min. Reverse Young's inequality gives

\[
B(x)\succeq B_T,
\quad
\xi^T B_T\xi=
\frac{a_{\min}s}{s+a_{\min}R^2}
\|\xi_h-h_0\xi_r\|_{G^{-1}}^2+f b_{\min}\xi_r^2.
\]

Indeed, with z=ξ_h-h0 ξ_r and δ=h-h0, take
ε=a_min R²/(s+a_min R²). The inequality
||z-δ ξ_r||² ≥ (1-ε)||z||²-(1/ε-1)||δ||² ξ_r²
loses at most s ξ_r² after multiplication by a_min. If R=0,
the conclusion follows directly without dividing by ε.

The old enclosure uses f=1/2. The candidate f=3/4 retains more vertical
energy while reducing the horizontal coefficient. Both the mass bound and
the CR metric-diameter scalar bounds must be recomputed with this same B_T.
The old scalar ledger cannot be used unchanged with the new matrix.

`tuned_error_bounds.py` performs this exact rational scalar replay on every
leaf of the frozen plan, with κ²=1661/15000. Each sigma summand is rounded
upward to a multiple of 2⁻⁴⁰; gamma is an exact maximum. The mesh, plan and
historical ledger are preserved. The resulting scalar ledger is a prerequisite
for the finite matrix check, not evidence of its positivity.
