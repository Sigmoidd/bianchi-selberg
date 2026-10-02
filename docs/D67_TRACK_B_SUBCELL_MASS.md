# Sharper finite mass integration on the same CR space

The interpolation-error constant continues to use the original safe m_T
upper bound. The finite mass matrix may use a sharper enclosure independently.

Partition each reference tetrahedron T into integration subtetrahedra S.
Let m_S be a rational pointwise upper bound for the pullback density on S.
For the affine parent CR basis phi define

\[
M_T^{sharp}=\sum_{S\subset T}m_S\int_S\phi\phi^T.
\]

Then the matrix upper bound holds for every affine w, even when individual
CR-basis products have negative integrals. For volume V and four vertex
value vectors f_i of the parent basis, the exact integral is

\[
\int_S\phi\phi^T=\frac V{20}
\left[(\sum_i f_i)(\sum_i f_i)^T+\sum_i f_i f_i^T\right].
\]

For red subdivisions all f_i are rational dyadics and V=|T|/8^depth.
Take H_min as the minimum of the four vertex heights, r_min as the minimum
vertex r coordinate, D_max=4-H_min and
g_min=(1-r_min)H_min+4r_min. Then m_S=D_max/(2g_min²) is a safe upper bound.
Concavity of H and monotonicity of g in H and r prove this bound.

Subcell bounds are no larger than the old parent bound, so
M_T^sharp <= m_T integral(phi phi^T). Subdivision changes neither the finite
space nor P, gamma or the inherited direct-index scalar budget. Young's mass
inequality uses M_h^sharp(Iv), while M(e) still uses the old m_T error envelope.

The finite matrix must still be assembled and independently enclosed with
exact rational or rigorous interval arithmetic. adaptive_matrix_probe.py is
a floating diagnostic and cannot supply those enclosures.
