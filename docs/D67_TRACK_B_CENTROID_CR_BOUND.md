# A self-contained sharp CR interpolation bound on arbitrary tetrahedra

This supplements the frozen threshold note; it changes no historical
certificate or error ledger. The derivation uses only the convex-domain
Poincare bound of Bebendorf, Theorem 3.2, DOI 10.4171/ZAA/1170.
The resulting constant agrees with Carstensen--Puttkammer §2.4.1,
arXiv:2203.01028. Neither a shape-regularity assumption nor a published
CR interpolation theorem is used in the derivation below.

Let T be any nondegenerate tetrahedron, its vertices p_0,...,p_3,
centroid c=(p_0+...+p_3)/4, and diameter h. For the CR interpolation
error e, every face integral is zero. Therefore the divergence theorem
for (x-c)e gives

    3 integral_T e = -integral_T (x-c).grad(e).

Indeed (x-c).n_F is constant on each flat face F, and its product with
the face integral of e vanishes. This applies to H1 functions by density
and continuity of the trace. For the element mean e_bar,

    |T| e_bar^2 <= [ integral_T |x-c|^2/(9|T|) ] ||grad(e)||_T^2.

The exact barycentric moments give

    integral_T |x-c|^2/|T|
      = sum_i |p_i-c|^2/20
      = sum_(i<j) |p_i-p_j|^2/80
      <= 3 h^2/40.

Thus |T| e_bar^2 <= h^2 ||grad(e)||_T^2/120. Orthogonality of the
mean and the mean-zero part, followed by convex Poincare, proves

    ||e||_T^2 <= (1/pi^2+1/120) h^2 ||grad(e)||_T^2.

The centroid is crucial: a vertex-based vector field gives a weaker
mean bound. The argument applies after any invertible affine change of
coordinates, so the metric-diameter variant in the threshold theorem
has the same constant. The elementary rational lower bound pi>25/8
already checked in that theorem yields the rational upper constant

    kappa^2 <= 64/625+1/120 = 1661/15000.

Compared with the frozen 8/45 bound, the exact ratio is 4983/8000.
Multiplying the frozen gamma^2 and sigma^2 upper bounds by this ratio
is safe: all local terms are proportional to kappa^2, and scaling an
upward-rounded nonnegative sum remains an upper bound. No stored hash
or historical error ledger needs to be changed to use this improvement.

Primary sources checked:
- Bebendorf's corrected all-dimensional convex Poincare theorem:
  https://ems.press/content/serial-article-files/35407?nt=1
- The matching CR constant in §2.4.1 (comparison only):
  https://arxiv.org/pdf/2203.01028
