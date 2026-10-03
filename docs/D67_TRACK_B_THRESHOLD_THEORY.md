# d=67 threshold-index exclusion and mapped-core CR lower bounds

Status: **conditional theorem; spectral exclusion OPEN**. This note supplies
an alternative exhaustive criterion to the frozen s-window design. It does
not supersede the geometry evidence or claim that its finite checks pass.

## 1. A single threshold form suffices

Use the exact glued core K from D67_TRACK_B_GEOMETRY.md and Y=2. Divide
all forms by A=sqrt(67)/2. Let t(v)=integral_P v(z,Y) da db and
b(v)=integral_P v(z,Y)^2 da db. Define

    C(v)=Q_K(v)/A-M_K(v)/A + alpha b(v) -(alpha+1/4)t(v)^2,
    alpha=1/2.

The core form domain consists of H1 functions with every quotient face
trace identification. No artificial Dirichlet condition is imposed.

**Threshold theorem.** If C has negative index at most one on this domain,
then the full level-one d=67 quotient has no L2 eigenvalue in (0,1).

Proof of cusp inequality. Write r=log(y/Y) and u_mu(y)=exp(r) w_mu(r)
for each torus Fourier coefficient. Integration by parts gives

    (Q_cusp-M_cusp)/A
      = Y^-2 sum_mu [ integral_0^infty (|w_mu'|^2
          + (2 pi |mu| Y)^2 exp(2r) |w_mu|^2) dr - |w_mu(0)|^2 ].

This identity holds first for compactly supported smooth cusp functions;
the resulting inequality extends to the energy form domain by density.
For mu=0 the integral is nonnegative. For mu!=0 put k=2 pi |mu|Y.
Using exp(2r)>=1 and completing the square,

    integral (|w'|^2+k^2|w|^2) dr >= k |w(0)|^2.

The dual lattice minimum is |mu|=2/sqrt(67), hence k_min=8 pi/sqrt(67)>3.
Thus (k_min-1)/Y^2>1/2, and Parseval gives

    (Q_cusp-M_cusp)/A >= -t^2/4 + (1/2)(b-t^2).

For a constant and any exceptional eigenfunction, their span E has
strictly negative global Q-M (orthogonality and lambda<1). Restriction
to K is injective on E: a vector vanishing there has top trace zero and
nonnegative cusp Q-M, contradicting strict negativity. The displayed
cusp bound implies C is strictly negative on the two-dimensional image
of E. This contradicts negative index at most one. The same argument
counts all eigenvalues below 1, including residual eigenfunctions.
There is no assumption about the behavior of a trial at an s-midpoint
or about excluding endpoints by numerical convergence. Essential spectrum
starts at 1 under the standard finite-volume spectral decomposition;
the argument itself excludes L2 eigenvalues in the requested open interval.

The constant function is a known negative direction for C. Define

    L(v)=integral_K v dV/A + t(v)/4,
    L(1)=Vol(K)/A+1/4.

Then C(1,v)=-L(v), and

    C(v-L(v)/L(1)) = C(v)+L(v)^2/L(1).

Consequently C is nonnegative on {L=0} if and only if
C+L^2/L(1) is nonnegative on the whole glued domain. Either condition
makes its negative index exactly one. It is enough to prove
C+rho L^2>=0 with any rho>=1/L(1) that actually passes the inequality.
The value of rho does not carry a mathematical assumption; its finite
matrix and error bounds do. Taking rho>1/L(1) allows a strict finite
margin in the constant direction.

Requiring positivity on {t=0} would also suffice but is much stronger:
a core function can have nearly zero top mean while approximating a
constant in the bulk. The first top-mean-penalty diagnostic fails that
stronger criterion. It is not evidence that C has two negative directions
or that the quotient has an exceptional eigenvalue.

A rational proof of k_min>3 uses pi>25/8, obtained by 4 times the
102-term alternating lower sum for arctan(1). Then
64(25/8)^2=625>9*67=603. No decimal value of pi enters the criterion.

## 2. Exact reference domain, without a floor lift

On each exact projected Ford patch set H=q_(c,l)(a,b), D=4-H, and

    Phi(a,b,r)=(a+b/2, sqrt(67)b/2, sqrt(g)),
    g=(1-r)H+4r, 0<=r<=1.

The patch prisms form a flat reference domain P x [0,1]. H>=2/67>0
and H<=1<4; Phi is continuous and bi-Lipschitz, piecewise smooth. Its
inverse is r=(y^2-H)/(4-H). Adjacent patch formulas agree on their
common edges. Thus pullback identifies H1(K) with H1 of the reference
box, including the matching quotient boundary traces.

The normalized mass weight is exactly

    m=D/(2g^2).

Let G^-1=(1/67)[[68,-2],[-2,4]], h=(1-r)(H_a,H_b)/D. The normalized
energy coefficient is the rational symmetric positive matrix

    B=D/(2g) [[G^-1, -G^-1 h],
               [-h^T G^-1, h^T G^-1 h + 4g/D^2]].

These formulas follow directly from the Jacobian of Phi. They require
neither numerical differentiation nor a curved-domain sliver estimate.
Although Phi contains square roots, B and m are rational functions of
rational coordinates and the field-specific quadratic H.

At r=0 the Ford face map is an affine real-plane isometry:

    z' = a_matrix/c - conjugate(c z+d_matrix)/c.

It preserves H and projected area. After the required lattice translation,
its Jacobian determinant has absolute value one in (a,b). At the vertical
sides the map is lattice translation, preserves H and r, and also preserves
reference face area. A compatible reference tetrahedral mesh can therefore
identify boundary face means exactly. It is essential that complete face
triangles map to complete face triangles; matching only vertices, midpoint
samples, or unrelated rows in a matrix does not establish this property.

The 47 refined floor cells and their centroid fans in reference_mesh.py
produce 228 exact rational triangles. Each maps bijectively onto a triangle
under its Ford pairing; periodic edge subdivisions agree. Uniform midpoint
refinement commutes with every affine pairing. Extruding identical r layers
and using a consistent three-tetrahedron prism split supplies the reference
mesh. A separate verifier must check coverage, all incidences, orientations,
pairing bijections and projected metric identities before certification.

## 3. Guaranteed mapped CR criterion

Let T be a reference tetrahedron in one patch prism. Obtain exact or
interval-verified constant matrices B_T>0 with B(x)>=B_T on T, and
scalars m_T>=m(x)>0. On the face-mean CR space with all paired boundary
face means tied, assemble

    Q_h(w)=sum_T integral_T grad(w)^T B_T grad(w),
    M_h(w)=sum_T m_T integral_T w^2,
    b_h(w)=sum_top_faces |F| w_F^2,
    t_h(w)=sum_top_faces |F| w_F.

Here w_F is the face mean, the CR degree of freedom. The positive boundary
term is the squared piecewise constant projection of the trace, NOT the
full tangentially varying CR trace. Using the latter without an error
bound would be unsound.

For the CR interpolant I and e=v-Iv, each element satisfies integral_T
grad(e)=0. Hence the constant-matrix energy splits exactly into Q_h(Iv)
and Q_h(e). On each reference tetrahedron, the standard self-contained
CR estimate has kappa^2=1/pi^2+1/15. Apply it after the linear transform
x=B_T^(1/2) xi to obtain

    integral_T e^2 <= kappa^2 h_(T,B)^2
                       integral_T grad(e)^T B_T grad(e),
    h_(T,B)^2=max_edges (p_i-p_j)^T B_T^-1(p_i-p_j).

Use the rational upper constant kappa^2<=8/45, following pi>3.
Let gamma^2=max_T m_T (8/45) h_(T,B)^2. For eta>0,

    M(v)<= (1+eta)M_h(Iv)+(1+1/eta)gamma^2 Q_h(e).

Jensen gives b(v)>=b_h(Iv); top face means give t(v)=t_h(Iv).
Paired affine area-preserving faces ensure Iv belongs to the tied CR space.
These properties remain valid at edges and stabilizers, which are measure
zero; their pointwise geometry must nevertheless be verified for the
quotient interpretation.

The core mean functional a(v)=integral m v has interpolation error

    |a(e)| <= sigma sqrt(Q_h(e)),
    sigma^2=sum_T m_T^2 |T| (8/45) h_(T,B)^2.

This follows by the local CR bound and Cauchy--Schwarz over elements;
the exact top mean has no error, so |L(e)| satisfies the same bound.
The finite functional L_h(w)=integral m w+t_h(w)/4 must use the actual
weight m, with independently enclosed element integrals. Replacing it
by a midpoint quadrature functional without controlling the error is
not permitted. For theta in (0,1), reverse Young gives

    L(v)^2 >= (1-theta)L_h(Iv)^2
                  -(1/theta-1)sigma^2 Q_h(e).

**Finite criterion.** For some eta,rho>0 and 0<theta<1, verify

    c_e=1-(1+1/eta)gamma^2-rho(1/theta-1)sigma^2 >= 0,
    Q_h +(1/2)b_h -(3/4)t_h^2
          +rho(1-theta)L_h^2 -(1+eta)M_h >= 0.

Then C(v)+rho L(v)^2>=0 for every glued-core v. Indeed subtract the
mass bound from the energy splitting, apply Jensen to the positive top
term and reverse Young to the positive rank-one term: the finite matrix
contributes nonnegatively, and so does c_e Q_h(e). On {L=0} this makes
C nonnegative. Any strictly negative subspace therefore has dimension
at most one, proving the threshold theorem. The extra condition
rho>=1/L(1) is a useful consistency check (testing the constant forces
rho L(1)>=1), not a substitute for the finite inequalities.

No field-specific matrix or scalar check is asserted here. Required final
evidence is reconstruction of a validated mesh; enclosures of B_T and m_T;
the computed gamma and sigma bounds; enclosures of L_h; a rigorous
positive-definiteness certificate for the displayed finite matrix; and
independent replay tied to the exact field, group, geometry and theorem
identities. The floating probe is solely a diagnostic and must not
promote spectral_exclusion_certified.

## 4. Exact coefficient envelopes

coefficient_bounds.py constructs a rational pointwise Loewner lower
matrix, avoiding a floating eigensolve for coefficient admissibility.
On the projected convex hull of each tetrahedron, concavity puts the
minimum of H at a vertex. Its maximum is attained at the stationary
center if contained, or at an edge stationary point or endpoint. The
finite clamp checks are exhaustive. Bounds on r and these H extrema
give D_min,D_max,g_min,g_max. Hence

    a=D/(2g)>=a_min=D_min/(2g_max),
    b=2/D>=b_min=2/D_max,
    m<=m_T=D_max/(2g_min^2).

Rational interval products enclose h=(1-r)(H_a,H_b)/D. With interval
midpoint h0 and component radii (r_a,r_b), set

    R^2=(68 r_a^2+4 r_a r_b+4 r_b^2)/67,
    a_T=a_min b_min/(b_min+2 a_min R^2), b_T=b_min/2.

The bound delta^T G^-1 delta<=R^2 and reverse Young prove

    B(x)>=B_T,
    xi^T B_T xi=a_T |xi_h-h0 xi_r|_(G^-1)^2+b_T xi_r^2.

B_T is strictly positive, and its inverse metric diameter is exactly

    max_edges [ (d_a^2+d_a d_b+17 d_b^2)/a_T
                   +(d_r+h0_a d_a+h0_b d_b)^2/b_T ].

All scalar arithmetic is rational. The accumulated sigma bound rounds
each nonnegative term upward to a multiple of 2^-40; this limits integer
size while preserving the bound. Float renderings of gamma and sigma
are diagnostics only; gate comparisons use the exact rational values.

Adaptive red refinement in adaptive_error_bounds.py establishes only
local coefficient/error bounds. Its leaf tetrahedra may have hanging
faces. These local bounds do not yet prove a globally compatible CR
interpolation space: nested face-moment aggregation, periodic and Ford
face propagation, and their matrix maps require separate proof and
reconstruction. An adaptive scalar pass must not be promoted to a mesh
or spectral certificate.

## 5. Hanging-face moment space (conditional mesh extension)

For a nonmatching reference mesh, first prove an exact common partition
of every two-sided or paired face into affine face fragments. Identify
matching fragments and assign each class one mean degree of freedom.
For a tetrahedron face F define its local mean by the area-weighted sum
of its fragments' means. The interpolation vector of an actual H1
function uses its exact mean on each fragment; continuity and affine
area-preserving gluing make that vector single-valued on all fragment
classes. Summing its fragment integrals reproduces every tetrahedron
face mean. Thus the local CR gradient identity, energy split, and error
bounds above survive unchanged, with a rational prolongation P from
fragment means to local face means.

The finite matrices must be P^T Q_local P and P^T M_local P; the top
terms and exact weighted L functional must use the same P. If P has a
kernel, remove it by an exact rank/basis calculation before a strict
positivity test. No projection map may be guessed from midpoint equality.
For recursive red refinement the face fragments are nested quarter
triangles, whose area weights are rational powers of 1/4. One still
must reconstruct their coverage, complete pairing orbits (including
self-pairings), incidences, and P. Local adaptive coefficient bounds
alone do not establish any of these conditions.
