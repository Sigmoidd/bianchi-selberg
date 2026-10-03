# d=43 Track B exclusion design — frozen v1

Status: DESIGN ONLY; no exclusion certificate. Base: `86048a4` on
`codex/d19-certificate`. Work branch: `track-b/d43-design`.
The fields d=43,67,163 remain open. No numerical campaign is authorized by
this design freeze. Geometry is the first unproved input.

## 1. Target and conventions

Let tau=(1+i sqrt(43))/2, O=Z[tau], tau^2=tau-11, D=-43, and
Gamma=SL_2(O)/{+I,-I}. The exact software identity is
`GroupKey(43, "full", (1,0,1))`; no congruence subgroup, induced bundle,
symmetry sector or six-copy lift is substituted for this full group.

On H^3={(x1,x2,y):y>0}, use curvature -1 metric
ds^2=(dx1^2+dx2^2+dy^2)/y^2 and dV=y^-3 dx1 dx2 dy.
The positive self-adjoint scalar Laplacian is
Delta=-y^2(partial_x1^2+partial_x2^2+partial_y^2)+y partial_y.
Its closed energy and mass forms are Q(u)=integral |grad_euc u|^2/y
and M(u)=integral |u|^2/y^3. Functions on the quotient satisfy all face
identifications and local orbifold invariance. Complex eigenfunctions
may be reduced to a nonzero real or imaginary part.

The target is every discrete L^2 eigenvalue in (0,1), including residual
and cuspidal eigenfunctions. Constants give the simple eigenvalue zero.
The essential-spectrum threshold is 1, conditional here on the standard
finite-volume spectral theorem and the geometry proof below.
Parameterize the exceptional interval by lambda=1-s^2, 0<s<1.
Neither lambda=1+r^2 for real r nor a scan of cusp forms covers this target.

## 2. Method choice and final criterion

Choose compact-core cusp reduction followed by guaranteed nonconforming
finite-element lower bounds. Hejhal probes may later guide refinement;
they are diagnostic and are not part of the exclusion theorem.

Let F be a proved measurable fundamental cell, with disjoint interiors
modulo Gamma, and choose its infinity end so that F above Y is exactly
P x (Y,infinity), where
P={a+b tau: -1/2<=a,b<=1/2} is a translation cell of O.
Freeze Y=2 for the first design; changing Y requires a new design version.
Set K=F intersect {y<=2}, and C=P x (2,infinity). This specifies the
truncation but does NOT assert that the as-yet-unconstructed F or K has
proved geometry. Require a finite cell decomposition of K with positive
height lower bound, piecewise Lipschitz boundary and a conforming top
torus. Multiple cells and identified faces are allowed. A single spherical
floor graph over P is not assumed.

Write A=area(P)=sqrt(43)/2, t(v)=integral_P v(x,Y) dx,
a(v)=integral_K v dV, and

    L_s(v)=a(v)+t(v)/((1+s)Y^2),
    B_s(v)=Q_K(v)-(1-s^2)M_K(v)-(1-s)t(v)^2/(A Y^2).

Let V be the H^1 form domain on the actual compact quotient, including
paired-face traces and stabilizer invariance. An explicitly proved larger
space may replace V for a sufficient relaxation; the Gaussian Neumann
relaxation is not presumed to work for d=43.

FINAL EXCLUSION CRITERION: for every s in (0,1),
B_s(v)>0 for every nonzero v in V with L_s(v)=0.
A sufficient unconstrained pencil is, for a fixed nu>1 and rho>0,

    Q_K(v)+rho L_s(v)^2 >= nu[(1-s^2)M_K(v)
                                    +(1-s)t(v)^2/(A Y^2)]

for all v in V and all s in (0,1). Parameters nu,rho, the mesh and
error constants remain unset; no Gaussian values are adopted.

## 3. Arithmetic and cusp derivations

The norm is N(a+b tau)=a^2+ab+11b^2=(a+b/2)^2+43b^2/4.
It is a positive integer on nonzero O, so |c|>=1 for c!=0.
Solving N=1 gives only units +/-1. For a cusp stabilizer matrix c=0,
ad=1 forces a=d=+/-1; in PSL its action is z->z+ell, ell in O.
Thus there is no rotational quotient and no area factor 1/2 beyond
A=sqrt(43)/2 itself. Translation matrices are
T_1=[[1,1],[0,1]] and T_tau=[[1,tau],[0,1]], both determinant one.
They are cusp generators, not a proposed generating set for Gamma.

The reduced-form enumeration for D=-43 has a<=floor(sqrt(43/3))=3.
The reduction, parity and divisibility conditions leave only (1,1,11).
Via the classical ideal-class/binary-form correspondence, O is a PID.
This implication is an external arithmetic theorem, not Euclidean
division. Each projective cusp has a primitive pair (p,q); the unit-ideal
identity supplies r,t with pt-qr=1, yielding a matrix sending infinity
to p/q. Therefore there is one cusp. Replay should independently check
the enumeration and each supplied Bezout witness; no nearest-point
division is allowed.

For gamma=[[a,b],[c,d]], the exact height formula is
y(gamma P)=y/(|cz+d|^2+|c|^2 y^2).
If c!=0 and y>1 it is <=1/y<1. Hence the horoball y>1 is precisely
invariant under Gamma_infinity, and Y=2 is admissible for its cusp chart.
This does not prove a finite fundamental domain for its complement.

Under Euclidean pairing Re(mu conjugate(z)), the dual frequencies are

    mu_(m,n)=(m,(2n-m)/sqrt(43)),  (m,n) in Z^2,
    |mu|^2=m^2+(2n-m)^2/43.

Indeed pairing with 1 gives m and pairing with tau gives n.
If m=0 and n!=0 this is >=4/43, with equality for n=+/-1;
if m!=0 it is >=1>4/43. The sharp shortest squared frequency is 4/43.
No Gaussian integer-frequency cutoff transfers.

For nonzero modes, the energy defect at lambda<1 is bounded below by
(4 pi^2 (4/43)Y^2-lambda) times their cusp mass, which is positive at
Y=2 (64 pi^2/43>1). This follows from y^-1>=Y^2 y^-3 and dropping
nonnegative vertical energy. Fourier normalization uses a full torus
of area A; coordinates (a,b) have Jacobian A and the inverse Gram
metric, not the Cartesian unit-square metric.

For lambda=1-s^2 the zero mode solves y^2 f''-y f'+lambda f=0.
L^2 excludes y^(1+s), leaving f=c(y/Y)^(1-s), c=t/A.
Its cusp mass is A c^2/(2sY^2); its integral against the constant is
A c/((1+s)Y^2); its energy defect is
-A(1-s)c^2/Y^2=-(1-s)t^2/(A Y^2).
These identities produce precisely L_s and B_s above.
For an eigenfunction with lambda>0, orthogonality to constants gives
L_s=0; the full energy identity and nonnegative nonzero-mode defects
give B_s<=0. Unique continuation gives a nonzero restriction to K.
The final criterion contradicts this. It covers residual spectrum too.
No omitted Fourier tail or assumed coefficient envelope enters this
reduction: the inequalities apply to the full H^1 cusp expansion.

## 4. Dependency ledger and geometry gate

| Input | Current status | Required evidence |
|---|---|---|
| Field identity and PID | exact norm and reduced-form enumeration above; correspondence dependency | arithmetic replay and pinned correspondence theorem |
| Full finite-volume quotient and cusp | conditional on classical Bianchi lattice theory | pinned theorem or field geometry proof |
| Fundamental cell | UNPROVED | exhaustive Ford/Dirichlet construction; finite termination, coverage and nonoverlap |
| Pairings and isotropy | UNPROVED | all face matrices in SL_2(O), inverse/cycle relations, edge/vertex links and stabilizers |
| Precisely invariant cusp and dual lattice | derived above | exact replay and rigorous pi comparison |
| Core exclusion reduction | conditional proof above | F/K hypotheses, self-adjoint energy/orthogonality, unique continuation |
| Essential threshold 1 | standard dependency; not source-verified in this freeze | spectral theorem with precise hypotheses and reference location |
| CR lower-bound theorem | template only | adaptation to actual cell topology, identified traces, curved boundaries and field weights |
| Mesh and quadrature | absent | validated maps, volumes, inclusion/complement bounds, error constants |
| All-window positivity | absent | continuum lower-bound bridge and verified finite matrix inequalities |
| Independent reconstruction | absent | regenerate geometry and entries from exact inputs; reject identity mismatches |

A finite list of hemispheres or determinant checks alone proves neither
coverage nor completeness. For non-norm-Euclidean O, do not assume the
usual translations and inversion generate the group or that unit-radius
hemispheres cover the floor. A Ford construction must bound omitted
unimodular rows and establish termination. If this cannot be done, keep
the gate false and continue arithmetic and source audits only.

The Gaussian `independent_exclusion/m3_certify.py` hardcodes area 1/2,
floor sqrt(1-x1^2-x2^2), a rectangular domain and y_min=1/sqrt(2).
Changing a field label does not adapt any of them. Its asserts are not
accepted as fail-closed guards in a future verifier under python -O.

Additional audit blocker in `lower_bound_theory.md`, Lemma G, verbatim:
"≤ d²/4 (maximized at an edge midpoint of the longest edge)".
This asserted bound for barycentric variance is false on an arbitrary
triangle: an equilateral triangle of side d has variance d^2/3 at its
centroid. A universal triangle bound follows from
sum_i theta_i |p_i-p|^2=sum_(i<j) theta_i theta_j |p_i-p_j|^2
<=d^2(1-sum_i theta_i^2)/2<=d^2/3.
Consequently the corresponding universal Taylor sag bound is H d^2/6,
not H d^2/8. A stronger bound might hold for a particular triangulation,
but it requires its own proof. No historical file or certificate is
altered by this observation, and this freeze does not settle the historical
Gaussian implementation. D43 must prove its own floor inclusion.

## 5. Finite certification and replay contract

After geometry passes, use an identified CR space or a proved sufficient
relaxation. Audit the arbitrary-tetrahedron CR interpolation lemma,
weighted energy projection, trace reproduction, volume-functional errors,
and all curved-core remainder estimates. CR face-mean interpolation
reproduces t exactly only when the actual top face is exactly tiled;
orbifold and paired-face constraints need compatible interpolation.
Conforming Galerkin positivity supplies upper bounds and cannot close
the target. The continuum-to-discrete lower-bound theorem is mandatory.

Use finitely many rational closed s-windows whose union is [0,1]. On
[s_lo,s_hi], enclose lambda from above by 1-s_lo^2, the zero-mode
coefficient by (1-s_lo)/(A Y^2), and the changing constraint coefficient
1/((1+s)Y^2) by an interval/Young-shift bound about a fixed rational
reference. Account for this constraint change; checking midpoints alone
is invalid. Endpoint windows cover both lambda->0 and lambda->1; endpoint
identities need not assert L^2 membership of the s=0 zero-mode solution.
Strictness in the target follows from (nu-1)lambda M>0 for 0<lambda<1.

For every window report all scalar/error checks and verified matrix
positivity, including failures. Reconstruct all entries from field,
geometry, mesh and quadrature definitions; do not trust stored matrices.
Require explicit exceptions for missing inputs, unsupported fields,
incomplete window cover, invalid geometry and mismatched hashes. Bind
field=43, D=-43, unit level, metric, Y, design hash, arithmetic inputs,
geometry, pairings, mesh, theorem version, windows, precision, software
and any IEEE/LAPACK assumptions. A serialized `certified=true` is never
itself evidence. Exit nonzero on a failed required condition.

Sensitivity checks will delete a face, corrupt a pairing/stabilizer,
change area or dual basis, omit the constant constraint or boundary term,
understate an error radius, remove an endpoint window, substitute field
67 or Gaussian artifacts, and run with python -O. Each must either fail
before a claim or demonstrably invalidate the criterion. Freeze a passing
baseline before running mutations. The trace B harness remains unchanged.

## 6. Preservation and carry-forward

Use `artifacts/track_b/d43/` for new records; d67 and d163 will get
separate directories and branches. Existing certificates, historical
reports, Gaussian JSON/JSONL files and known failing tests remain intact.
Baseline commands and exact stdout/stderr belong in `baseline/`; a
timeout is a recorded failure, never a pass. Dependency installation is
from `requirements-trace.txt`. No mutation or expensive FEM run follows
a failing frozen replay.

To extend to d67 or d163, reuse only the general cusp/energy argument,
audited local interpolation theory and validated positivity machinery.
Reprove field geometry, pairing completeness, lattice and shortest dual
frequency, all core/mesh/error constants and all-window positivity.
For d=67,163 the corresponding squared shortest dual frequencies are
4/67 and 4/163; even though Y=2 makes the elementary nonzero-mode bound
positive for both, their designs must verify this explicitly. Neither
the d43 cell nor a successful d43 matrix transfers. Current target
statuses: d43 OPEN, d67 OPEN, d163 OPEN.
