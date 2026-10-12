# Ford vertices and global odd trials

The exact vertex constraints are now explicit. They do not produce a spectral
certificate. `ford_vertex_equivariance.py` regenerates
`ford_vertex_equivariance.json` using rational arithmetic. The JSON contains
every stabilizer matrix, its two fiber-permutation conventions, exact
Jacobian, local Reynolds projector, simultaneous odd projector, and a finite
atlas incidence/transition ledger. No existing partition, floor, mass, or
defect certificate is changed.

The canonical source commit is
`b5f29bd1b314e4d385e5af01805729ba27e3f3e2`. The JSON records both raw SHA-256
and explicitly versioned `sha256-utf8-lf/v1` hashes for the geometry source,
geometry ledger, partition source, partition result, and primary/refined
partition audit leaves. These are new provenance records; historical
certificate digests and their policies are unchanged.

## Conventions and exact geometry

Write `t=y²`. The 24-dimensional first jet is

```
(F_0,...,F_5, dx1 F_0,dx2 F_0,dt F_0,...,dx1 F_5,dx2 F_5,dt F_5).
```

These coordinates keep every entry rational at all six vertices. A gradient
covector has squared hyperbolic norm `t*(dx1²+dx2²)+4*t²*dt²`, and
`dt F=(dy F)/(2y)`. The jet metric is the direct sum of the ordinary value
metric and six copies of `diag(t,t,4t²)`. Complex jets are obtained by extending
these rational matrices complex linearly; dimensions below are complex
dimensions (equivalently rational dimensions of the matrices).

For `g=[[a,b],[c,d]]`, put `D=|cz+d|²+|c|²t` and
`N=(az+b)conj(cz+d)+a*conj(c)*t`. Then `(z,t)` maps to
`(N/D,t/D²)`. Differentiating this formula gives the serialized rational
Jacobian. Determinants have norm one for all matrices used here.

The six points and a deterministic reduction `d_v Qv = w` are:

| Vertex | z | t | Qz | w | d_v |
|---|---|---|---|---|---|
| `v_00` | 0 | 1 | 0 | `v_00` | I |
| `v_0h` | i/2 | 3/4 | -1/2 | `v_mh0` | I |
| `v_mh0` | -1/2 | 3/4 | -i/2 | `v_0h` | R |
| `v_ph0` | 1/2 | 3/4 | i/2 | `v_0h` | I |
| `v_mhh` | (-1+i)/2 | 1/2 | (-1-i)/2 | `v_mhh` | T_i |
| `v_phh` | (1+i)/2 | 1/2 | (-1+i)/2 | `v_mhh` | I |

Here `R=diag(i,-i)`, `T_i=[[1,i],[0,1]]`, and `T_1=[[1,1],[0,1]]`.
The two additional closed-cell aliases are `T_1 v_mh0=v_ph0` and
`T_1 v_mhh=v_phh`. In particular, Q exchanges the two distinct order-six
Picard vertex classes; imposing an independent odd eigenspace of the fiber
permutation at each such vertex would be wrong.

All point stabilizers are exhaustively enumerated using the existing
Gaussian-entry denominator bound, not a word-search cutoff. At a fixed point,
`|cz+d|²+|c|²t=1` and `|-cz+a|²+|c|²t=1`; hence `|c|²<=2` and
`|a|,|d|<=2` at these vertices. The determinant determines `b` when `c!=0`;
when `c=0` the fixed-point equation determines it. Enumeration modulo `±I`,
exact fixed-point tests, closure, and all induced permutations are replayed.
Every stabilizer also has an explicit generator word recovered from the
existing depth-six BFS ledger and re-evaluated exactly against its matrix.
No new word search is used to infer completeness. Reducing words in the
table are literal `I`, `R`, or `Ti`; all six reducing transformations belong
to H_pi (their lower-left entry is zero). The words for the two aliases are
`T1`.

| Vertex | Stabilizer order | Fiber orbit sizes | Stabilizer value dimension | Stabilizer gradient dimension | Allowed odd value dimension at vertex | Allowed odd gradient dimension at vertex |
|---|---:|---|---:|---:|---:|---:|
| `v_00` | 4 | 2+2+2 | 3 | 3 | 1 | 1 |
| `v_0h` | 6 | 3+3 | 2 | 2 | 2 | 2 |
| `v_mh0` | 6 | 3+3 | 2 | 2 | 2 | 2 |
| `v_ph0` | 6 | 3+3 | 2 | 2 | 2 | 2 |
| `v_mhh` | 12 | 6 | 1 | 1 | 0 | 0 |
| `v_phh` | 12 | 6 | 1 | 1 | 0 | 0 |

## Simultaneous value and first-gradient equations

The repository's forward permutation is `p_g(c)=c.g`; its pullback
representation satisfies `(rho(g)v)_c=v[p_g(c)]`. The stored GLUE convention
is the inverse permutation. For `gp=p`, define the jet operator

```
A_g J = jet of rho(g)^(-1) F(gx) at p.
```

Its value block uses `p_g^(-1)` and its covector block uses the transpose of
the Jacobian. The stabilizer Reynolds projector is `P_G=sum A_g/|G|`.
It is idempotent and self-adjoint for the jet metric. This is the derivative
version of the orbifold averaging in `theorem_DK_sixcopy.tex`.

Let `q=(0,1,4,2,5,3)` and `A=d_v Q`, so `Av=w`. Oddness is

```
J_v + T_v J_w = 0,
T_v = jet pullback by A with fiber permutation p_d^(-1) o q.
```

This order matters. We use `rho_q rho(d)^(-1) F(dQx)`, not an involution
of the fiber alone. For a translation alias `gv=w`, the equation instead is
`J_v - A_g J_w=0`.

The simultaneous system has three components:

| Closed-cell vertices in component | Ambient jet dimension | Odd value dimension | Odd gradient dimension |
|---|---:|---:|---:|
| `v_00` | 24 | 1 | 1 |
| `v_0h,v_mh0,v_ph0` | 72 | 2 | 2 |
| `v_mhh,v_phh` | 48 | 0 | 0 |

Thus the three order-six entries in the first table are linked copies of
one four-dimensional jet space, not twelve independent parameters.
Every simultaneous projector is constructed as
`B*(B^T M B)^(-1)*B^T M`, where the columns of `B` are an exact nullspace
basis and `M` is the jet metric. The replay checks `P²=P`, `MP=P^T M`, and
`CP=0` for the complete constraint matrix `C`. The serialized local odd
projector describes the restriction of this simultaneous space to a single
vertex; cross-vertex equations remain necessary when using those local
projectors separately.

## Scalar nodal lines versus six-fiber jets

For the scalar H-component `F_0=f`, oddness and H-periodicity force
`f=0` and its **full** first gradient to vanish on every vertical line through
`z=0,(-1+i)/2,(1+i)/2`. At zero use Q. At the negative corner use `T_i Q`;
at the positive corner use `T_1 Q`. Each is an affine quarter turn fixing
that vertical line, has trivial transport on component zero, and acts on
tangent vectors with eigenvalues `i,-i,1`. Its differential has no `-1`
eigenvector. Differentiating `f(Ax)=-f(x)` therefore kills the entire first
jet on the line, for every positive height.
The two corner lines are H-translation-equivalent under T_1; these three
representatives therefore describe two quotient nodal geodesics.

This scalar argument does not kill the entire six-vector at `v_00`.
For example an allowed value vector there is `(0,0,1,-1,-1,1)`. An allowed
nonzero gradient jet has `dx1 F_2=-1`, `dx2 F_3=-1`, `dx2 F_4=1`,
`dx1 F_5=1`, with all other entries zero.

At the two corner vertices, the stronger full six-fiber vanishing **does**
hold after imposing the order-twelve isotropy as well. The unique
stabilizer-invariant value direction and the unique invariant gradient
direction are both even for reduced Q, so their odd intersections vanish.
This statement concerns the corner vertices; the scalar nodal-line argument
alone gives no such six-fiber assertion along the whole vertical lines.

## Which constraints are already in the Fourier basis?

The current odd coefficient basis has alternating signs on each C4 orbit at
both cusps. It already gives the exact bundle odd relation and the parabolic
relations (`T_1,R,T_iR`) through the two exact cusp lattices. Consequently the
Q-edge equations and the translation aliases add **no** independent Fourier
constraints. They should be checked as identities, not inserted as tiny
floating rows in a rank calculation.

For component zero, the scalar nodal jets above are also redundant
mode-by-mode. At zero and either half-integer corner, the four rotated
integer frequencies have the same phase. The alternating sum of these four
phases is zero, as is the alternating vector sum of their frequencies.
Their radial factors and radial derivatives agree. This proves value,
horizontal-gradient, and vertical-gradient vanishing without a numerical
rank tolerance.

The full isotropy equations involve inversion and are not implied by
cusp Fourier parity. For example take a single odd orbit of zero-cusp
frequency `beta=(1,0)`. At `v_mhh`, component one has value
`2*(cos(pi/5)-cos(3*pi/5))*radial = sqrt(5)*radial`, not an identically zero
function of the spectral parameter. Full corner vanishing is therefore a
real automorphy restriction. Genuine added diagnostic rows can be assembled
as `(I-P_G) J_v` on the existing odd coefficient basis. These finite point
and first-jet conditions are necessary for global smooth automorphy; they
are not sufficient for it, and their coefficient-space ranks depend on
cutoff and spectral parameter.

There is also an exact universal redundancy bound. Replace each local full
stabilizer by its subgroup of matrices with lower-left entry `c=0`, while
retaining Q and translation equations. Its simultaneous jet space has
dimension 20: dimensions 4, 12, and 4 on the three components above.
The full system leaves dimensions 2, 4, and 0. Thus the genuinely added
operator has rank **14** on this ambient parabolic odd jet space (component
ranks 2, 8, and 4), and its restriction to any existing odd Fourier cutoff
has rank at most 14. This statement is independent of the spectral
parameter. The replay computes the nullspace, applies the full constraint
matrix to that basis, and checks rank 14 with exact rational elimination.
A larger numerical rank indicates inconsistent jet coordinates, a trial
map failing one of these identities, or an unsuitable numerical threshold;
it cannot be attributed to new independent vertex equations.

`reduced_vertex_constraint_rows()` selects exactly 14 isotropy rows whose
restrictions to that 20-dimensional structural jet space are independent.
These sparse rows act on the full 144-vector of vertex jets and are
equivalent to all full vertex constraints modulo the structural identities.
The test compares the exact row-reduced systems. Numerical assembly should
use these rows before equilibration; row-normalizing roundoff-sized zero
rows from the unreduced system can fabricate extra numerical ranks.

Hard vertex constraints can exhaust a small truncated coefficient space.
That means no nonzero trial in that finite space satisfies those exact
jets at the tested parameter. It does not prove cutoff instability, exclude
an infinite Fourier solution, or produce a spectral certificate.

Two explicit kernel columns sharpen the upper bound at the small cutoffs.
The infinity orbit of `(1,1)` has horizontal factor
`-4 sin(2*pi*x1) sin(2*pi*x2)`. Its value and full first gradient vanish at
all six half-integer-coordinate vertices. The zero-cusp orbit
`beta=(3,1)=conj(pi)*(1+i)` has integer frequency `(1,1)` and the same
vanishing in every zero-cusp fiber, because integer translations do not
change it. They are independent coefficient columns and both occur by M=2.
At M>=4, the integer `(2,0)` orbit gives another such zero column at each
cusp. The implementation verifies these facts directly in Q(zeta_10),
checking every fiber and vertex; common radial factors and radial derivatives
make the conclusions independent of the spectral parameter.

| M | Odd coefficient count | Identically zero full-jet columns | Exact upper bound for added-constraint rank |
|---|---:|---|---:|
| 2 | 11 | 1, 9 | 9 |
| 3 | 13 | 1, 9 | 11 |
| 4 | 20 | 1, 2, 10, 19 | 14 |
| 5 | 25 | 1, 2, 12, 21 | 14 |
| 6 | 29 | 1, 2, 12, 21 | 14 |

Column indices are zero-based in `global_quarter_turn.sector_basis(M,-1)`.
These are upper bounds; equality requires a separate lower-rank proof.
The infinity diagonal column has zero coefficient in the prescribed
normalization mode `(1,0)`, so its nonzero kernel does not by itself supply
an admissibly normalized trial. Exact-zero columns must not acquire rank
from floating residuals followed by column normalization.

## Global projection and the atlas ledger

Pointwise equivariance of the existing partition is unnecessary to define
an odd global section. If a valid gluing construction has already produced
`W in Dom(Delta)` on the quotient, then

```
W_- = (W-U_Q W)/2
```

is globally well-defined, smooth where W is smooth, in `Dom(Delta)`, and
exactly odd. Normalization of Gamma and H by Q supplies compatibility on all
chart transitions. The standard domain/commutation argument is in
`GLOBAL_CUSP_SYMMETRY_THEOREM.md`, Theorems B1-B4 and F1-F3.

Equivalently, double a valid atlas using `U_j` and `Q^-1 U_j`, with weights
`chi_j/2` and `(chi_j o Q)/2`. Their sum is one. The transported candidate
on the second chart is `-U_Q F_j`; the resulting partition sum is exactly
`P_- W`. Isometry covariance gives the same gradient and Laplacian bounds
for `chi_j o Q`; halving each of the paired weights preserves the old
conservative sum-of-suprema bound. No numerical re-integration of the
existing 192-bit partition derivative certificate is needed for this
isometric transfer statement. It does **not** transfer a residual or mass
certificate to a different trial.

More explicitly, write `B0=sum_j |Delta chi_j|` and
`B1=2 sum_j |grad chi_j|`. On the doubled atlas,
`B_l,double(x)=(B_l(x)+B_l(Qx))/2` for l=0,1, with gradient norms understood
hyperbolically. This need not equal the original pointwise B_l(x).
Nevertheless a uniform bound is preserved, and the L2 norm is no larger
by the triangle inequality and the unitary isometry Q. The existing
certificate conservatively sums a per-label derivative majorant over all
33 labels, including labels that Reynolds averaging can activate. Doubling
gives 66 labels, each original/transformed member carrying factor 1/2;
their total majorant is still the same sum over 33 original bounds.
Stabilizer Reynolds sums have coefficients 1/|G| and are convex averages,
so they introduce neither an additional stabilizer-order factor nor extra
integration leaves. We do not assume an unchanged active-label count.
Local finiteness is preserved by taking a finite union of a valid cover and
its isometric image. These statements concern derivative bounds and a
valid input cover; they do not fill missing original or mixed overlap
domains in the current transition-defect ledger.

The JSON atlas ledger records all five face transitions, all eleven listed
one-dimensional singular strata, all six exact vertex groups, their fiber
transport, and all Q-conjugates. The finite permutation intertwiner is
checked on every serialized transition. These are incidence and word data;
they do not supply a complete set of continuum overlap boxes. Mixed
original/Q overlaps in the doubled atlas also require coverage if that
atlas is used directly in a new commutator calculation.

Accordingly `complete_continuum_transition_cover`,
`all_transition_defects_covered`, `odd_trial_mass_certified`, and
`odd_trial_residual_certified` are false. The public
`require_complete_transition_coverage` gate rejects this ledger. One may
instead bound the residual of `P_- W` by that of an already certified global
W, by contraction; its projected mass still needs a same-trial certificate.
No spectral enclosure, spectral gap, nonzero odd trial, newform assertion,
or existence flag follows from the vertex algebra or the atlas theorem.

## API and reproduction

- `vertex_projector(name)`: exact 24-by-24 stabilizer projector.
- `allowed_odd_vertex_projector(name)`: exact projector onto the allowed local odd jet image.
- `q_relation(name)`: target name and exact transport in `J_source+T J_target=0`.
- `alias_relation(source,target,matrix)`: transport in `J_source-T J_target=0`.
- `constraint_rows(names)`: sparse rational simultaneous constraints in the supplied vertex order.
- `simultaneous_projector(tuple_names)`: exact simultaneous odd projector.
- `reduced_vertex_constraint_rows()`: 14 independent added sparse rows on all 144 vertex jets.
- `structurally_zero_odd_jet_columns(M)`: exact Fourier columns whose full vertex jets vanish identically.

Run `python ford_vertex_equivariance.py` and
`python -m pytest -q tests/test_ford_vertex_equivariance.py`.
The tests replay serialization, exact stabilizer dimensions, off-stratum
Jacobian chain rules, all Q reductions, Reynolds invariance, simultaneous
projector identities, scalar Fourier identities, and fail-closed coverage.
