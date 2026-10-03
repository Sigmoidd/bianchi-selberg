# d=43 full-group geometry certificate

Status: **geometry gate proved**, for the exact full group
`PSL_2(Z[(1+sqrt(-43))/2])` and the translation cell fixed in
`D43_TRACK_B_DESIGN.md`. This certifies the fundamental region, its
lower and vertical faces, explicit face maps, and pointwise orbifold
stabilizers. It does **not** certify the spectral exclusion target;
compact-core FEM theory and all-window positivity remain open.

Reproduce from the repository root:

```sh
python artifacts/track_b/d43/ford_exact.py
python artifacts/track_b/d43/geometry_certify.py
python artifacts/track_b/d43/verify_geometry.py
python -O artifacts/track_b/d43/verify_geometry.py
python artifacts/track_b/d43/negative_checks.py
```

`cover_witnesses.json` has 322 exact rational squares and their
unimodular hemisphere witnesses. `geometry_result.json` has the
complete 19-face lower boundary ledger (every exact matrix and polygon),
56 edges with pointwise stabilizers, and the 38 vertex stabilizer groups
(all matrices modulo ±I).
The verifier reconstructs all 1,011 candidate rows, clips the polygons
against every competing hemisphere, checks image pieces of each paired
face against all rows, and independently enumerates vertex stabilizers.
All decisions are rational and remain active under `python -O`.

## Fundamental-domain argument

Write `tau=(1+sqrt(-43))/2` and `N(a+b tau)=a²+ab+11b²`.
Let `P={z=a+b tau: -1/2<=a,b<=1/2}`. For each unimodular pair
`(c,l)` in `O²`, `c!=0`, define

    q_(c,l)(z) = 1/N(c) - |z-l/c|².

Let `H(z)=max(0,sup q_(c,l)(z))`. The proposed region is

    F={ (z,y): z in P, y>0, y²>=H(z) }.

An orbit has a point of maximal height: for a fixed `P=(z,y)` and
`c!=0`, the image height is at most `1/(N(c)y)`. For any positive
height threshold this bounds `N(c)`; then
`|cz+d|²+N(c)y²` bounds `d`, leaving finitely many lower rows.
Translations do not alter height. At a maximal point, the height formula
for every lower row gives precisely every Ford inequality defining F;
a translation moves its horizontal coordinate into P. Thus every orbit
meets F. If two interior points of F are related by a matrix with
`c!=0`, the strict Ford inequality at the source would strictly reduce
height, contradicting maximality at the target. A matrix with `c=0`
is a translation, and P has disjoint translated interiors. Thus the
translates of F cover H³ with disjoint interiors. Boundary multiplicity
is represented by the face maps below. The source of the classical
Bianchi-Humbert/Swan construction and its strict fundamental-polyhedron
formulation is [Rahm et al., §2](https://arxiv.org/html/1207.7133v1).
The height/maximality argument above also supplies the needed coverage
and interior-disjointness proof for this concrete P.

## Exhaustive finite cutoff

The exact dyadic cover establishes, for every `z in P`,

    max_(322 recorded witnesses) q_(c,l)(z) >= 1/24.

On each square a single recorded `q` exceeds `1/24` at all four
corners. Its Hessian is negative definite, so its minimum on the square
is attained at a corner. The independent verifier proves the leaf squares
form a disjoint quadtree partition of P, tests primitivity by the gcd of
integer minors of `(c,c tau,l,l tau)`, and checks the four exact values.

Every unrecorded denominator of norm at least 25 has sphere radius
squared at most `1/25<1/24`, so it cannot touch F. Norm 24 does not
occur: `b=0` would require `a²=24`; `|b|>=2` implies norm at least 43;
and `|b|=1` would require `a²±a+11=24`, whose discriminant 53 is
not a square. Thus enumerate every `0<N(c)<=23`.

At any point of P, `|z|²<=13/4`. If a sphere with `N(c)<=23`
reaches the lower region, `|cz-l|<=1`, whence

    |l| <= sqrt(23*13/4)+1 < 10, so N(l)<100.

These exact bounds give 1,011 primitive `(c,l)` rows after identifying
simultaneous signs. Integer ranges in both producer and verifier contain
every such row; no Euclidean nearest-point division or Gaussian
generator assumption occurs. The largest omitted radius is therefore
strictly smaller than the certified floor. This is a concrete finite
termination proof, stronger than merely invoking the general finiteness
theorem in [Swan's criterion as presented by Rahm et al., §6](https://arxiv.org/html/1207.7133v1).

## Faces and pairings

In `(a,b)` coordinates every `q` equals `-N(z)` plus an affine
function. The region on which a candidate reaches the upper envelope
is therefore a rational convex polygon: intersect P with its affine
comparison halfspaces against all 1,010 other rows. Exact clipping
leaves 19 positive-area polygons with total area exactly 1. Their
individual vertices, areas and centers are in `geometry_result.json`.
The minimum among all 38 distinct vertices is `y²=2/43`.
Each polygon's `q` is concave, so that is the global minimum of the
finite floor, not just a sampled value.

For a lower face indexed by `(c,l)`, the ledger gives a matrix
`g=[[a,b],[c,-l]]` with exact determinant 1. The bounded search for
`a,b` is used only to produce the ledger; the verifier reconstructs the
rows, checks every determinant, and checks the image of every face.
On its hemisphere, the projected action is the affine isometry

    z -> a/c - conjugate(c*z-l)/c,     y -> y.

The verifier clips each image polygon into translated copies of P,
checks area conservation and tests at every vertex that the image is on
the reconstructed upper envelope. Since all comparisons are affine,
vertex inequalities prove them throughout each image piece. The
inverse pairing follows exactly from determinant 1. Some faces split
across the chosen vertical cell boundary; these pieces are checked
separately. The four vertical sides pair by `T_1` and `T_tau` and their
inverses; they commute because O is an additive lattice. This identifies
all geometric boundary faces, including those introduced by clipping.

## Stabilizers and cusp

For `P=(z,y)` fixed by `g=[[a,b],[c,d]]`, the height equation bounds
`N(c)<=1/y²<=43/2`, hence `N(c)<=21`; the equations for g and its
inverse bound `|a|,|d|<=sqrt(21*13/4)+1<10` in this cell.
For each exact vertex `(z,y²)`, the verifier independently enumerates
all such `c,d`, imposes the exact fixed-height equation, obtains
`a=c*z+conjugate(c*z+d)`, and solves `b=(a*d-1)/c` by integral
divisibility. It compares every matrix with the ledger and checks
inverse and multiplication closure modulo ±I. The 38 local vertex
orders in PSL are: twelve of order 1, six of order 2, sixteen of order
3, two of order 6, and two of order 12. An open face has trivial
pointwise stabilizer; for an edge segment, its pointwise stabilizer is
the intersection of the exact groups at its two endpoints. The ledger
records all 56 edges: 34 have two incident lower faces, while 22 lie on
vertical boundaries and have one. The verifier checks each incidence,
boundary classification and group intersection. These rules
also apply after subdivision for a rigid cell complex. They do not
claim a classification of setwise symmetries of the un-subdivided
polygonal faces.

At infinity the units ±1 give only translations in PSL. The cusp
lattice is exactly `Z+Z tau` of covolume `sqrt(43)/2`, with no rotational
quotient. For `c!=0`, a point above height 1 has image height at most
`1/(N(c)y)<1`; hence the horoball above `Y=2` is precisely invariant.
For every `z in P`, the finite lower floor is below 1, so the region
above 2 is exactly `P×(2,infinity)`. The dual frequencies and metric
weights remain as frozen in the design.

## Scope of the result

The geometry gate is closed for d=43. The next required theorem input is
a compact-core lower-bound argument and a validated mesh for this
19-patch floor and its face identifications. No spectrum statement is
claimed for d=43,67,163. Neither this polyhedron nor its constants
transfer to the other two fields.
