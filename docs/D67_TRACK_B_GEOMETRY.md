# d=67 full-group geometry certificate

Status: **geometry gate proved**, for the exact full group
`PSL_2(Z[(1+sqrt(-67))/2])` and the translation cell fixed in
`D67_TRACK_B_DESIGN.md`. This certifies the fundamental region, its
lower and vertical faces, explicit face maps, and pointwise orbifold
stabilizers. It does **not** certify the spectral exclusion target;
compact-core FEM theory and all-window positivity remain open.

Reproduce from the repository root:

```sh
python artifacts/track_b/d67/geometry_certify.py
python artifacts/track_b/d67/verify_geometry.py
python -O artifacts/track_b/d67/verify_geometry.py
python artifacts/track_b/d67/geometry_negative_checks.py
```

`cover_witnesses.json` has 712 exact rational squares and their
unimodular hemisphere witnesses. `geometry_result.json` has the
complete 37-face lower boundary ledger (every exact matrix and polygon),
102 edges with pointwise stabilizers, and the 66 vertex stabilizer groups
(all matrices modulo ±I).
The verifier reconstructs all 1,575 candidate rows, clips the polygons
against every competing hemisphere, checks image pieces of each paired
face against all rows, and independently enumerates vertex stabilizers.
All decisions are rational and remain active under `python -O`.

## Fundamental-domain argument

Write `tau=(1+sqrt(-67))/2` and `N(a+b tau)=a²+ab+17b²`.
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

    max_(712 recorded witnesses) q_(c,l)(z) >= 1/34.

On each square a single recorded `q` exceeds `1/34` at all four
corners. Its Hessian is negative definite, so its minimum on the square
is attained at a corner. The independent verifier proves the leaf squares
form a disjoint quadtree partition of P, tests primitivity by the gcd of
integer minors of `(c,c tau,l,l tau)`, and checks the four exact values.

Every denominator of norm at least 36 has sphere radius squared at
most `1/36<1/34`, so it cannot touch F. The completed-square norm
`(a+b/2)²+67b²/4` bounds `|b|<=1` for norms below 36; exact enumeration
gives `1,4,9,16,17,19,23,25,29`, with no norms 30 through 35.
Thus enumerate every `0<N(c)<=29`.

At any point of P, `|z|²<=19/4`. If such a sphere reaches the lower
region, `|cz-l|<=1`, whence

    |l| <= sqrt(29*19/4)+1 < 13, so N(l)<169.

The integer enumeration ranges are `-6<=a<=6,-1<=b<=1` for c and
`-14<=a<=14,-3<=b<=3` for l. The completed-square identity proves
these contain all permitted elements.

These exact bounds give 1,575 primitive `(c,l)` rows after identifying
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
comparison halfspaces against all 1,574 other rows. Exact clipping
leaves 37 positive-area polygons with total area exactly 1. Their
individual vertices, areas and centers are in `geometry_result.json`.
The minimum among all 66 distinct vertices is `y²=2/67`.
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
`N(c)<=1/y²<=67/2`, hence `N(c)<=33`; the equations for g and its
inverse bound `|a|,|d|<=sqrt(33*19/4)+1<14` in this cell, so their
norms are strictly below 196. The c range `[-7,7]×[-1,1]` and d range
`[-16,16]×[-3,3]` are exhaustive by the completed-square norm.
For each exact vertex `(z,y²)`, the verifier independently enumerates
all such `c,d`, imposes the exact fixed-height equation, obtains
`a=c*z+conjugate(c*z+d)`, and solves `b=(a*d-1)/c` by integral
divisibility. It compares every matrix with the ledger and checks
inverse and multiplication closure modulo ±I. The 66 local vertex
orders in PSL are: forty of order 1, six of order 2, sixteen of order
3, two of order 6, and two of order 12. An open face has trivial
pointwise stabilizer; for an edge segment, its pointwise stabilizer is
the intersection of the exact groups at its two endpoints. The ledger
records all 102 edges: 76 have two incident lower faces, while 26 lie on
vertical boundaries and have one. The verifier checks each incidence,
boundary classification and group intersection. These rules
also apply after subdivision for a rigid cell complex. They do not
claim a classification of setwise symmetries of the un-subdivided
polygonal faces.

At infinity the units ±1 give only translations in PSL. The cusp
lattice is exactly `Z+Z tau` of covolume `sqrt(67)/2`, with no rotational
quotient. For `c!=0`, a point above height 1 has image height at most
`1/(N(c)y)<1`; hence the horoball above `Y=2` is precisely invariant.
For every `z in P`, the finite lower floor is at most 1, so the region
above 2 is exactly `P×(2,infinity)`. The dual frequencies and metric
weights remain as frozen in the design.

## Scope of the result

The geometry gate is closed for d=67. The next required theorem input is
a compact-core lower-bound argument and a validated mesh for this
37-patch floor and its face identifications. No spectrum statement is
claimed for d=43,67,163. The relaxed compact-core criterion is now
rigorously obstructed for d=67; see `D67_TRACK_B_EXCLUSION_STATUS.md`. Neither this polyhedron nor its constants
transfer to the other two fields.
