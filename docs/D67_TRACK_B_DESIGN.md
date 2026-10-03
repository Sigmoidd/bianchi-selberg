# d=67 Track B design — frozen initial geometry scope

Status: the full level-one spectral exclusion target is OPEN. This
branch begins d=67 independently of the d=43 geometry artifact.

Set `tau=(1+sqrt(-67))/2`, `O=Z[tau]`, `tau²=tau-17`, and
`Gamma=PSL_2(O)`. The field discriminant is -67. The norm
`N(a+b tau)=a²+ab+17b²=(a+b/2)²+67b²/4` gives only units ±1.
The reduced positive form list at discriminant -67 is `(1,1,17)`;
the ideal-class/form correspondence gives class number one and hence
one cusp. No norm-Euclidean division is used.

Use the curvature -1 metric `(dx1²+dx2²+dy²)/y²`, measure
`dx1 dx2 dy/y³`, and positive Laplacian
`-y²(partial_x1²+partial_x2²+partial_y²)+y partial_y`.
The exceptional target is all discrete eigenvalues `lambda=1-s²`,
`0<s<1`, including cuspidal and residual. The essential threshold is 1.
Set `P={a+b tau:-1/2<=a,b<=1/2}`, `A=sqrt(67)/2`, and `Y=2`.
The infinity stabilizer in PSL consists of translations by O with
no rotational quotient. The horoball `y>1` is precisely invariant:
for any matrix with `c!=0`, its image height is at most
`1/(N(c)y)<1` there.

The Euclidean dual under `Re(mu conjugate(z))` has
`mu_(m,n)=(m,(2n-m)/sqrt(67))`, so the least nonzero squared
frequency is `4/67`. At Y=2 its cusp energy coefficient satisfies
`4 pi²*(4/67)*Y²>64*9/67>1` using `pi>3`.

The chosen exclusion method remains the compact-core criterion, with

    L_s(v)=integral_K v dV + t(v)/((1+s)Y²),
    B_s(v)=Q_K(v)-(1-s²)M_K(v)-(1-s)t(v)²/(A Y²),
    t(v)=integral_P v(z,Y) dx1 dx2.

Prove `B_s>0` for every nonzero `v` in the correct compact-core form
domain with `L_s(v)=0`, uniformly for `0<s<1`, using guaranteed
nonconforming lower bounds. The zero cusp mode is exactly
`c(y/Y)^(1-s)` and gives the displayed boundary and mean terms.
This is the field-independent derivation from the d=43 design with
the new cusp area and dual lattice substituted and checked.

For geometry define the Ford height floor as the upper envelope of
`q_(c,l)(z)=1/N(c)-|z-l/c|²` for primitive pairs `(c,l)`.
An exact cover certificate now proves it is at least `1/34` on P.
The possible norms below 36 are: `1,4,9,16,17,19,23,25,29`;
there are no norms 30 through 35. Thus every omitted hemisphere
with `N(c)>=36` has radius squared at most `1/36<1/34` and cannot
cut the floor. The finite candidate set is bounded by `N(c)<=29`;
if such a hemisphere meets the region then
`|l|<=sqrt(29*19/4)+1<13`, hence `N(l)<169`.
This is a field-specific finite termination bound, not yet a
complete face/edge arrangement.

An initial exact-rational face probe enumerates 1,575 primitive candidate
rows inside these bounds. It reports 37 positive-area projected patches,
area sum 1 and minimum enumerated vertex height squared `2/67`.
`artifacts/track_b/d67/face_probe.py` reproduces this diagnostic.
Independent reconstruction of its polygons, pairing maps and cell
stabilizers has not yet been done, so the geometry gate remains open.

The next geometry work is to enumerate all primitive rows inside
those bounds, exact-clip every lower face against its competitors,
produce determinant-one pairing matrices and vertical translations,
and independently enumerate pointwise stabilizers with a separate
verifier. Only then adapt the compact-core mesh and positivity
certificate. No d=43 face, matrix, height or mesh constant transfers.
