# d=67 exclusion status and certified obstruction

**The spectral target is OPEN.** The geometry gate is closed by
`D67_TRACK_B_GEOMETRY.md`. No mesh, discretization-error ledger, or
all-window matrix-positivity certificate for the glued d=67 core exists.
The frozen initial design remains unchanged; this note records its next
dependency and a proved obstruction to one proposed relaxation.

## Compact-core reduction

Let `K={(z,y):z=a+b tau in P, sqrt(H(z))<=y<=2}` with the exact Ford
floor H. Its volume measure is `dx1 dx2 dy/y³`, energy measure is
`|grad_euc v|² dx1 dx2 dy/y`, and cusp area is `A=sqrt(67)/2`.
The relevant form domain consists of H¹ functions whose boundary traces
match under all face identifications. Pointwise stabilizers impose the
orbifold interpretation; no artificial Dirichlet boundary is introduced.

For an L² eigenfunction at `lambda=1-s²`, `0<s<1`, its cusp constant
mode is `c(y/2)^(1-s)`. The growing solution is not L². Integration of
this mode gives cusp mass `A c²/(2s*4)` and energy `(1-s)²` times this
mass. With `t=integral_P v(z,2) dx1 dx2=A c`, its contribution to
`Q-lambda M` is `-(1-s)t²/(4A)`. Orthogonality to the constant function
gives `L_s=integral_K v dV+t/(4(1+s))=0`.

The dual squared frequencies are `m²+(2n-m)²/67`, with least nonzero
value `4/67`. At y>=2 their horizontal energy/mass coefficient is at
least `64 pi²/67>1`. Thus the nonzero cusp modes contribute nonnegative
energy minus lambda times mass. The global eigenfunction identity then
gives `B_s(v)<=0`, where

    B_s(v)=Q_K(v)-(1-s²)M_K(v)-(1-s)t(v)²/(4A).

Consequently strict positivity of B_s on nonzero matching-trace v with
L_s=0, for every s in (0,1), is sufficient for exclusion. Unique
continuation rules out v identically zero on the core. The standard
finite-volume spectral decomposition places essential spectrum at 1;
this reduction treats both residual and cuspidal L² eigenfunctions.
These field-independent analytic inputs are the same dependencies as
`independent_exclusion/DESIGN.md`; its Gaussian numerical data do not
supply the required d=67 positivity.

## Why the unconstrained H¹(K) relaxation fails

The Gaussian implementation tests a sufficient inequality on all H¹(K),
relaxing the face identifications. This larger-domain criterion cannot
succeed for the present d=67 K. This is proved by an explicit trial,
not inferred from a failed mesh or approximate eigenvalue search.

Take `v(z,y)=b`. The candidate set is exactly invariant under
`(c,l)->(c,-l)`, hence `H(-z)=H(z)`. Both core mean and top integral of
v vanish by central symmetry. Therefore L_s=t=0 for every s.
The inverse horizontal metric gives `|grad_euc b|²=4/67`, so

    Q_K(v)/A=(4/67) integral_P log(2/sqrt(H)) da db,
    M_K(v)/A=integral_P b² (1/(2H)-1/8) da db.

The floor cover gives H>=1/34. Since
`exp(5)>sum_(k=0)^11 5^k/k!=58918733/399168>136`, we have
`log(2 sqrt(34))<5/2`, proving `Q/A<10/67`.

The independent obstruction script partitions P into 16×16 exact
squares. For every candidate q on each square it obtains the maximum
by the stationary interior point (if present) and the clamped stationary
points on its four edges. Concavity makes this exhaustive. The largest
of these maxima bounds H from above on the square. Multiplying
`1/(2 H_max)-1/8` by the exact integral of b² gives a lower mass bound.
The exact sum exceeds `6/25`. At `s=1/2` it follows that

    B_(1/2)(v)/A < 10/67 - (3/4)(6/25) = -103/3350 < 0.

This v fails the quotient boundary identifications: in particular its
traces differ by 1 under translation by tau. It is not a quotient trial
and gives no upper bound for the quotient's first eigenvalue.
The result obstructs only the all-H¹(K) sufficient criterion. It does
not prove impossibility of a glued-core criterion or of trace tests.
The initial 8×8 bound was insufficient; the 16×16 exact bound passes.

Replay from repository root:

```sh
python artifacts/track_b/d67/relaxation_obstruction.py
python -O artifacts/track_b/d67/relaxation_obstruction.py
```

## First remaining input

A guaranteed finite-element lower-bound theorem and validated mesh that
retain the d=67 face identifications are required. The existing Gaussian
floor inclusion proof in `lower_bound_theory.md` uses the false general
triangle variance bound `d²/4`; an equilateral triangle has variance
`d²/3` at its centroid. Keep that historical text and its certificates
unchanged. Any new inclusion lemma must prove a valid bound (the general
variance bound `d²/3` suffices for the corresponding Taylor sag estimate).
The curved paired faces also need compatible trace/interpolation and
sliver estimates; imposing identifications only on a numerical matrix
without proving their continuum correspondence is insufficient.

After these inputs, an independently reconstructed, interval-certified
pencil must cover the whole parameter interval, including both endpoints
as limits. No such field-specific matrices or certified windows are
present. The theorem claim stops here.
