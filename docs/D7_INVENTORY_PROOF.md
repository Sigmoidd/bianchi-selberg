# Complete elliptic inventory and spectral certificate for d=7

Set s=√−7, τ=(1+s)/2, K=Q(s), O=Z[τ], Γ=PSL₂(O).
Coordinates (a,b) mean a+bτ; τ²=τ−2. This is a self-contained
field-specific proof, with general inputs Minkowski's ideal-class bound,
local valuation theory, and the standard trace formula
[Friedman, arXiv:math/0612807v1](https://arxiv.org/abs/math/0612807v1),
Theorem 4.1.1. No conjugacy-count table or unproved search cutoff is used.
`python -m groups.d7_inventory` replays the finite identities. The arguments
below supply completeness; the replay is not a formal proof assistant.

## Inventory theorem

There are exactly two elliptic **element** classes, both non-cuspidal.
The trace formula counts element classes, not cyclic subgroup classes.

| Class | Representative R | PSL order | Full maximal finite centralizer size | Primitive N(T₀) |
|---|---|---|---|---|
| A₂ | (0,−1;1,0) | 2 | 2 | 8+3√7 |
| A₃ | (0,−1;1,1) | 3 | 3 | (23+5√21)/2 |

Here matrix notation lists its two rows. Inverse order-3 elements are
conjugate in Γ. The involution has no endpoint flip in its centralizer.
Thus

\[
 C_{\rm ell}=\frac{\log(8+3\sqrt7)}8
       +\frac{\log((23+5\sqrt{21})/2)}9.
\]

These differences from d=2 are proved, rather than inferred from its data.

## D7-L1: base ring, cusps and all possible lattices

The norm is a²+ab+2b²=(a+b/2)²+7b²/4. To approximate a complex number
by O, first round its imaginary coefficient b and then its real coordinate.
The residual squared modulus is at most 7/16+1/4=11/16<1. Hence O is
norm-Euclidean, a PID, and O×={±1}. Bezout applied to a primitive pair
(p,q)∈O² gives an SL₂ matrix mapping infinity to p/q. Thus there is one
cusp, rotation index GG=1, and no cuspidal elliptic sector.

An elliptic SL lift has real trace strictly between −2 and 2. Since
O∩R=Z, its trace is 0 or ±1; negate to choose 0 or 1. The polynomials
X²+1 and X²−X+1 are irreducible over K. Let r=i or α=ζ₆, respectively.
As in D2-L1, an O² lattice carrying this action is a rank-one lattice
I in L=K(r), stable under A=O[r]. Two actions are GL₂(O)-conjugate iff
the lattices differ by multiplication by L×: a K-intertwiner is L-linear
on a one-dimensional L-space. Conversely every such lattice is O-free
of rank two because O is a PID. The multiplier ring is {x:xI⊂I}, with
multiplication determinant N_{L/K}(x).

## D7-L2: maximal orders and exhaustive ideal classes

The orders A₂=O[i] and A₃=O[α] have Z-basis (1,τ,r,τr).
Their trace-pairing determinants are 784=28² and 441=21².
Both are maximal. Indeed the quadratic discriminants −7 and −4
(respectively −7 and −3) are coprime. At every rational prime one
quadratic integer ring is unramified (possibly split). Its tensor with
the other integer ring is a product of unramified extensions of DVRs,
and is integrally closed. The tensor order is therefore maximal locally
everywhere, hence globally. The tensor spans the biquadratic field, so
this proves maximality of the displayed rank-four order.

For a quartic CM field Minkowski's bound is 3√|disc|/(2π²).
For A₂ it is 42/π²<5; only ideal norms 1,2,3,4 need consideration.
At 2, X²−X+2 has the two simple roots 0,1 mod 2, so the base
completions are Q₂. In each, i gives a totally ramified quadratic
extension: (Y+1)²+1=Y²+2Y+2 is Eisenstein at 2. Consequently there
are exactly two norm-2 primes and **no norm-4 prime**. Both norm-2
primes are principal: τ+i and 1−τ+i have absolute norm 2, and vanish
at different residue pairs (τ,i)=(1,1) and (0,1). At 3, X²−X+2 has
no F₃ root, so there is no norm-3 prime. Any ideal of norm 4 is a
product of the principal norm-2 primes. Every possible Minkowski
representative is principal, proving h(A₂)=1.

For A₃ the bound is 63/(2π²)<4. At 2, X²−X+1 has no F₂ root;
at 3, X²−X+2 has no F₃ root. Thus there is no norm-2 or norm-3
prime, and h(A₃)=1. Every A-stable lattice is a fractional ideal of
its maximal order. These class-number arguments prove exactly one
GL lattice type for each trace; there are no nonmaximal-order types.

## D7-L3: complete units from proved finite boxes

Write x=u+vr, with u=a+bτ, v=c+dτ. Relative conjugation σ fixes K.
A unit is exactly an element of relative norm ±1: necessity follows
from O×, and sufficiency from its integral inverse ±σx.

| Order | Relative norm | Expanding generator ε | Norm of ε | |ε|² | Torsion root |
|---|---|---|---|---|---|
| A₂ | u²+v² | (1+τ)+(2−τ)i | +1 | 8+3√7 | i |
| A₃ | u²+uv+v² | τ+α−1 | −1 | (5+√21)/2 | α |

In the chosen embedding these generators equal
(3+√7)(1+i)/2 and i(√7+√3)/2. Their relative norms and integral
inverses are checked exactly. Put E=|ε|>1. Multiplying any unit by a
power of ε reduces it to 1≤|x|<E, while |σx|=|x|⁻¹. Solving
v=(x−σx)/(r−σr) and u=(rσx−σr x)/(r−σr) bounds both squared
coefficient moduli strictly by (E+E⁻¹)²/|r−σr|². This is 9/2 for
A₂ and 7/3 for A₃. Since |a+bτ|²=(a+b/2)²+7b²/4, these bounds
force |b|,|d|≤1 and |a|,|c|≤2. This is a **proved exhaustive box**.

The exact absolute squares used to filter the box are

\[
 |u+vi|^2=N(u)+N(v)+(bc-ad)\sqrt7,
\]
\[
 |u+v\alpha|^2=
 \frac{2(N(u)+N(v))+2ac+ad+bc+4bd+(bc-ad)\sqrt{21}}2.
\]

Integer comparisons of radicals suffice; no float decides membership.
Among 12 norm-equation candidates for A₂, its four reduced units are
exactly the powers of i. Among 18 candidates for A₃, its six reduced
units are exactly the powers of α. Thus the entire unit groups are
μ₄ ε₂^Z and μ₆ ε₃^Z. Their determinant images are {1} and {±1}.
In particular **A₂ has no norm-minus-one unit**.

## D7-L4: SL splitting and PSL element multiplicities

A GL orbit splits into SL orbits indexed by O× modulo the determinant
image of its unit centralizer. A₂ therefore gives two SL lift classes.
D=diag(1,−1), of determinant −1, sends R₂ to −R₂. The two SL classes
are these opposite lifts and become the same element in PSL. There is
exactly one PSL involution class.

A₃ gives one SL class, since its unit determinant image is {±1}.
Explicitly J=(1,1;0,−1) represents relative conjugation and sends R₃
to R₃⁻¹. Multiplication U by ε₃ has determinant −1 and centralizes R₃.
Hence K=JU has determinant 1, square −I, and K R₃=R₃⁻¹ K. The
inverse elements merge already in SL. Passing to PSL introduces no
extra trace-one identification: conjugation cannot change trace 1 to −1.
This proves precisely the two element classes in the table.

Their eigenvalues are outside K, so neither fixes a K-rational cusp:
a K eigenline for a K matrix would give a K eigenvalue. This also follows
from the single cusp's trivial rotation quotient in L1.

## D7-L5: primitive translations and full centralizers

The SL centralizer is the norm-one unit group. Its smallest expanding
exponent is one for A₂ and two for A₃. Multiplication by torsion does
not affect eigenvalue modulus. Thus the primitive translations are

\[
 T_2=(1+\tau,-2+\tau;2-\tau,1+\tau),
 \qquad N(T_2)=8+3\sqrt7,
\]
\[
 T_3=(-2-\tau,1-2\tau;-1+2\tau,-3+\tau),
 \qquad N(T_3)=(23+5\sqrt{21})/2.
\]

The latter is multiplication by ε₃² and has trace −5. Each determinant
is one and each commutes with its class representative.

For R₂, a PSL centralizer lift could satisfy XR₂=±R₂X. All
anti-commuting maps have the form relative conjugation (determinant −1)
times a commuting unit. L3 proves every such unit has determinant +1,
so none is in SL. There is no endpoint flip. The centralizer's finite
subgroup is μ₄/{±1}, order 2; every remaining unit translates the axis.
For R₃, trace 1 excludes the minus sign, so its finite centralizer is
μ₆/{±1}, order 3. K from L4 is a **normalizer**, not a centralizer;
it must not double the orbital denominator.

The cylinder integral in D2-L7 applies unchanged: an endpoint-preserving
centralizer with rotation size q and primitive length ℓ contributes
g(0)ℓ/(4q sin²(π/m)); a full-centralizer endpoint flip would divide
this by a further factor two. Here both endpoint indices are one,
ℓ=log N(T₀), and q=2 or 3. This gives exactly the displayed coefficient.
D2-L7 includes the derivation from the spherical/Fourier transforms,
so the normalization uses an integral, not a patched legacy value.

## D7-gap: analytic enclosure and parameter comparison

Use the shared level-one backend: one cusp, GG=1, no cuspidal elliptics,
volume from ζ_K(2), systole from trace τ. The complete systole proof and
witness are in [SYSTOLES.md](SYSTOLES.md); the scattering reductions and
rigorous two-sided tails are in [ANALYTIC_DERIVATIONS.md](ANALYTIC_DERIVATIONS.md).
The exporter replays L1–L5, binds every record to this exact full group,
checks support inside the systole in Arb, and requires the upper B bound <1.

The selected sinc⁴ test uses k=2, R=256, support fraction 1: its exact
binary δ is rounded down and **still satisfies** 4δ≤ℓ₀ in Arb.
The frozen ball and its outward endpoints are in
[`d7-k2.json`](../certificates/d7-k2.json). The nonnegative real-axis
h(r)=sinc⁴(δr), with h(iσ)>1 for 0<σ<1, then excludes every discrete
Laplace eigenvalue in (0,1), after subtracting the constant eigenfunction.
The continuous spectrum starts at 1. This proves λ₁≥1, not a stronger gap.

[`d7-search.json`](../certificates/d7-search.json) records all 55 evaluations:
k=2,…,8 at seven support fractions and R=40, a five-point k=2
refinement at R=128, and the final k=2, fraction=1 evaluation at R=256.
The final point minimizes the **certified upper bound on this finite set**.
The larger cutoff sharpens tail errors; it is not a different test function.
This does not prove a continuous/global optimum or optimize over all
admissible functions. Both the manifest and replay say so explicitly.

```sh
python -m groups.d7_inventory
python examples/d7_certificate.py --output /tmp/d7-k2.json
python examples/d7_search.py --output-dir /tmp/d7-search
python -m unittest discover -s tests -v
```
