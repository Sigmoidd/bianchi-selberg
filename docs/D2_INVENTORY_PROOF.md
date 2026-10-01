# Complete elliptic inventory for PSL₂(Z[√−2])

This note proves the field-specific arithmetic completeness used by the d=2
certificate. It uses no published Bianchi conjugacy-count table and no
bounded search for matrix conjugators. `python -m groups.d2_inventory`
replays the finite arithmetic steps. The local arguments, ideal-class bound,
lattice correspondence, and orbital integral below are mathematical proofs;
the replay is not a proof-assistant formalization.

Write s=√−2, K=Q(s), O=Z[s], and Γ=PSL₂(O). All displayed matrices have
entries in O. Arithmetic uses the integral basis (1,s), so an entry (a,b)
in the JSON means a+b s. The trace formula is the standard external input:
[Friedman, arXiv:math/0612807v1](https://arxiv.org/abs/math/0612807v1),
Theorem 4.1.1, printed pp.41–42. Minkowski's ideal-class bound and elementary
local valuation theory are the general number-theory inputs. Everything
specific to these orders is established here.

## The inventory and coefficient

There are four Γ-conjugacy classes of elliptic **elements**:

| Label | Representative R | PSL order m | Maximal finite centralizer size | Minimal N(T₀) |
|---|---|---|---|---|
| A2 | (0,−1; 1,0) | 2 | 4 | 17+12√2 |
| S2 | (1,s; s,−1) | 2 | 4 | 3+2√2 |
| α | (0,−1; 1,1) | 3 | 3 | 5+2√6 |
| α⁻¹ | (1,1; −1,0) | 3 | 3 | 5+2√6 |

The matrix notation `(a,b; c,d)` means rows (a,b) and (c,d).
All four are non-cuspidal. The two order-3 classes are inverse elements of
one cyclic subgroup; the trace formula sums element classes, so both count.
The non-cuspidal coefficient is

\[
 C_{\rm ell}=
 \frac{\log(17+12\sqrt2)+\log(3+2\sqrt2)}{16}
 +\frac{2\log(5+2\sqrt6)}9
 =\frac3{16}\log(3+2\sqrt2)+\frac29\log(5+2\sqrt6).
\]

The denominator 16 for each involution uses the full PSL centralizer.
Section "Orbital normalization" derives this factor directly. Existing
Eisenstein input records remain frozen; this proof does not re-certify them.

## Reduction to quadratic lattices

O is norm-Euclidean: rounding the real coefficient and s-coefficient of
any element of K gives error norm at most 1/4+2/4=3/4<1. Therefore O is a PID
and its only units are ±1. Each cusp in P¹(K) can be written p/q with
coprime p,q∈O. Bezout supplies an SL₂(O) matrix with first column (p,q),
mapping infinity to that cusp, so there is exactly one cusp. Its stabilizer
has diagonal entries in O×={±1} and is unipotent after passing to PSL₂.
The cusp rotation quotient has size |O×|/2=1 and has no elliptic
sector. The same class number is checked by reduced forms in `SYSTOLES.md`.

An elliptic SL₂(O) lift has real trace in (−2,2). Since O∩R=Z, its trace is
0 or ±1. Negating the lift lets us take 0 or 1. Its polynomial is respectively
X²+1 or X²−X+1; both are irreducible over K. Thus an O² lattice with this
action becomes a rank-one lattice I in the quadratic extension L/K stable
under A=O[root]. Conversely such a lattice is O-free of rank two (O is a
PID), and an O-basis gives a matrix of the required trace and determinant.
Two matrices are GL₂(O)-conjugate exactly when their lattices differ by a
factor in L×. This follows by extending an intertwiner to K: it is an
L-linear map between one-dimensional L spaces, hence multiplication by
that factor. The integral commuting ring is

\[
 \{x\in L:xI\subset I\},
\]

and the determinant of its multiplication action is N_{L/K}(x). This
correspondence includes lattices that are not invertible over A.

## Trace zero: two lattice types, including the nonmaximal order

Put z=ζ₈=(1+i)/√2. Then s=z+z³, z²=s z+1, and i=1+s z. Define

\[
 S_2=O[z]=Z[z],\qquad A_2=O[i]\subset S_2.
\]

Their Z-bases are (1,s,z,sz) and (1,s,i,si). Direct trace-pairing
determinants give discriminants 256 and 1024. The inclusion has index two:
i=1+sz and si=s−2z.

**Maximality of S₂.** Its power basis equals the displayed basis, with
unimodular change of basis. The polynomial of z is X⁴+1. At 2 the shifted
polynomial (Y+1)⁴+1 is Eisenstein, so π=z−1 generates a totally ramified
extension of Q₂ of degree four and has valuation 1/4, with v(2)=1.
For x=Σ_{j=0}³ c_j π^j, c_j∈Q₂, the nonzero summands have distinct
valuation fractional parts. Thus v(x)≥0 forces v₂(c_j)≥0 for every j.
This proves the full local integer ring is Z₂[π]. No odd prime can divide
the index, since the discriminant is 2⁸. Hence S₂ is maximal.

**Class number of S₂.** Minkowski gives an integral ideal of norm at most

\[
 (4!/4^4)(4/\pi)^2\sqrt{256}=24/\pi^2<3
\]

in each ideal class. The only possible nontrivial norm is 2. The unique
prime above 2 is principal: the absolute norm of 1−z is 2. Thus S₂ has
class number one.

**All A₂-stable lattices.** Let f=s S₂. Its norm is 4 and it lies in A₂:
its basis generators are s, −2, sz=i−1, and −2z=si−s. Also f is an S₂
ideal. In the quotient,

\[
 S_2/f=\mathbb F_2[z]/((z+1)^2),\qquad A_2/f=\mathbb F_2.
\]

For any lattice I, the extended S₂-ideal S₂I is principal. Rescale I so
S₂I=S₂. Then f⊂I⊂S₂: indeed f S₂I⊂I because f S₂⊂A₂.
Consequently I/f is an F₂-linear subspace of the four-element quotient,
and its S₂-span must be the whole quotient. The full list of subspaces is
0, ⟨1⟩, ⟨z⟩, ⟨1+z⟩, and the whole quotient. Only ⟨1⟩, ⟨z⟩, and the
whole quotient have full span. The integral unit z exchanges the first
two lines, so they give the same lattice type. Their lifts are A₂ and
z A₂; the remaining lift is S₂. These are exactly two GL₂ lattice types.
Their multiplier rings are A₂ and S₂: for either order I containing 1,
xI⊂I implies x∈I, and conversely every element of I multiplies I into I.

Multiplication by i on the bases (1,i) and (1,z) gives the two listed
involutions. They also cannot merge after negating a lift: modulo the ideal
(s), R_A is non-scalar whereas R_S is scalar. Both properties are invariant
under conjugation, and the residue characteristic is two.

## Trace one: one GL₂ lattice type

Put α=ζ₆=(1+√−3)/2 and A₃=O[α], with polynomial α²−α+1=0.
The basis (1,s,α,sα) has discriminant 576. This order is maximal:
the maximal quadratic orders of discriminants −8 and −3 have coprime
discriminants. Locally at any rational prime at least one factor is
unramified (a product of unramified extensions is allowed). Tensoring
with it produces unramified extensions of the other factor's integer
rings, which are integrally closed DVRs. Hence the tensor order is
maximal at every prime.

Minkowski's bound is 36/π²<4. There is no norm-two prime, since
X²−X+1 has no F₂ root. At 3, reduction gives s=±1 and α=−1, so there are
exactly two norm-three primes. They are principal: β=s+α and
β̄=−s+α have absolute norm 3, and reducing them distinguishes the two
primes. Therefore A₃ has class number one. Every A₃-stable lattice is
an ideal of this maximal order, so up to rescaling there is exactly one
GL₂ lattice type. Multiplication by α gives R_α.

## Complete unit groups from finite, proved bounds

The following computation proves the unit groups needed for both the
SL₂ class splitting and the primitive translation lengths. It is not a
search with an unproved cutoff. Write x=u+v r, u=a+b s, v=c+d s, where
r is i, z, or α in A₂, S₂, or A₃. Relative conjugation σ fixes K.
A unit must have relative norm ±1, and this condition is sufficient
because σ(x) is integral in each order.

| Order | Relative norm | Expanding generator ε | N_{L/K}(ε) | Torsion generator |
|---|---|---|---|---|
| A₂ | u²+v² | 1−si = 1+√2 | −1 | i, order 4 |
| S₂ | u²+suv−v² | 1−s+2z = 1+√2 | −1 | z, order 8 |
| A₃ | u²+uv+v² | s+2α−1 = i(√2+√3) | +1 | α, order 6 |

Let E=|ε|>1. Multiplying any unit by an integer power of ε reduces it
to 1≤|x|<E. Since |N_{L/K}(x)|=1, |σx|=|x|⁻¹. Both r and σr have
absolute modulus one. Solving for u,v gives

\[
 v=\frac{x-\sigma x}{r-\sigma r},\qquad
 u=\frac{r\sigma x-\sigma r\,x}{r-\sigma r}.
\]

Thus both coefficient moduli are strictly less than
(E+E⁻¹)/|r−σr|. The squared bounds are exactly 2 for A₂ and 4 for
S₂ and A₃: |r−σr|² is respectively 4,2,3, whereas (E+E⁻¹)² is 8,8,12.
Consequently a²+2b² and c²+2d² are <2 or <4. Every coefficient lies in
{−1,0,1}; in A₂, b=d=0. This is an exhaustive finite box.

The checker filters this box by relative norm ±1 and the reduced modulus
condition, using exact integer comparisons of a+b√2 or a+b√6. The
absolute squares used in that comparison are

\[
\begin{array}{ll}
 A_2:& a^2+2b^2+c^2+2d^2+2(bc-ad)\sqrt2,\\
 S_2:& a^2+2b^2+c^2+2d^2+2(bc-ad)+(ac+2bd)\sqrt2,\\
 A_3:& a^2+2b^2+c^2+2d^2+ac+2bd+(bc-ad)\sqrt6.
\end{array}
\]

| Order | Norm-equation candidates in the box | Reduced units | Exact reduced list |
|---|---|---|---|
| A₂ | 4 | 4 | powers of i |
| S₂ | 16 | 8 | powers of z |
| A₃ | 10 | 6 | powers of α |

The polynomial equations and integral inverses of ε are checked too.
It follows that each full unit group is precisely the product of its
listed torsion group and ε^Z. In particular the relative determinant
images are {±1}, {±1}, and {1}, respectively. No norm-minus-one unit
exists in A₃.

## GL₂ to SL₂ to PSL₂: four element classes

The determinant map GL₂(O)→O×={±1} is onto. Within a GL₂ conjugacy
orbit, SL₂ orbits are the determinant cosets modulo the determinant
image of the stabilizer, which is the relative norm image of the
multiplier-order unit group. Thus each trace-zero GL₂ type gives one
SL₂ class, while the trace-one GL₂ type gives two SL₂ classes.

For trace one the relative conjugation matrix J=(1,1; 0,−1) has
determinant −1 and sends R_α to R_α⁻¹. If there were an SL₂ conjugator,
its difference from J would be a centralizer unit of determinant −1,
contrary to the unit calculation. Thus the two inverse classes are
distinct. Passing to PSL₂ cannot merge them: their chosen lifts have
trace one, so the minus sign in a possible PSL conjugacy would force
trace minus one. For trace zero, allowing a sign introduces no further
types and the two multiplier types remain distinct modulo (s).

There are no cuspidal elliptics: their eigenvalues lie in L\K, whereas
fixing a cusp in P¹(K) would give a K eigenline and thus a K eigenvalue.
Alternatively the single cusp stabilizer has trivial rotation quotient.
This closes completeness for all elliptic elements of Γ.

## Primitive translations and maximal finite centralizers

SL centralizer elements are exactly norm-one units in the corresponding
multiplier order. The unit description proves the following translation
matrices have minimal positive translation; no finite enumeration of
loxodromic matrices is needed:

\[
\begin{array}{ll}
 T_A=3I-2sR_A=(3,2s;-2s,3), &N(T_A)=17+12\sqrt2,\\
 T_S=(2,1+s;1+s,s), &N(T_S)=3+2\sqrt2,\\
 T_\alpha=sI+2R_\alpha-I=(-1+s,-2;2,1+s),
   &N(T_\alpha)=5+2\sqrt6,\\
 T_{\alpha^{-1}}=sI+2R_{\alpha^{-1}}-I=(1+s,2;-2,-1+s),
   &N(T_{\alpha^{-1}})=5+2\sqrt6.
\end{array}
\]

For A₂, torsion units have norm +1 and ε has norm −1, so the smallest
positive exponent is 2. For S₂, z has norm −1, so zε already has norm
+1 and exponent 1. For A₃, ε has norm +1, hence exponent 1 is minimal.
Multiplication by torsion preserves the absolute eigenvalue modulus.
The inverse-class witness has a contracting eigenvalue in the chosen
embedding; its other eigenvalue gives the same expanding N(T₀).

In PSL₂, centralizing R means XR=±RX on lifts. For trace-one lifts,
trace invariance excludes the minus sign. Their torsion centralizer is
μ₆/{±1}, of size three, and all other units translate the axis. Hence the
maximal finite centralizer size is three.

For each involution the norm-one torsion is μ₄/{±1}, of size two, but
there are also endpoint flips. Specifically R_A R_S=−R_S R_A and both
squares are −I. The four PSL classes of I,R_A,R_S,R_A R_S form a
Klein-four subgroup of each full centralizer. To prove maximality,
restrict the centralizer to its invariant geodesic axis. Its kernel is
the order-two rotation subgroup. A finite subgroup of the isometry
group of a line has order at most two. Therefore any finite subgroup
of the full centralizer has size at most 2·2=4, which the witness attains.
An endpoint-reversing lift is off-diagonal in an eigenbasis for R, so has
trace zero and square −I. It has zero translation and is finite; it
cannot give a smaller loxodromic norm than the orientation-preserving
unit calculation.

## Orbital normalization from the cylinder integral

Let C⁺ be the subgroup preserving both endpoints, q its rotation-kernel
size, and e=[C:C⁺] (one without flips and two with flips). Let
ℓ=log N(T₀) be the shortest positive axis translation. In cylindrical
coordinates (ρ,θ,t) about the axis,

\[
 ds^2=d\rho^2+\sinh^2\rho\,d\theta^2+\cosh^2\rho\,dt^2,
 \qquad dV=\sinh\rho\cosh\rho\,d\rho\,d\theta\,dt.
\]

For rotation angle 2ϕ=2πk/m, its displacement obeys
u=cosh d(P,RP)=1+2sin²ϕ sinh²ρ. A C⁺ fundamental cylinder has angular
width 2π/q and axis length ℓ (a screw twist identifies its ends but
preserves volume). The Selberg transform in the repository's Fourier
normalization satisfies

\[
 g(0)=2\pi\int_1^\infty k(u)\,du.
\]

Indeed the radial spherical transform is
h(r)=4π/r ∫₀^∞ k(cosh ρ)sin(rρ)sinh ρ dρ; integrating by parts shows
its inverse Fourier transform is g(t)=2π∫_{cosh t}^∞ k(u)du. Consequently

\[
 \int_{C^+\backslash\mathbb H^3} k(\cosh d(P,RP))\,dV
 =\frac{2\pi\ell}{q}\frac1{4\sin^2\phi}\int_1^\infty k(u)\,du
 =\frac{g(0)\ell}{4q\sin^2\phi}.
\]

The full centralizer quotient divides this by e. In these cases the
maximal finite subgroup has size q e, proved above, so the coefficient is
exactly g(0)log N(T₀)/(4|E(R)|sin²ϕ), as in Friedman's theorem. For an
involution q=2,e=2, not just q=2. For order three q=3,e=1 and
sin²ϕ=3/4 for both inverse element classes. This proves the normalization
used in the inventory table.

## Certificate and replay boundary

`GroupData.require_inventory()` invokes the verifier for the complete d=2
records before full assembly. Export invokes it again and checks the
systole/support and B<1. It compares every class, matrix, norm, finite
centralizer size, and multiplicity against the arithmetic witnesses;
merely changing a status string cannot admit another field. Fields d=7,
11,19 still have incomplete inventories and remain blocked.

The certificate combines this proof with the systole proof, analytic
reductions and rigorous two-sided quadrature tails documented in
`SYSTOLES.md` and `ANALYTIC_DERIVATIONS.md`. Positivity of
h(r)=sinc⁴(δr) on R and h(iσ)=(sinh(δσ)/(δσ))⁴>1 for 0<σ<1 then excludes
every discrete eigenvalue in (0,1) when the bound after subtracting the
constant eigenfunction is below one. The continuous spectrum begins at
one. No claim about a higher spectral gap is needed for M2.
