# Complete elliptic inventory for d=11

Let s=√−11, τ=(1+s)/2, K=Q(s), O=Z[τ], and Γ=PSL₂(O).
Coordinates (a,b) represent a+bτ; τ²=τ−3. The proof uses ordinary
Minkowski ideal-class bounds, local valuation theory, the rank-one lattice
correspondence described below, and the trace formula and orbital integral
already derived in [D2_INVENTORY_PROOF.md](D2_INVENTORY_PROOF.md).
The exact finite identities replay in `groups/d11_inventory.py`; the
completeness arguments are mathematical arguments, not a formalization.
No published elliptic class-count table or unproved search cutoff is used.

## Inventory theorem

There are exactly three noncuspidal elliptic **element** classes:

| Class | SL representative | PSL order | Full maximal finite centralizer size | Primitive N(T₀) |
|---|---|---|---|---|
| A₂ | (0,−1;1,0) | 2 | 4 | 199+60√11 |
| A₃ | (0,−1;1,1) | 3 | 3 | 23+4√33 |
| A₃ inverse | (1,1;−1,0) | 3 | 3 | 23+4√33 |

The involution has an endpoint flip. The two order-3 inverses do **not**
merge in SL or PSL. Thus the coefficient is

\[
 C_{ell}=\frac{\log(199+60\sqrt{11})}{16}
          +\frac{2\log(23+4\sqrt{33})}{9}
 =\frac{\log(10+3\sqrt{11})}{8}
          +\frac{2\log(23+4\sqrt{33})}{9}.
\]

These facts differ from d=7 and must not be copied from its inventory.

## D11-L1: base PID, cusp and lattice correspondence

The norm is a²+ab+3b²=(a+b/2)²+11b²/4. Rounding the imaginary
coefficient b and then the real coordinate leaves squared error at most
11/16+1/4=15/16<1. Consequently O is norm-Euclidean, hence a PID.
The norm equation gives O×={±1}. Bezout completes any primitive cusp
pair to an SL₂ matrix; there is one cusp and rotation index GG=1.

A noncentral elliptic lift has integral real trace 0 or ±1. Negating the
lift reduces to X²+1 or X²−X+1. These polynomials are irreducible over
K. Write r=i or α=ζ₆ and A=O[r]. An action on O² corresponds to an
A-stable O-lattice in the one-dimensional L=K(r) space. A K-intertwiner
is multiplication by L×, so GL₂(O) classes are lattices up to that
multiplication. The lattice is O-free because O is a PID. A unit's
multiplication determinant is N_{L/K}(x).

Neither root is in K, so no elliptic fixes a K-rational cusp; there are
no cuspidal elliptics. This is also consistent with the trivial cusp
rotation quotient.

## D11-L2: maximal orders and all ideal/lattice classes

For A₂=O[i] and A₃=O[α], the basis (1,τ,r,τr) has discriminant
1936=44² and 1089=33². These orders are maximal: the pairs of quadratic
integer discriminants (−11,−4) and (−11,−3) are coprime. At each rational
prime at least one factor is unramified (possibly split); its tensor with
the other integer ring is a product of unramified extensions of DVRs,
and integrally closed. The tensor spans L, proving maximality locally
and globally. Every A-stable lattice is therefore a fractional ideal
of this maximal order; there are no conductor or nonmaximal lattice types.

The quartic CM Minkowski bound is 3√|disc|/(2π²).

### A₂: bound 66/π²<7

Every ideal class has an integral representative of norm at most 6.
The following list treats every prime ideal norm that can divide such
an ideal; composite representatives then are products of principal primes.

- At 2, τ²−τ+3 is irreducible modulo 2, giving an unramified base
  completion with residue F₄. Adjoining i is totally ramified: i−1
  satisfies Y²+2Y+2, Eisenstein over that completion. There is exactly
  one prime of norm 4 and no norm-2 prime. The element 1+i has absolute
  norm 4 and generates that prime.
- At 3, τ has the two simple residues 0,1. On either F₃ factor, X²+1
  is irreducible, so the primes have norm 9; no norm-3 prime occurs.
- At 5, τ has simple residues 2,4 and i has residues 2,3. These give
  exactly four norm-5 primes. The elements τ+i, τ−i, 1−τ+i, and
  1−τ−i each have absolute norm 5 and vanish at exactly one distinct
  residue pair. They therefore generate all four primes.

No additional prime ideal norm ≤6 is possible. In particular norms 6
would require a norm-2 and norm-3 factor, neither of which exists. Thus
h(A₂)=1.

### A₃: bound 99/(2π²)<6

Every ideal class has a representative of norm at most 5.

- At 2, the base residue field is F₄ and α²−α+1 has its two distinct
  roots α=τ and α=τ+1. The relative extension is unramified and split.
  There are two norm-4 primes and no norm-2 prime. The elements τ−α
  and τ+α−1 have relative norm −2, absolute norm 4, and vanish on
  different residue factors. Both norm-4 primes are principal.
- At 3, the base splits with τ=0 or 1. In each base completion α is
  ramified: Y=α+1 satisfies Y²−3Y+3, Eisenstein at 3. There are exactly
  two norm-3 primes. The elements −2+τα and −2+(1−τ)α have absolute
  norm 3 and vanish respectively at τ=1 and τ=0, α=2; both are principal.
- At 5, the base splits but α²−α+1 has no F₅ root. Its primes have
  norm 25, so there is no norm-5 prime.

This exhausts all prime/composite norms ≤5. Thus h(A₃)=1. Both traces
have one GL lattice type, not an inferred count from h(K).

## D11-L3: complete unit groups with proved finite bounds

Write x=u+vr with u=a+bτ, v=c+dτ. Relative conjugation σ fixes K.
An element is a unit iff its relative norm is ±1: O×={±1}, and the
inverse then is ±σx.

| Order | Relative norm form | Expanding generator ε | N_{L/K}(ε) | |ε|² |
|---|---|---|---|---|
| A₂ | u²+v² | (1+τ)+(2−τ)i | −1 | 10+3√11 |
| A₃ | u²+uv+v² | (−3+2τ)+4α=s+2√−3 | +1 | 23+4√33 |

These displayed units and integral inverses are checked exactly. Set
E=|ε|>1. Multiplication by an integral power of ε reduces any unit to
1≤|x|<E. Since |σx|=|x|⁻¹, solving for u and v bounds both squared
coefficient moduli strictly by (E+E⁻¹)²/|r−σr|². The bounds are 11/2
for A₂ and 16 for A₃.

For A₂, |a+bτ|²<11/2 forces |b|≤1 and |a|≤2. For A₃,
|a+bτ|²<16 forces |b|≤2 and |a|≤4. The same applies to c,d.
Thus the complete coordinate boxes are respectively
[-2,2]×[-1,1]×[-2,2]×[-1,1] and
[-4,4]×[-2,2]×[-4,4]×[-2,2]. The strict inequalities handle the
boundary; no float decides unit membership.

The exact modulus filters are

\[
 |u+vi|²=N(u)+N(v)+(bc-ad)\sqrt{11},
\]
\[
 |u+v\alpha|²=
 \frac{2(N(u)+N(v))+2ac+ad+bc+6bd+(bc-ad)\sqrt{33}}{2}.
\]

Integer radical comparisons in the proved boxes give 12 norm-equation
candidates and exactly four reduced units for A₂; 10 candidates and exactly
six reduced units for A₃. The reduced units are precisely μ₄ and μ₆.
It follows that the **entire** unit groups are μ₄ ε₂^Z and μ₆ ε₃^Z.
Their determinant images are {±1} and {1}. In particular, A₃ has no
relative norm-minus-one unit anywhere outside the box either, since every
unit reduces to a torsion unit by powers of the norm-one generator.

## D11-L4: SL/PSL multiplicity and inverse classes

A GL orbit splits into SL orbits indexed by O× modulo the determinant
image of its unit centralizer. A₂ therefore gives one SL orbit, hence one
PSL involution class. Opposite involution lifts cannot create an extra class.

A₃ gives two SL orbits. Relative conjugation J=(1,1;0,−1) has
determinant −1 and sends R₃=(0,−1;1,1) to R₃⁻¹=(1,1;−1,0).
Any other intertwiner between them differs from J by a commuting unit.
All such units have determinant +1 by L3, so the determinant obstruction
cannot be removed. The two inverse element classes are distinct in SL.
Passing to PSL does not merge them: a trace-one lift cannot be conjugate
to a trace-minus-one lift. Thus these are exactly the two trace-one
PSL element classes, rather than a single cyclic subgroup contribution.

## D11-L5: primitive translations and full centralizers

For the involution, the SL centralizer consists of norm-one units; its
first expanding generator is ε₂². The primitive norm is therefore
(10+3√11)²=199+60√11. Rotation torsion is μ₄/{±1}, of order 2.

Let J₂=diag(1,−1) and U₂ be multiplication by ε₂. Both have determinant
−1, so X=J₂U₂ belongs to SL, has X²=−I, and anticommutes with R₂.
It is an endpoint flip that centralizes R₂ in PSL. The subgroup generated
by R₂ and X has four elements in PSL, checked by exact closure. It is
maximal finite: the orientation-preserving kernel has size 2 and the
axis-orientation quotient has size at most 2. The full finite centralizer
size is **4**, not 2. Any endpoint-reversing lift is conjugation times a
commuting unit; the norm-minus-one unit supplies exactly this possibility.

For R₃, trace one forbids XR₃=−R₃X. Its centralizer consists of norm-one
units, whose first expanding generator is ε₃ itself. The primitive norm
is 23+4√33, and finite torsion μ₆/{±1} has order 3. The inverse class
uses conjugation by J to transport this translation; it has the same
primitive norm and finite centralizer size.

The exact matrices are produced by `expected_classes()` and checked for
SL determinant, commutation, order, flip, subgroup closure, and inverse
intertwining. Minimality and maximality come from L3 and the axis argument,
not from finding matrices in an arbitrary bounded search.

The cylinder integral in D2-L7 gives log N(T₀)/(4 q sin²(π/m)), using
the **full** maximal finite centralizer size q. The involution endpoint
flip therefore gives denominator 16; each order-3 element class has
sin²(π/3)=3/4 and denominator 9. This proves the coefficient in the theorem.
The derivation and record-bound arithmetic replay protect this normalization;
B<1 cannot detect omitted positive contributions or an incorrectly lowered
coefficient. Independent trace-length coefficient checks are separate tests.

## D11-gap: scope of the analytic interface

The arithmetic inventory is complete for the full level-one group. It
registers through the existing inventory interface; the shared one-cusp
backend supplies volume, cusp, and scattering constants. The global shortest
trace and the rigorous finite cutoff are already proved in
[SYSTOLES.md](SYSTOLES.md), with d=11 witness τ and length
acosh((3+√45)/4). The support gate uses a rigorous lower endpoint of
that systole and downward rounding of the binary δ.

The analytic formulas and full-line tail bounds are in
[ANALYTIC_DERIVATIONS.md](ANALYTIC_DERIVATIONS.md). A spectral claim requires
a full Arb evaluation with B<1 after subtracting h(i), not just this
inventory theorem. Parameter-search results and any resulting certificate
are recorded separately; a finite comparison is not a global optimum.
