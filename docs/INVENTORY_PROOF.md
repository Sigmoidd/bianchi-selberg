# Self-contained elliptic inventory: required proof and present boundary

The user's chosen standard is a self-contained arithmetic proof for every
new group. Published classification tables may be cross-checks, but cannot
replace the completeness argument. The existing Picard/Eisenstein data are
explicitly marked `legacy`; this change does not retroactively satisfy the
new standard for them.

## What is already established

The field constants, one-cusp checks, unit counts, systoles and mechanical
Arb screens are available. For d=2,7,11,19 the only possible nontrivial
finite orders in PSL are 2 and 3: an elliptic SL lift has real trace in
O_K intersect R=Z, so its noncentral trace is 0 or +/-1. Their characteristic
polynomials are X^2+1 and X^2-X+1 (up to negating the lift). Neither has a
root in these K, so the fixed endpoints are not K-rational cusps.
There are no cuspidal elliptics because the only units are +/-1.
This restricts the types; it does **not** count their conjugacy classes.

`EllipticClass` records coefficients separately for each element-conjugacy
class. `QuadraticRing` provides exact multiplication in each integral basis.
The exact Klein-four witness addresses one local centralizer issue.
The complete d=2 list and proofs are now in
[D2_INVENTORY_PROOF.md](D2_INVENTORY_PROOF.md). There are two trace-zero
lattice types and one trace-one GL type which splits into two SL types.
The exact replay includes complete unit groups, primitive translations,
full PSL finite centralizers, and the four element-class records.

The complete d=7 proof is in [D7_INVENTORY_PROOF.md](D7_INVENTORY_PROOF.md):
two maximal orders, class number one, complete units and two PSL element classes.

## Proof obligations (closed for d=2 and d=7; next, d=11)

1. **All integral embeddings.** For each of the two characteristic
   polynomials, describe the associated quadratic extension L/K and every
   relative order that can occur as the multiplier order of an O_K lattice
   stable under the root. Prove that the list includes nonmaximal orders
   and all local conductor possibilities, especially at primes above 2 or 3.
2. **Every global lattice class.** Establish the embedding/lattice
   correspondence and enumerate the relevant proper ideal or lattice
   classes with a proved finite reduction bound. Supply representatives
   and the exact arithmetic needed to check exhaustiveness. A class number
   of K alone does not count these relative embeddings.
3. **SL2 rather than GL2 conjugacy.** Track determinants and the norm map
   on centralizer units. Determine which GL2 classes split into distinct
   SL2 classes. Prove identifications using determinant-one witnesses;
   never silently discard determinant obstructions.
4. **PSL element multiplicity.** Account for +/- lifts and powers. An
   unoriented singular axis, a finite subgroup and an individual elliptic
   conjugacy class are different objects. Prove which inverse classes
   merge and sum the trace formula over exactly its required classes.
5. **Individual integral centralizers.** For each representative M, solve
   XM=MX and, where allowed, XM=-MX, with det(X)=1 and X integral. The
   commuting algebra is K[M], but its integral intersection can be larger
   than O_K[M]. Do not substitute a shared per-trace unit search.
6. **Finite subgroup and primitive translation.** Prove the finite
   subgroup order and the minimum loxodromic norm for each class by a unit
   group computation with certified completeness. Record witnesses as
   exact algebraic data; a shortest element in a bounded search is not a
   primitive-unit proof.
7. **Orbital normalization.** Derive the involution orbital factor using
   the full PSL centralizer, including any endpoint exchange. Compare the
   result with the cited theorem's convention and resolve the Eisenstein
   witness discrepancy before transferring a formula to new groups.

Each resulting class must carry its representative, centralizer equations,
finite subgroup witnesses, primitive-unit proof, conjugacy disposition,
and precise links to the completeness derivation. A final theorem ties
those records to the trace-formula sum. Only then mark the inventory
`self-contained` and each normalization `proved`.

The code admits d=2 and d=7 only after replaying their complete records,
and rejects full d=11,19 assembly until their proofs are supplied. An empty
inventory does not mean an elliptic-free group. The legacy search remains
available for exploratory witnesses and is now labeled as bounded evidence.
Neither stable component counts nor agreement with a Weyl coefficient
replaces these arithmetic proof obligations.
