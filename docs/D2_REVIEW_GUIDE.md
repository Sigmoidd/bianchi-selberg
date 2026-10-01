# Reviewing and reproducing the d=2 proof

The result is for Γ=PSL₂(Z[√−2]) at level one. There are four elliptic
element classes, with coefficient

\[
 C_{\rm ell}=\frac3{16}\log(3+2\sqrt2)+\frac29\log(5+2\sqrt6).
\]

The full Arb evaluation gives B∈[0.42455184,0.42957479]<1. This excludes
discrete Laplace eigenvalues in (0,1). The derivation is
[D2_INVENTORY_PROOF.md](D2_INVENTORY_PROOF.md); this page is its review map.

## What supplies each input

| Input | Source | Role |
|---|---|---|
| Cofinite Kleinian trace formula | Friedman, [math/0612807v1](https://arxiv.org/abs/math/0612807v1), Theorem 4.1.1, printed pp.41–42 | Standard external identity |
| Minkowski bound and local valuation theory | General algebraic number theory; applied explicitly in L2–L3 | General tools, not a Bianchi classification table |
| Integral orders, lattices, units, and element multiplicity | Lemmas D2-L1–L6 | Self-contained field-specific arithmetic proof |
| Endpoint-exchange denominator | Lemma D2-L7 | Direct cylinder orbital integral |
| Systole and loxodromic annihilation | [SYSTOLES.md](SYSTOLES.md) and runtime Arb support check | Proved shortest trace; omitted loxodromic sum is zero |
| Field/scattering constants, integration, tails | [ANALYTIC_DERIVATIONS.md](ANALYTIC_DERIVATIONS.md) | Analytic reduction and full-line bounds |
| Numerical evaluation | python-flint 0.9.0 / Arb | Enclosed special functions, quadrature and inequalities |

No published class-count table, private OCR extract, or bounded conjugator
search supplies completeness. The mathematical arguments use ordinary
number theory and exact finite checks; this is not a proof-assistant
formalization.

## Read and replay the load-bearing steps

| Lemma | Mathematical conclusion | Replay in `groups/d2_inventory.py` |
|---|---|---|
| D2-L1 | O is a PID; one cusp; elliptics reduce to trace-zero or trace-one lattices | The Euclidean/Bezout/lattice arguments are written in the proof; base class number is independently checked by `groups.systoles.verify_all` |
| D2-L2 | Maximal S2 has class number one; every A2-stable lattice has one of two multiplier types, including the nonmaximal order | `verify_orders_and_class_numbers`, `verify_lattice_types` check discriminants, principal small primes, and every normalized conductor subspace |
| D2-L3 | A3 is maximal and has class number one; one trace-one GL lattice type | `verify_orders_and_class_numbers` checks the coprime discriminants and every possible small ideal norm |
| D2-L4 | Complete unit groups; determinant images {±1}, {±1}, {1} | `verify_unit_reduction` enumerates boxes whose bounds are proved in L4; integer radical comparisons select the reduced units |
| D2-L5 | Two involutions and two distinct inverse order-3 PSL element classes | Determinant-coset argument in the proof; `verify_class_witnesses` checks the determinant-minus-one inverse conjugator and the distinct involution reductions |
| D2-L6 | Individual minimal loxodromic norms and finite full-centralizer sizes 4,4,3,3 | Complete units establish minimality; `verify_class_witnesses` checks determinant, commutation, norms, flips and finite subgroup closure; the axis argument proves maximality |
| D2-L7 | Each involution uses denominator 16; both inverse order-3 classes use denominator 9 | Cylinder integral in the proof; independent coefficient identity in `tests/test_d2_inventory.py` |
| D2-gap | B<1 after subtracting the constant eigenfunction; no exceptional eigenvalues | `core.assemble.evaluate`, registered proof and analytic backends, `core.certificate.certificate_payload` |

The replay checks finite identities. For example, an expanding matrix alone
does not prove it is primitive: minimality uses the complete unit-group
argument. A four-element subgroup alone does not prove maximality: the
axis-action bound supplies that step. Read those arguments with the checks.

## Exact conventions to keep during review

- Coordinates (a,b) mean a+b s for s=√−2. Matrices are four such pairs in
  row-major order. Group-level ideals use the standard integral basis.
- Classes are conjugacy classes of **elements**, rather than cyclic
  subgroups or unoriented singular axes. The inverse order-3 class is a
  separate contribution because no norm-minus-one A3 unit removes its
  determinant obstruction.
- N(T₀) is the expanding eigenvalue modulus squared. The A2 norm is
  17+12√2, whereas the S2 norm is 3+2√2. They must not share a per-trace norm.
- E(R) is a maximal finite subgroup of the **full PSL centralizer**. The
  involution flips commute in PSL and anticommute on SL lifts.
- Every claimed new inventory is bound to its group key and proof id.
  Changing d, subgroup, level ideal or class witnesses invalidates that
  binding. The analytic backend independently checks its supported groups.

## Reproduce from a checkout

```sh
python -m pip install -r requirements-trace.txt
python -m groups.d2_inventory
python examples/group_certificate.py 2 --output /tmp/d2-certificate.json
python -m unittest discover -s tests -v
```

The shorter `python examples/d2_certificate.py` command delegates to the
same uniform CLI. The v1 [frozen certificate](../certificates/d2-k2.json)
remains unchanged. The [v2 certificate](../certificates/d2-k2-v2.json) adds
exact group identity and backend metadata; the numerical parameters,
arithmetic manifest and both endpoint balls are identical.

The original Picard and Eisenstein inputs retain legacy status. The
Eisenstein normalization is still frozen at the user's request. Neither
this proof nor the new interface asserts inventories for d=11,19 or a
new congruence-level trace formula.
