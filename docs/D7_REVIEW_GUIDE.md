# Reviewing the d=7 proof

Start with [D7_INVENTORY_PROOF.md](D7_INVENTORY_PROOF.md). The complete
argument uses maximal orders, their Minkowski ideal representatives,
proved unit reduction boxes, determinant splitting, and the full PSL
centralizer. The matrix checks alone would not prove completeness.

| Obligation | Mathematical argument | Executable evidence |
|---|---|---|
| Base PID, one cusp, GG=1 | L1 Euclidean error 11/16; Bezout | Ring arithmetic and registered cusp records |
| All trace-zero/trace-one lattices | L1 correspondence; L2 maximality and class numbers | `verify_orders_and_lattices` |
| Complete unit groups | L3 archimedean bounds; reduced units are torsion | `verify_units` |
| Element class multiplicity | L4 determinant images and sign/inverse identifications | `verify_witnesses_and_splitting` |
| Primitive lengths | L5 norm-one exponent 1 or 2 | Unit powers and individual norm witnesses |
| Full centralizer coefficient | L5 no anti-commuting SL lift; D2-L7 cylinder integral | Finite subgroup closure; determinant images |
| Systole/support | `SYSTOLES.md`; Arb support inequality | Shared systole verifier and certificate gate |
| Analytic enclosure | `ANALYTIC_DERIVATIONS.md` | `certificates/d7-k2.json` |
| Parameter comparison scope | Last section of proof | 55-point manifest; global optimality flag is false |

`tests/test_d7_inventory.py` exercises class/denominator/centralizer tampering,
exact rational radicals and certificate support. d=2 and historical-field
regressions remain in the full suite. The d=7 proof backend registers through
the same interfaces as d=2; no new branch was added to trace assembly.

The discrete finite images in `quotients/` are a separate output. They do
not turn this level-one trace theorem into an expansion theorem. See
[FINITE_QUOTIENTS.md](FINITE_QUOTIENTS.md) for the operator contract.
