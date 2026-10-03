# Current trace-formula proof ledger

This is the current entry point. `old_RIGOR_GAPS.md` is a historical record,
including its original Eisenstein normalization, and is retained unchanged.

| Input or claim | Current status | Evidence |
|---|---|---|
| Picard/Eisenstein numerical baselines | Reproduced; historical group data retained | `certificates/legacy-pre-m0.json`, `certificates/m0-regression.json` |
| `g''(0)` for all integer k >= 2 | Exact polynomial derivative; previous stencil was only exact at k=2 | `core/bspline.py`, `docs/ANALYTIC_DERIVATIONS.md` |
| Loxodromic annihilation | Closed-form systoles, exhaustive finite trace check and analytic tail; Arb support comparison at runtime | `groups/systoles.py`, `docs/SYSTOLES.md` |
| Characters, prime splitting, field constants | Field-specific arithmetic, independent splitting tests | `fields/quadratic.py`, `core/terms.py` |
| Quadrature and analytic tails | Arb evaluation, both half-lines included; elementary digamma majorant | `core/assemble.py`, `docs/ANALYTIC_DERIVATIONS.md` |
| Eisenstein full PSL centralizer | Exact Klein-four witness established; normalization remains open | `scripts/feasibility/flip_check.py`, `docs/NORMALIZATION_ISSUE.md` |
| Historical elliptic class counts | Existing external classification inputs retained; not replaced by a new proof | `old_RIGOR_GAPS.md`, `docs/REFERENCES.md` |
| d=2 inventory completeness | Four element classes; complete lattice and unit proofs; full PSL orbital factor | `docs/D2_INVENTORY_PROOF.md`, `groups/d2_inventory.py` |
| d=2 spectral certificate | B in [0.42455184, 0.42957479], excluding discrete eigenvalues in (0,1) | `certificates/d2-k2.json`, `examples/d2_certificate.py` |
| Uniform identity and proof/analytic dispatch | Exact field/subgroup/ideal binding; full-group formulas reject proper subgroups | `docs/ADDING_GROUPS.md`, `tests/test_group_interfaces.py` |
| d=7 inventory and certificate | Complete two-class proof; B in [0.36394687, 0.36400056] | `docs/D7_INVENTORY_PROOF.md`, `certificates/d7-k2.json` |
| d=11,19 inventory completeness | **Open. Self-contained arithmetic proofs required.** | `docs/INVENTORY_PROOF.md` |
| Other new-field spectral certificates | **Blocked**, even if the mechanical bound is < 1 | `GroupData.require_inventory`, `core/certificate.py` |

An Arb enclosure certifies the numerical evaluation of supplied inputs. It
does not certify a missing conjugacy classification or an unsettled group
normalization. Mechanical screens omit positive elliptic contributions and
cannot establish a spectral gap. The new certificate exporter refuses both
mechanical screens and historical inventories.

The historical normalization is left frozen at the user's request. The
existence of a flip must not be presented as a completed re-certification.
The d=2 obligations are now closed by the local arithmetic proof and exact
replay. d=7 is also closed by its independent proof. Future work proceeds to d=11 and d=19. These are ordinary
mathematical proofs with executable checks, not proof-assistant formalizations.
The current d=2 schema-v2 artifact is `certificates/d2-k2-v2.json`; its
arithmetic manifest and numerical endpoints match the frozen v1 artifact.

## d=7 closure and finite-quotient boundary (2026-10-01)

The full d=7 inventory, primitive translations and centralizers are proved
in `docs/D7_INVENTORY_PROOF.md` and bound to `d7-arithmetic-v1` replay.
`certificates/d7-k2.json` gives B in [0.36394687, 0.36400056] <1.
The 55-point search is reproducible and only claims the smallest upper bound
among sampled sinc parameters, not a global optimum. d=11,19 remain open.
`quotients/` implements exact finite images and adjacency actions; graph
expansion bounds and congruence-cover Laplace bounds remain future work.
