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
| New-field inventory completeness | **Open. Self-contained arithmetic proof required.** | `docs/INVENTORY_PROOF.md` |
| New-field spectral certificates | **Blocked**, even if the mechanical bound is < 1 | `GroupData.require_inventory`, `core/certificate.py` |

An Arb enclosure certifies the numerical evaluation of supplied inputs. It
does not certify a missing conjugacy classification or an unsettled group
normalization. Mechanical screens omit positive elliptic contributions and
cannot establish a spectral gap. The new certificate exporter refuses both
mechanical screens and historical inventories.

The historical normalization is left frozen at the user's request. The
existence of a flip must not be presented as a completed re-certification.
Future work proceeds in the order d=2, d=7, then d=11 and d=19, after the
inventory obligations are resolved.
