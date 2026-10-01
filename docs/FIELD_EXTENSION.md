# Field-extension execution record

The selected decisions are:

- Preserve the historical Eisenstein elliptic normalization and note the issue.
- Require self-contained new-field inventory proofs; cited classifications
  serve as cross-checks only.
- Replace unavailable OCR-file dependencies with stable citations and local
  derivation notes.

| Milestone | State | Exit condition |
|---|---|---|
| M0: interface and safety net | Implemented | Two legacy enclosures reproduced within numerical quadrature variation; exact derivative; closed-form systoles; support gate; regression tests |
| M1: complete exact class inventory | Closed for d=2; per-field proof replay implemented | Four element classes; nonmaximal lattice included; SL splitting, complete units, primitive norms, full PSL orbital factor proved in `D2_INVENTORY_PROOF.md` |
| M2: d=2 | Full certificate frozen | `certificates/d2-k2.json`, B in [0.42455184, 0.42957479]; `examples/d2_certificate.py` |
| M3: d=7 | Mechanical screen ready; spectral certificate blocked | Same, after d=2 |
| M4: d=11,19 | Mechanical screens ready; later work | Same, after d=2 and d=7 |

The new layout is `groups/` for explicit group inputs and exact arithmetic,
`fields/` for characters and field constants, and `core/` for analytic
evaluation and certificate gates. The original `bianchi_omega_arb.py` import
and command remain compatible. Its Eisenstein elliptic coefficient remains
log(7+4 sqrt(3))/8, labeled with the unresolved normalization issue.
`picard_stf.py` remains the independent historical quadrature implementation.

The uniform interfaces are now documented in [ADDING_GROUPS.md](ADDING_GROUPS.md):
canonical field/subgroup/ideal keys, registered arithmetic replay and an
independent analytic backend. The d=2 review map is
[D2_REVIEW_GUIDE.md](D2_REVIEW_GUIDE.md). Higher-level backends can reuse the
shared test functions, support gate, elliptic/identity assembly and export;
their cusp and scattering proofs must be supplied separately.

Install the trace dependencies and run:

```sh
python -m pip install -r requirements-trace.txt
python -m unittest discover -s tests -v
python bianchi_omega_arb.py
python picard_stf.py
python -m groups.systoles
python -m groups.d2_inventory
python examples/group_certificate.py 2 --output /tmp/d2-certificate.json
python scripts/feasibility/flip_check.py --bound 2
python scripts/feasibility/gpp_check.py
python scripts/feasibility/screen.py 1 3 2 7 11 19
python scripts/feasibility/budget.py
python examples/field_screen.py 2 7 11 19
```

The d=2 inventory has two involution classes with multiplier orders A2 and
S2, and two inverse order-3 element classes. Its primitive norms are
17+12 sqrt(2), 3+2 sqrt(2), and 5+2 sqrt(6), respectively. Both involutions
have maximal finite full PSL centralizer of order four. The norm map on
A3 units has image {1}, so determinant-minus-one GL conjugacy does not
merge the inverse order-3 classes in SL or PSL. All class records are
bound to an exact replay, at assembly and at certificate export.

The new arithmetic proof does not replace the historical Picard/Eisenstein
classification inputs. Their regression bounds and coefficients remain
unchanged. Future fields require their own self-contained replay before
the full analytic assembly admits them.

The four attached feasibility scripts are in `scripts/feasibility/`. The
mpmath screen and budget remain exploratory, use truncated quadrature, and
are clearly separated from the Arb screen. `flip_check.py` discovers the
checkout from its own path and accepts `--repo`; it runs from an unrelated
working directory and exits nonzero on failed exact identities.

Frozen pre-M0 reports preserve the original numerical inputs and enclosure
at the base commit. Current regression reports use the closed-form systole
and Arb sinc, which slightly change floating point inputs and quadrature
radii. This is why equality of every endpoint is not an appropriate test:
the mathematical constants and enclosing intervals, rather than a particular
adaptive partition, are the reproducibility criteria.

For d=43,67,163, the mechanical lower endpoint already exceeds 1 for the
tested k=2, frac=0.999 parameters. This rules out those particular tests,
not every member of the sinc^(2k) family or every possible admissible
function. A universal impossibility claim needs a separate optimization
argument. Positive NCE budgets for the smaller fields also do not prove
feasibility without the inventory.
