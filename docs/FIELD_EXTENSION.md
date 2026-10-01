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
| M1: complete exact class inventory | **Open** | All obligations in `INVENTORY_PROOF.md`, including full PSL normalization, resolved |
| M2: d=2 | Mechanical screen ready; spectral certificate blocked | Self-contained inventory, then frozen full certificate and example |
| M3: d=7 | Mechanical screen ready; spectral certificate blocked | Same, after d=2 |
| M4: d=11,19 | Mechanical screens ready; later work | Same, after d=2 and d=7 |

The new layout is `groups/` for explicit group inputs and exact arithmetic,
`fields/` for characters and field constants, and `core/` for analytic
evaluation and certificate gates. The original `bianchi_omega_arb.py` import
and command remain compatible. Its Eisenstein elliptic coefficient remains
log(7+4 sqrt(3))/8, labeled with the unresolved normalization issue.
`picard_stf.py` remains the independent historical quadrature implementation.

Install the trace dependencies and run:

```sh
python -m pip install -r requirements-trace.txt
python -m unittest discover -s tests -v
python bianchi_omega_arb.py
python picard_stf.py
python -m groups.systoles
python scripts/feasibility/flip_check.py --bound 2
python scripts/feasibility/gpp_check.py
python scripts/feasibility/screen.py 1 3 2 7 11 19
python scripts/feasibility/budget.py
python examples/field_screen.py 2 7 11 19
```

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
