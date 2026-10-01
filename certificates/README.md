# Frozen baseline and current reports

These are numerical reports, not new self-contained theorem certificates.
The exporter makes that distinction explicit with `spectral_certificate: false`.

| Report | Contents |
|---|---|
| `legacy-pre-m0.json` | Unmodified engine at source commit `7c96c49fe4d4241347c219d7bc6ba48c286ff920`; original systole input and elliptic coefficient frozen |
| `m0-regression.json` | Current two-field assembly; exact systoles, Arb sinc, and unchanged historical elliptic normalization |
| `mechanical-k2.json` | Arb mechanical bounds for d=2,7,11,19,43,67,163; elliptic terms deliberately omitted |

All endpoint strings retain Arb radii. Hexadecimal delta strings preserve
the exact binary floating point test-function input.

| Field | Pre-M0 enclosure (approximate) | M0 enclosure (approximate) |
|---|---|---|
| Picard | [0.3045761085, 0.3178103502] | [0.3045784077, 0.3178080510] |
| Eisenstein | [0.5252159567, 0.5436742803] | [0.5252155528, 0.5436746842] |

The difference reflects the exact systole and certified quadrature
implementation, not a changed Eisenstein elliptic coefficient. No claim
that the full centralizer issue is settled is attached to these reports.

Regenerate the current reports:

```sh
python examples/trace_report.py i omega --output certificates/m0-regression.json
python examples/trace_report.py 2 7 11 19 43 67 163 --mechanical --output certificates/mechanical-k2.json
```

Reproducing `legacy-pre-m0.json` uses the original module at its recorded
base commit. The original `old_RIGOR_GAPS.md` remains unchanged. The current
proof ledger is `RIGOR_GAPS.md`; the new-field completeness boundary is in
`docs/INVENTORY_PROOF.md`.
