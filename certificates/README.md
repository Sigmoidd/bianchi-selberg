# Frozen baseline and current reports

Legacy regression reports and mechanical screens have
`spectral_certificate: false`. The new d=2, d=7 and d=11 results have
`spectral_certificate: true` and include their arithmetic proof manifests.

| Report | Contents |
|---|---|
| `legacy-pre-m0.json` | Unmodified engine at source commit `7c96c49fe4d4241347c219d7bc6ba48c286ff920`; original systole input and elliptic coefficient frozen |
| `m0-regression.json` | Current two-field assembly; exact systoles, Arb sinc, and unchanged historical elliptic normalization |
| `mechanical-k2.json` | Arb mechanical bounds for d=2,7,11,19,43,67,163; elliptic terms deliberately omitted |
| `d2-k2.json` | Full d=2 certificate, complete exact inventory and centralizers, with B in [0.42455184, 0.42957479] |
| `d2-k2-v2.json` | Same d=2 numerical result and arithmetic manifest, with canonical group identity, volume and backend metadata |
| `d7-k2.json` | Complete d=7 inventory and full certificate: B in [0.36394687, 0.36400056] |
| `d7-search.json` | All 55 sinc-family evaluations and selected parameters; no global optimum claim |
| `d11-k2.json` | Complete d=11 inventory and full certificate: B in [0.52924365, 0.52928221] |
| `d11-search.json` | All 55 d=11 evaluations; best-of-55, no global optimum claim |

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
python examples/trace_report.py i omega --output /tmp/current-regression.json
python examples/trace_report.py 2 7 11 19 43 67 163 --mechanical --output /tmp/current-screens.json
python examples/group_certificate.py 2 --output /tmp/d2-k2-v2.json
python examples/d7_certificate.py --output /tmp/d7-k2.json
python examples/d7_search.py --output-dir /tmp/d7-search
python examples/d11_certificate.py --output /tmp/d11-k2.json
python examples/d11_search.py --output-dir /tmp/d11-search
```

Reproducing `legacy-pre-m0.json` uses the original module at its recorded
base commit. The original `old_RIGOR_GAPS.md` remains unchanged. The current
proof ledger is `RIGOR_GAPS.md`; the new-field completeness boundary is in
`docs/INVENTORY_PROOF.md`.

The d=2, d=7 and d=11 proofs are in their `docs/D*_INVENTORY_PROOF.md` notes;
the d=19 obligations remain open. The frozen d=2 file includes all four elliptic
element representatives, primitive translation matrices, flip witnesses,
unit reduction bounds and class-number checks. Re-exporting it reruns the
arithmetic proof and the full Arb assembly. Its theorem statement excludes
discrete Laplace eigenvalues in (0,1) at level one.
