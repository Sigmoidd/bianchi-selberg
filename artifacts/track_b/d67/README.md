# d67 start record

The d67 spectral target and geometry gate are OPEN. The frozen initial
design is [D67_TRACK_B_DESIGN.md](../../../docs/D67_TRACK_B_DESIGN.md).
Run from the repository root:

```sh
python artifacts/track_b/d67/verify_cover.py
python -O artifacts/track_b/d67/verify_cover.py
python artifacts/track_b/d67/negative_checks.py
python artifacts/track_b/d67/face_probe.py
```

The independent cover verifier checks 712 rational leaf squares and
unimodular witnesses, proves `H(z)>=1/34` across the full translation
cell, and proves that omitted denominators have norm at least 36.
Thus the finite candidate limits are `N(c)<=29` and `N(l)<169`.
The exact face probe reports 37 positive-area patches, total projected
area 1, and minimum enumerated vertex height squared `2/67`.
It remains diagnostic until face pairings and stabilizers have an
independent field-specific replay. `cover_verification.json` and
`negative_checks.json` state the precise passed and open gates.

The optimized interpreter reproduces the cover result byte-for-byte;
all four negative controls are rejected. The frozen d19 replay still
passes (`d19_replay.log`). Historical artifacts are unchanged.
