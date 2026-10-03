# d=11 predictions, recorded before execution

Run the existing `run.py` and `translation_consumers.py` algorithms with d=11
selected; only the field/parameter dispatch and CLI group selection are added.
Baseline parameters: k=2, frac=1, R=256; frozen d11-k2.json. Stop if any frozen
baseline B or endpoint ball fails exact replay. No production math or frozen
artifact changes.

| Existing mutation | Expected B direction | Reason |
|---|---|---|
| Drop each of the three elliptic classes | Decrease, visibly | Each positive log-norm coefficient contributes positively to NCE; native arithmetic should reject before B, diagnostic bypass should move B |
| Centralizer denominator ×2 (formula or records) | Decrease | Coefficient is inverse in the denominator; record edits should be rejected |
| Centralizer denominator ×1/2 (formula or records) | Increase | Positive elliptic coefficient doubles; record edits should be rejected |
| Remove each primitive translation witness or all | Native rejection; diagnostic unchanged | Arithmetic consumes the witness but numerical coefficient uses its separate stored norm |
| Suppress each translation's executing coefficient or all | Decrease, visibly | Removes the positive contribution at the actual NCE consumer; counter must execute |
| Volume +1% | Increase | Positive identity contribution is linear in volume |
| Volume −1% | Decrease | Same positive linear identity term |
| Shrink trace box [-6,6] to [-1,1] | Unchanged at frozen support | Proven shortest trace survives; no geodesic lies in support; this cannot validate truncated completeness |
| Omit cusp/scattering block | Increase | Baseline aggregate CE+Ch0+PARg0+PSI+PHIINT is negative |
| Flip cusp/scattering sign | Increase | Negating that negative aggregate raises B twice as much |

A native arithmetic rejection before B is unavailable B, not a zero-movement
certificate. Any gate-bypassed B is a diagnostic, never a certificate. Uncaught
undercounts can satisfy B<1, which cannot establish completeness or coefficients.

For the existing class-1 inside-support probe, shorten the first omitted
trace length to support/2: the complete verifier should reject it and the
truncated verifier should miss it (counter zero). This synthetic inconsistent
length fixture does not claim an actual shorter group geodesic.
