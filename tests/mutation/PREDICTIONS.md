# Predictions written before mutation execution

Baseline parameters stay frozen: d=2 k=2, frac=.999, R=40; d=7 k=2,
frac=1, R=256. These predictions were written after baseline regression
replay and before any mutation execution.

| Mutation | Expected change in B | Reason |
|---|---|---|
| Drop each elliptic class | Decrease | Each noncuspidal class adds a positive coefficient times g(0). |
| Multiply the centralizer denominator by 2 | Decrease | Divides the positive NCE contribution by 2. |
| Multiply the centralizer denominator by 1/2 | Increase | Doubles the positive NCE contribution. |
| Remove each primitive translation witness; remove all witnesses (shrink set to empty) | Zero numerical change, arithmetic rejection | Coefficients use the separately stored primitive norm, not the witness matrix; arithmetic verification requires the matrix. |
| Volume +1% | Increase | Identity term is positive and linear in volume. |
| Volume -1% | Decrease | Identity term is positive and linear in volume. |
| Truncate exhaustive loxodromic trace box from [-6,6]^2 to [-1,1]^2 | Zero numerical change if accepted; completeness obligation violated | The shortest witness survives; the analytic pipeline has no enumerated geodesic sum because support is at/below the systole. |
| Omit cusp/scattering block | Increase | Baseline CE+Ch0+PARg0+PSI+PHIINT is strictly negative for both groups. |
| Flip cusp/scattering block sign | Increase | Replaces the strictly negative baseline block by its positive opposite. |

The cusp/scattering block includes CE, Ch0, PARg0, PSI and PHIINT; the
`prime` diagnostic is already contained in PHIINT and is not added twice.

Native checks mean the production pipeline's arithmetic inventory replay,
systole completeness verifier, geometry positivity checks, and certificate
B<1 gate, all invoked by evaluate/certificate_payload. A comparison against
the frozen numerical output is a separate regression check, not a native
mutation rejection. A passed mutated exporter is not evidence for a valid
theorem after trusted backend code has been monkeypatched.

If arithmetic verification rejects a changed record before B is computed,
a separate diagnostic evaluate call bypasses only require_inventory, in
memory, to measure the affected term. It does not export a certificate.
The results must distinguish these diagnostic intervals from intervals
returned by an accepted native pipeline.

Centralizer scaling is exercised both by altering the stored finite-centralizer
orders and by scaling the coefficient denominator in memory with the exact
class records intact. Both have the same predicted numerical direction; the
record mutation is expected to be rejected by arithmetic binding.
