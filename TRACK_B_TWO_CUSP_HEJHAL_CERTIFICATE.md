# Track B physical two-cusp Hejhal certificate

## Certified result

The first correctness-grade physical two-cusp system is certified at 192 Arb
bits with Fourier cutoff `M=1`.  This deliberately small cutoff establishes
the assembly, scaling, solve, and independent-verification pipeline before
larger numerical trials are attempted.

```text
physical matrix                         96 x 24 complex
normalized square system                24 x 24 complex
physical/reduced unknowns               24 / 24
all four cusp blocks                    certified
direct resolved Arb K evaluations       40 unique
unresolved Bessel fallbacks             0
verified contraction upper              1.084e-11
certified sigma_min lower               1.519e-4
maximum physical row residual upper     1.892e-5
physical residual L2 upper              2.581e-5
normalization residual upper            6.548e-16
independent reconstruction              certified
rung4_certified                         false
```

The independent verifier rebuilds all 96 rows and every Fourier-Whittaker
entry.  It does not read serialized matrix entries from the solver.

## Exact cusp conventions

The scaling matrices are

```text
sigma_infinity = I
sigma_0        = S = [[0,-1],[1,0]].
```

In PSL2, `S^-1=S`.  The cusp-infinity lattice is `Z[i]`.  The cusp-zero
lattice is `(2+i)Z[i]`, with dual lattice

```text
(1/conj(2+i)) Z[i] = ((2+i)/5) Z[i].
```

For a Gaussian index `a+bi`, the zero-cusp frequency is represented exactly
as

```text
u = (2a-b)/5,
v = (a+2b)/5.
```

The direct upper-half-space matrix action and the specialized cusp-zero
formula

```text
S(x1,x2,y)=(-x1/rho,x2/rho,y/rho)
```

agree in Arb balls on every exact rational collocation point.  The largest
recorded difference is below `3e-57`.

## Physical matrix

The assembly uses the exact six-copy identities

```text
F_c(P) - F_{pi_gamma(c)}(gamma P) = 0
```

for `gamma` in `T1,R,TiR,S` and all six exact coset copies.  `F_0` is the
infinity-cusp expansion.  `F_1,...,F_5` are the cusp-zero expansion evaluated
at the exact translations `0,...,4`.

Each nonzero mode is

```text
y K_(i r)(2 pi |mu| y) exp(2 pi i <mu,x>).
```

Every target component is evaluated at the true transformed horizontal
coordinate and height.  No source-height copying, height-matched proxy,
modeled inversion block, or synthetic scattering matrix remains.

The 96-row ledger contains

```text
infinity -> infinity     12 rows
infinity -> zero          4 rows
zero -> infinity          4 rows
zero -> zero             76 rows
```

The coefficient reduction is the exact sparse identity map.  The two cusp
coefficient systems remain independent unknowns; their consistency is
imposed by the physical cross-cusp rows rather than by an unsupported
coefficient identification.

## Validated Bessel backend

Resolved entries use `python-flint`'s `acb.bessel_k` at complex order `i r`.
The backend accepts an entry only when the enclosure is finite and its
mathematically real value independently contains zero in the imaginary
component.  Failure raises immediately.

The interface separately counts direct, asymptotic, shifted-order recurrence,
majorant, and failed evaluations.  This run used only direct resolved Arb
values.  A validated shifted-order derivative implementation is present for
later defect-gradient assembly.  No analytic majorant is substituted for a
resolved matrix entry.

As a regression against the old `kve` failure, the same resolved radial entry
at `r=6` and `r=8` gives disjoint Arb intervals.

## Scaling and verified solve

The homogeneous system is normalized by the exact physical equation

```text
a_infinity,(1,0) = 1.
```

Twenty-three independent physical rows, including every cusp block, and this
normalization form an explicit 24-by-24 system.  No Tikhonov regularization,
identity blend, random right-hand side, Gram normal equation, or unexplained
perturbation is used.

The only preconditioning is an explicit power-of-two diagonal map

```text
V_tilde = D_r V D_c,
a_tilde = D_c^-1 a,
a = D_c a_tilde.
```

`D_r`, `D_c`, and both inverses are serialized.  The coefficient vector is
mapped back to physical coordinates before any theorem-facing residual is
computed.

With a floating midpoint inverse `C`, Arb proves

```text
||I-C V_tilde||_infinity < 1.084e-11 < 1.
```

The Neumann/Rump contraction therefore proves invertibility and encloses the
unique normalized solution.  The physical automorphy rows have certified rank
at least 23; the complete normalized square system has certified rank 24.

## Independent verification

`track_b_two_cusp_verify.py` reads only the result definition, coefficient
balls, exact cusp data, and collocation ledger.  It independently:

1. checks the exact row set and six-copy permutations;
2. recomputes every transformed point;
3. reconstructs each Fourier-Whittaker row with direct Arb calls;
4. verifies the scaling round trip;
5. verifies that every selected normalized equation contains zero;
6. reconstructs the full unscaled physical residual;
7. checks all deterministic hashes and published upper endpoints.

The independently reconstructed worst row is `p3:S:c0`, the
`infinity->zero` transition.

## Conditioning diagnostics

These floating values are diagnostics, not proof inputs:

```text
physical matrix condition estimate       3.240e2
normalized subsystem estimate            2.281e7
power-of-two scaled estimate              7.586e3
face-only numerical nullity               12
complete-system numerical nullity          0
```

The nullity comparison confirms the expected qualitative point: face-only
constraints do not determine the vector, while the complete two-cusp system
plus normalization does.  Certified invertibility comes from the interval
contraction, not from these SVD diagnostics.

## Reproduction

```powershell
& 'C:\Users\Admin\AppData\Local\Programs\Python\Python313\python.exe' `
  track_b_two_cusp_hejhal.py `
  --bits 192 `
  --r-interval '6.7439020359331625,6.7439020359331625' `
  --fourier-cutoff 1 `
  --collocation-points 4 `
  --assemble-physical `
  --verified-solve `
  --independent-verify `
  --collocation-ledger track_b_hejhal_rows.jsonl `
  --json-out track_b_two_cusp_result.json `
  --verification-json-out track_b_two_cusp_verification.json

& 'C:\Users\Admin\AppData\Local\Programs\Python\Python313\python.exe' `
  track_b_two_cusp_verify.py `
  --result track_b_two_cusp_result.json `
  --collocation-ledger track_b_hejhal_rows.jsonl `
  --json-out track_b_two_cusp_verification.json
```

## Remaining work

The finite `M=1` trial establishes physical assembly correctness, not a
complete global Hejhal defect.  Rung 4 remains false pending:

1. cutoff growth and a certified Fourier truncation tail;
2. both cusp residuals on complete continuum cusp regions;
3. all non-floor pairing-face residuals and first gradients;
4. reprojection error and the complete partition-weighted defect;
5. the final interval comparison with the admissible global defect budget;
6. the later defect-to-spectrum/counting bridge.

No eigenvalue-existence claim is made by this certificate.
