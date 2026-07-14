# Six-copy continuum defect enclosure

**Status:** rigorous Arb enclosure for one fixed finite trial; **not admissible to the current single-cusp Theorem D(K)**.  No hard-map flag changes.

## Certified object

The coefficient vector in `six_copy_hejhal_balanced_coeffs.json` is interpreted as an exact vector of decimal complex numbers.  It defines six analytic components over the Picard reference cell:

- one expansion at infinity with Gaussian modes `N(beta)<=8`;
- five translated expansions at cusp zero with dual-lattice modes `(2+i) beta/5`, `N(beta)<=40`;
- exact `P1(F5)` gluing permutations for `T1`, `R`, `TiR`, and `S`.

For every component and generator, `continuum_defect_arb.py` encloses

\[
F_c(p)-F_{\pi_\gamma(c)}(\gamma p)
\]

on a box cover of

\[
[-1/2,1/2]^2\times[1/\sqrt2,5/4].
\]

This rectangle contains the truncated Humbert core, so its supremum is conservative.  Generator pullbacks, phases, Bessel values, gradient bounds, and normalization are evaluated with python-flint at 160 bits.

## Enclosure method

Direct natural interval evaluation loses the modal cancellations and becomes unusably wide.  The production result therefore uses a centered first-order Taylor enclosure on each box:

\[
|R(p)|\le |R(p_c)|+
\|\nabla F_c\|_\infty |p-p_c|+
\|\nabla F_j\|_\infty |\gamma p-\gamma p_c|.
\]

The gradient majorants are analytic termwise bounds.  In particular,

\[
|K_{ir}(x)|\le \sqrt{\frac\pi{2x}}e^{-x},\qquad
|K'_{ir}(x)|\le K_1(x),
\]

the second inequality following directly by differentiating the integral representation and taking absolute values.  Exact point-center values use Arb `acb.bessel_k`.

All floats exported from Arb are rounded one ULP outward.

## Normalization

A positive `L2` lower bound is proved directly on the actual Picard half-torus

\[
[-1/2,1/2]\times[0,1/2]\times[1.001,1.249]\subset K_{5/4}.
\]

The planar integral uses the exact Fourier Gram matrix: full-period
orthogonality in (x_1) and the analytic half-period integral in (x_2).
Each Hermitian Gram block's smallest eigenvalue is enclosed by Arb's Rump
eigensolver. This avoids assuming that the recovered component is exactly
rotation invariant. The (y)-integral is split into 256 intervals. The
certified result is

\[
\|F\|_2^2\ge 0.0034510513425441325,\qquad
\|F\|_2\ge 0.058745649562704914.
\]

Only the infinity component is used, so omitted copies can only increase the true norm.

## Convergence record

| subdivision per coordinate | boxes | raw defect upper | normalized defect upper |
|---:|---:|---:|---:|
| 4 | 64 | 10404.587970610486 | 177112.8788415041 |
| 8 | 512 | 2965.556077348526 | **50482.0339085347** |

The enclosure decreases under refinement. A separate 128-bit subdivision-4
run before the normalization correction gave raw upper `10404.588022464937`,
agreeing with the 160-bit raw computation to about `5e-9` relative. The
normalized quotient remains an Arb ball through final endpoint extraction.

The present bound is far too pessimistic for certification: it sums absolute modal derivative contributions and therefore discards the cancellations responsible for the sampled defect near `1e-5`.  A higher-order Taylor model or validated spectral interpolation is required for a useful bound.

## PDE residual

For the explicit finite trial,

\[
\tau_{\rm core}=0
\]

analytically: each Whittaker character (yK_{ir}(2\pi|\mu|y)e^{2\pi i\langle\mu,z\rangle}) solves ((\Delta-(1+r^2))f=0).  This is a symbolic basis identity, not a floating residual estimate.

## Why this does not certify an eigenvalue

The current `theorem_DK.tex` explicitly treats the single-cusp group `PSL2(O_K)`.  The object above is a two-cusp, six-copy induced representation for `Gamma_0(2+i)`.  Before this defect can be passed to `defect_to_lambda_error`, the following theorem is still required:

1. a multi-copy/two-cusp version of D(K);
2. constants accounting for six-copy face overlap and both cusp collars;
3. proof that the vector-valued defect correction descends to the subgroup quotient;
4. the same explicit normalization and threshold bookkeeping.

Even after that extension, the current normalized upper bound `5.05e4` is vastly above `eta0`. Therefore:

```text
eta_le_eta0: false
width_lt_tol: false
counting_certified: false
rung4_certified: false
```

## Reproduce

```powershell
& 'C:\Users\Admin\AppData\Local\Programs\Python\Python313\python.exe' -B `
  continuum_defect_arb.py `
  --trial six_copy_hejhal_balanced_coeffs.json `
  --subdivisions 4,8 --bits 160

& 'C:\Users\Admin\AppData\Local\Programs\Python\Python313\python.exe' -B `
  -m unittest -v test_verified_hejhal.py
```

Machine record: `continuum_defect_arb_result.json`.  The trial file SHA-256 is embedded there so the certificate cannot be silently paired with different coefficients.

## Independent review

An adversarial sub-agent review independently confirmed the exact F5 gluing,
all four upper-half-space actions, the centered Lipschitz enclosure, the
derivative inequality, domain coverage, coefficient semantics, and the
termwise PDE identity. It identified two defects in the first draft:

1. full-torus Parseval was not a proved quotient-norm lower bound;
2. the final normalized quotient was performed in binary64.

Both are corrected in the recorded result: normalization now integrates the
actual area-1/2 half-torus through Rump-enclosed Fourier Gram eigenvalues, and
the quotient remains an Arb ball. The theorem-extension and
enclosure-pessimism objections remain open and are reflected in
`theorem_DK_admissible=false`.
