# Morning brief — repaired dual-certification path

**Date:** 2026-07-12  
**Scope:** `dual_certification_` follow-up after the residual audit  
**Hard map:** `pipeline_ok=true`, `rung4_certified=false` (unchanged)

## Headline

The old `eta≈0.037` iteration has been retired as evidence about the spectrum.  Its Bessel backend silently replaced (K_{ir}) by an (r)-independent upper bound, and its residual imposed full-Picard scalar invariance on one cusp instead of the exact six-copy Γ0(2+i) coset relations.

A repaired diagnostic now:

- evaluates signed (K_{ir}(x)) with 160-bit python-flint/Arb and fails closed;
- agrees with 60-digit mpmath at all cross-check points;
- retains every amplitude, column, and equilibration diagonal;
- uses the exact (P^1(ℝ_5)) gluing permutations;
- represents one infinity component and five translated zero-cusp components;
- uses Λ0-dual = `(2+i)/5 Z[i]`;
- balances physical frequency cutoffs with `M0=5*M_inf`;
- samples the coset identity in a Nyquist-resolved 3-D interior cloud.

This repaired operator detects a sharp residual minimum at the published engineering target (r≈6.62212).  That is the expected level-1 Picard eigenfunction lifted to the congruence subgroup, not yet a certified first eigenvalue.

## New evidence

| check | result | status |
|---|---:|---|
| Old model matrix change, `r=6→8` | exactly 0 | RED, retired |
| Arb-repaired model change, `r=6→8` | 0.967 relative | GREEN backend check |
| Arb vs mpmath (K_{ir}) | agrees at displayed precision | GREEN backend check |
| Equilibration ledger reconstruction | 1.35e-16 relative | GREEN diagnostic |
| Exact six-copy operator rank, M=28 | 176/176 | GREEN diagnostic |
| Minimum location | near `r=6.6221` | strong cross-check |
| M=28 asymmetric-cutoff RMS / max | 1.80e-5 / 1.42e-4 | YELLOW diagnostic |
| Balanced M=8 RMS / max | 1.20e-5 / 7.26e-5 | YELLOW diagnostic |
| Balanced M=12 RMS / max | 1.45e-5 / 1.43e-4 | YELLOW, non-monotone |

At the M=28 minimum, the recovered two-cusp vector has the structure of a lifted level-1 form:

- zero-cusp coefficient norm off the integer-frequency sublattice: `6.48e-5`;
- matched infinity/zero coefficients disagree by `1.56e-4` after one scalar;
- best scalar is `0.99986 + 2.2e-15 i`.

This is independent evidence that the corrected coset operator is seeing the intended lifted eigenmode.

## What remains uncertified

The current residual is still not Theorem D(K)'s eta.  It encloses Arb Bessel values at floating sample coordinates, but phases, coordinates, SVD, between-sample control, and normalization are not interval-certified.  The M sequence is not monotone, so no truncation limit has been established.

Therefore:

```text
eta_le_eta0: false
width_lt_tol: false
counting_certified: false
multi_copy_kappa_certified: false
rung4_certified: false
```

## Revised priority order

1. **Stabilize the exact six-copy Hejhal discretization.** Use nested point sets, balanced cusp cutoffs, and independent grids. Require monotone/stable singular-vector and residual behavior before increasing M further.
2. **Create two distinct residual objects.** `DiagnosticResidual` may use floating SVD; `CertifiedDefect` must contain interval face/core norms, tail, normalization lower bound, and projection error matching D(K).
3. **Interval-isolate the lifted eigenparameter.** Enclose an (r)-interval by a Krawczyk/interval singular system with a normalization pin. The five-decimal Then value is only a center.
4. **Prove the residual transfer.** Bound the continuum norm between interior samples using derivative/Sobolev estimates and Arb (K_{ir}), phase, and coordinate evaluations.
5. **Continue independent lower/counting route.** Prove the glued-space κ bridge and obtain a certified counting enclosure below the lifted eigenvalue. Existence near 44.85 does not prove it is first.

## Files

| file | purpose |
|---|---|
| `verified_kir.py` | fail-closed Arb (K_{ir}) enclosure |
| `verified_hejhal_phase1.py` | backend and complete scaling ledger audit |
| `six_copy_hejhal.py` | exact-coset, two-cusp diagnostic collocation |
| `test_verified_hejhal.py` | backend/lattice regression tests |
| `verified_hejhal_phase1_result.json` | repaired model/scaling evidence |
| `six_copy_hejhal_*_result.json` | scan and truncation evidence |

## Bottom line

We have moved from a false (r)-independent plateau to a mathematically appropriate coset operator that locates the lifted Picard mode.  We have **not** reached dual certification.  The immediate blocker is now a precise one: construct a convergent and interval-certified continuum defect from the six-copy operator, then combine it with certified counting and the multi-copy κ bridge.
