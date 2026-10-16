# Picard Q-odd quasimode: cancellation-aware 3D interval prototype

Date: 2026-10-15. Status: **LOCAL INTERVAL CHECKPOINT; GLOBAL NARROW RESONANCE/EIGENVALUE CERTIFICATE OPEN.**

## The new exact algebraic reduction

The fixed, real, `t`-even 40-orbit trial (source SHA in the shipped coefficient JSON) has only **21 nonzero real cosine-difference profiles**:

\[
f(x,t,y)=\sum_{a>b\ge0} c_{ab}F_{a^2+b^2}(y)
\big[\cos(2\pi a x)\cos(2\pi b t)-\cos(2\pi b x)\cos(2\pi a t)\big],
\qquad F_n(y)=yK_{ir}(2\pi\sqrt n\,y).
\]

The 21-term identity was checked against all 160 original signed Fourier coefficients *exactly as rationals*. It is an algebraic identity, not a fitted surrogate. All coefficients and the parameter `r=6.622119348586011` are frozen decimal rationals.

The radial second derivative is eliminated exactly by the Bessel ODE:

\[
F_n''(y)=4\pi^2nF_n(y)+F_n'(y)/y-(1+r^2)F_n(y)/y^2.
\]

This allows the existing MPFR-certified Chebyshev tables for `F_n,F_n'` to drive a full 3-variable interval jet carrying function value, 3 partials and all 9 second partials through each Möbius transform. The interval computations use `mpmath.iv` at 30 decimal digits and the certified interpolation error additions `±1e-20` and `±4e-20`. The input root/interval table provenance is reproduced in `inputs/radial/manifest_radial.json`.

For each nonidentity coset, construct the *centered* automorphy defect

\[
 D_i(X)=f(\gamma_iX)-f(X),\qquad
 D_i(X)\in D_i(m)+\nabla D_i(m)\cdot (X-m)
      +\tfrac12(X-m)^T H D_i(B)(X-m).
\]

The Hessian interval holds throughout the rational box `B`. This retains first- and second-order cancellation before forming interval hulls. The quotient correction `δ=W-f` and its gradient feed the identity

\[
R=\frac1G\sum_i[(L\psi_i)(D_i-\delta)-2y^2\nabla\psi_i\cdot(\nabla D_i-\nabla\delta)].
\]

`certify_local_mass_tiles.py` checks `|δ|<1e-6`; `local_residual_iv.py` checks `|R|<3.4e-10`. Both correctly return **UNDECIDED** on a box whose interval height straddles the interpolation domain boundary, rather than extrapolating outside `[3/5,5/3]`. Branch decisions about the quintic cutoff are made using interval comparisons; cutoff values use monotone endpoint enclosure.

## What actually passed

On the interior rational box centered at `(x,t,y)=(0.3,0.1,1.15)`:

| Interval check | Side length | Outward interval result | Decision |
|---|---:|---:|---|
| `W-f` | `0.02000` | `[-5.6785173e-7, +5.6785253e-7]` | strict local PASS |
| `R` | `0.00001` | `[-3.266161e-10, +3.175149e-10]` | strict local PASS |
| `R` | `0.00010` | `[-3.253344e-8, +3.252541e-8]` | UNDECIDED, NOT FAILED |
| `R` | `0` (point) | `[-4.623022e-12, -4.622594e-12]` | strict local PASS |

At the same center an independent 45-digit, **direct Bessel** evaluation (no Chebyshev tables, no interval-jet code) gives

`R=-4.6228078939377454744300978931470501e-12`, correctly inside the interval.

Changing the lowest radial Fourier coefficient by 1% yields a point residual near `-4.34e-7`, and the certificate gate **rejects** the mutation. Five source regression tests pass. The 21-profile identity uses exact rational coefficient matching, and is not based on floating-point comparisons.

## Why global certification did NOT close

The residual gate still loses substantial cancellation as soon as the box expands. At `side = 1e-4`, the interval is already approximately `±3.25e-8`, more than 95 times the allowable pointwise bound. Covering a positive-volume 3D Ford region by `side = 1e-5` cells is wholly infeasible. The second-order jet makes a few local cells provable, but is **not** a practical global certificate. Never extrapolate these local PASS records to the rest of the Ford cell.

The remaining theorem gates are unchanged:

* On `0<=t<=x<=1/2` and `1<=y<=5/3`, certify `sup|W-f|<1e-6` globally, thereby transferring the existing strip mass to `||W||>1.28e-5`.
* On the entire Ford wedge, certify `sup|(L-lambda)W|<3.4e-10`. Together with the exact residual-support volume `11/80`, this would give `R/||W||<1e-5`.

**Neither global bound has been obtained.** The spectral theorem flag must remain false.

## Recommended redesign: coefficientwise high-order Taylor models

Use order 12–16 **multivariate** Taylor polynomials for `D_i` on rational wedge cells. Compute *and subtract* the Taylor coefficients of `f(gamma_i X)` and `f(X)` **before** enclosing the remainder or taking a supremum. Certify the analytic remainder by a complex polydisc Cauchy bound or interval Taylor remainder with MPFR/Arb. Propagate the already-small Taylor model of `D_i` and `∇D_i` through the rational cutoff formulas. Split geometrically at the exact quadratic equations `psi_i(y)=0/1`. The checker needs an auditable full-cover record, all ambiguous cells recursively subdivided, and independent MPFR/Arb replay.

A rigorous **squared-residual integration** via exact vertical root splitting and analytic quadrature error is a viable alternative if interval suprema keep demanding microscopically small boxes. Neither path is completed in this checkpoint.

## Replay

Requires `mpmath`, `numpy`; no paid compute, no remote data. From this directory:

```sh
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_local_interval
python certify_local_mass_tiles.py
python local_residual_iv.py
python independent_mp_point.py
python seed_interval_prototype.py
```

`independent_mp_point.py` is **non-interval diagnostic only**. The default `seed_interval_prototype.py` does naive interval hulling, displayed only as a failure-control comparison. The full MPFR radial-table generator is available in the older archived radial certificate; this package includes the checked 24-norm coefficient table and its SHA manifest, not the full node-generation source.
