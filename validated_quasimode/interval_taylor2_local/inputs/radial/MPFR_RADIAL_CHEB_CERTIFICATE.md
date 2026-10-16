# Picard 40-orbit quasimode: interval-validated radial Chebyshev library

Date: 2026-10-15. **Status: RADIAL TABLE CHECKPOINT COMPLETE; 3D POINTWISE RESIDUAL CERTIFICATE STILL OPEN.**

For each of the 24 sum-of-two-squares integers `1<=n<=50`, the C/MPFR verifier produces outward-directed intervals for 121 Chebyshev-Lobatto nodes and all 121 coefficients of *two* degree-120 Chebyshev interpolants:

    F_n(y)=y K_{i r}(2 pi sqrt(n) y)
    G_n(y)=d F_n(y)/dy=K_{i r}(z)+z K'_{i r}(z), z=2 pi sqrt(n)y
    r=6622119348586011/1000000000000000 exactly.

The interval polynomial is represented in standard `numpy.polynomial.chebyshev` coefficient convention on `t=(15/8)*(y-17/15)`, `y in [3/5,5/3]`. Each midpoint coefficient is an exact rational decimal printed in the CSV; its interval bounds are outward decimal-enclosed MPFR values.

## Certified input quadrature

The standalone C code uses 512-bit MPFR with directed-rounding enclosures, loading and revalidating 132 Gauss-Legendre polynomial-root brackets and positive weights. The Bessel integral representation for real `r` and real `z>18/5` is

    K_ir(z)=int_0^infty exp(-z cosh u) cos(ru) du,
    K'_ir(z)=-int_0^infty cosh(u)exp(-z cosh u)cos(ru) du.

On the Bernstein ellipse `rho_u=13/10` of `[0,6]`, `|Im u|<=207/260<4/5`, `cos(Im u)>2/3`, `r<7`, and `z>18/5`. Thus the integrand and its `cosh(u)`-weighted counterpart are bounded by `<73`, with purely rational exponential comparison inequalities. The standard analytic Gauss-Legendre error bound is then

    abs(error) < 3000*(10/13)^264 < 10^-25.

The improper `u>6` tails are each `<10^-200` since `sinh 6>200` and `z>18/5`. Each computed K and K' is therefore enlarged by a *conservative* `+/-10^-24`, with outward MPFR arithmetic. The exact rational checks are in `audit_radial_theorem.py`.

## Global analytic interpolation bound

On the Bernstein ellipse `rho_y=3/2` around `[3/5,5/3]`,

    Re(y)>=5/9, |y|<=77/45<2, Re(z)>10/3.

Using `cosh(u)>=1+u²/2`, `cosh(u)<=exp(u²/2)`, the real-u Bessel integral, elementary Gaussian integration and Taylor bounds for `exp(10/3)>27`, one obtains uniform complex-ellipse envelopes

    |F_n(y)|<1/10, |G_n(y)|<4.

The usual degree-120 Chebyshev-Lobatto interpolation bound gives, for interpolants with **exact** node values,

    ||F_n-P_n||_infty < (4/5)*(2/3)^120 < 1e-20,
    ||G_n-Q_n||_infty < 32*(2/3)^120 < 3e-20.

C interval DCT-I turns the MPFR node enclosures into enclosures of all coefficients of these ideal interpolants. `check_radial_table.py` parses decimal endpoints as exact rationals via `decimal.Decimal`; over all 24 norms, the largest total coefficient-interval width is

    sum_k width(C_n,k for F_n) <=3.470e-22,
    sum_k width(C_n,k for G_n) <=1.573e-20.

Since `|T_k(t)|<=1` for real `|t|<=1`, the midpoint coefficient polynomials (read as exact-decimal polynomials) therefore satisfy, for **every real** `y in [3/5,5/3]` and every norm `n` listed,

    |F_n(y) - midpoint_P_n(y)| < 1e-20,
    |G_n(y) - midpoint_Q_n(y)| < 4e-20.

These constants are deliberately relaxed rational consequences of the printed interval widths and analytic interpolation estimates. Direct `mpmath` evaluations at held-out points agree numerically to ~1e-20, but are only an independent diagnostic.

## Replay

```bash
cd mpfr_radial
python audit_radial_theorem.py
# Required if regenerating node brackets:
# copy gl_n132.txt from the existing three-mode MPFR certificate replay
cc -O2 -std=c11 -o certify_radial_chebyshev_nodes_mpfr \
    certify_radial_chebyshev_nodes_mpfr.c -l:libmpfr.so.6 -lgmp -lm
./certify_radial_chebyshev_nodes_mpfr 0 24 radial_nodes_all.csv radial_coeffs_all.csv
OPENBLAS_NUM_THREADS=1 python check_radial_table.py
```

The source supports the standard MPFR header and a Linux x86_64 ABI fallback; the fallback is platform-specific. The exact root-bracket file is shipped unchanged. Zero theorem flags are promoted.

## Remaining decisive step

Use these validated radial polynomials and derivative polynomials inside **an outward-rounded 3D interval Taylor-model verifier** on `0<=t<=x<=1/2`. Cut at exact quadratic cutoff surfaces; certify **both**

    sup |W-f| < 1e-6 on the triangular strip,
    sup |(L-lambda)W| < 3.4e-10 on the triangular Ford core.

The first is an immediate mass-transfer gate. The second must preserve cancellation of the automorphy defects; naive interval numerator/denominator subtraction is likely too wide. The resulting `R/mu<1e-5` is **not yet proved**.
