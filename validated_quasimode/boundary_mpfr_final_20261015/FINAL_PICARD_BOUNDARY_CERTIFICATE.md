# Picard Q-odd cusp spectrum: completed boundary MPFR calculation

Date: 2026-10-15. Research branch: `Sigmoidd/bianchi-selberg`, draft PR #14. Status: **the outstanding 128×128 outward-rounded coefficient calculation is complete and passes**. The spectral conclusion follows from the explicit analytic lemmas and geometric/function-space argument below. This is a *computer-assisted, independently reproducible theorem candidate*, not a peer-reviewed or uniqueness result.

## The precise theorem being tested

Let \(\Gamma=\operatorname{PSL}_2(\mathbb Z[i])\) and \(L=-y^2(\partial_x^2+\partial_t^2+\partial_y^2)+y\partial_y\) be its positive hyperbolic Laplacian. Set

\[
 r=\frac{6622119348586011}{10^{15}},\quad \lambda_0=1+r^2.
\]

Assuming the standard Picard Ford fundamental half-domain and the exact reflection/Q symmetry identifications detailed below, the data in this package imply

\[
 \boxed{\sigma_{\mathrm{cuspidal},Q\text{-odd}}(L)\cap(\lambda_0-8\cdot10^{-6},\lambda_0+8\cdot10^{-6})\ne\varnothing.}
\]

This asserts **existence of some eigenvalue**, not uniqueness, index, absence of smaller eigenvalues, or a claim for \(\Gamma_0(2+i)\).

## The fixed interior-exact trial

On \(F=\{ -\tfrac12\le x\le\tfrac12,0\le t\le\tfrac12, y\ge\sqrt{1-x^2-t^2}\}\), let

\[
 f(x,t,y)=\sum_{(a,b,n,c)\in\mathcal T} c\,yK_{ir}(2\pi\sqrt n y)
 [\cos(2\pi ax)\cos(2\pi bt)-\cos(2\pi bx)\cos(2\pi at)].
\]

The 21-entry **exact rational** set \(\mathcal T\) is in `terms.txt`. `audit_exact_theorem.py` independently regenerates it from `analytic_quasimode_coefficients_50_reflection_even.json`, verifies quarter-turn oddness, reflection evenness, normalization and every term in the MPFR input.

The trial satisfies \((L-\lambda_0)f=0\) in the interior. It is periodic, even in \(x,t\), and odd under \((x,t)\mapsto(-t,x)\). Across the Ford hemispherical pairing \(S\), its value matches because on the sphere \(S(x,t,y)=(-x,t,y)\). Thus it is in the \(H^1\) form domain, but need not be in the operator domain; no normal-derivative matching is assumed. At the cusp the Fourier modes decay.

## Boundary-only weak residual

Let \(y_0=\sqrt{1-x^2-t^2}\), \(D=[-\tfrac12,\tfrac12]\times[0,\tfrac12]\). The **one-sided** outward hyperbolic flux and hyperbolic floor measure are

\[
 J=-y_0(xf_x+tf_t+y_0f_y),\qquad dS_H=\frac{dx\,dt}{y_0^3}.
\]

Green's identity and cancellation of all other face terms give

\[
 a(f,v)-\lambda_0(f,v)=\int_D J\overline{v}\,dS_H.
\]

The elementary collar trace bound is \(\|v\|_{L^2(\text{floor})}\le\frac92\|v\|_{H^1_\mathrm{form}}\). The floor area is \(<5/6\), and the existing independently replayed MPFR height-slab certificate gives \(\|f\|>13\cdot10^{-6}\). In the Q-odd reflection-even subspace the cusp zero Fourier mode is absent, yielding discrete spectrum. The standard form-dual spectral inclusion argument then shows that

\[
 \sup_D|J|<4\cdot10^{-12}
 \quad\Longrightarrow\quad \operatorname{dist}(\lambda_0,\sigma_{\mathrm{odd}})<10^{-5}.
\]

In fact the final stronger bound below suffices for \(8\cdot10^{-6}\), since \(7(9/2)\sup_D|J| < 13\cdot10^{-6}\cdot8\cdot10^{-6}\).

## The completed new machine calculation

**Reproduction of the previously missing gate:**

- `certify_radial_192.c`: standalone 1024-bit MPFR directed-rounding Bessel integrals, using **192 root-validated Gauss–Legendre points** on \([0,6]\), for all 24 radial norms and their derivatives. Uniform quadrature error \(<3000(10/13)^{384}<10^{-38}\) plus negligible truncated-tail error. All 192 root brackets, weights, and 24×121 Lobatto-node and coefficient intervals are written to CSV; positivity, ordering and Legendre sign changes are checked by MPFR.
- `certify_boundary_dct_mpfr.c`: standalone 256-bit MPFR *interval* evaluation of all **16,384** Chebyshev nodes. It evaluates trigonometric factors, \(y_0\), the degree-120 radial polynomials and the 21-term normal flux with outward rounding, then performs an interval tensor DCT-II in two passes. It **never trusts a floating-point DCT midpoint**.
- The full 128×128 run returned

\[
 \boxed{\sum_{k,l=0}^{127}\sup|c_{kl}|\le 3.054010740798494611699450816936\times10^{-12}<3.5\times10^{-12}.}
\]

Largest width of any nodal \(J\) interval was \(3.42563676\times10^{-33}\). The program returns `CERTIFICATE_DCT_L1_PASS` and exit code zero only if the strict DCT inequality holds. Smaller 8, 32 and 64 node tests pass. An independent, directly evaluated 45-decimal Bessel calculation agrees with the interval nodal enclosures on 3 check points.

### Uniform analytic remainders (separate from machine rounding)

1. For every required radial norm, the degree-120 Lobatto interpolants for \(F_n(y)=yK_{ir}(2\pi\sqrt n y)\) and \(G_n=F_n'\) satisfy \(\sup |F_n-P_n|,\sup |G_n-Q_n|<10^{-32}\) on \([3/5,5/3]\). This follows from the Bernstein ellipse \(\rho=2\), where \(\Re y\ge7/15\), \(|y|\le9/5\), \(\Re z>14/5\), \(|K_{ir}(z)|,|K'_{ir}(z)|<1/16\), \(|F_n|<1/8\), \(|G_n|<7\), and a standard Lobatto interpolation bound \(4M2^{-120}\). The full MPFR radial coefficient intervals enclose the ideal interpolation coefficients, so only this analytic remainder remains.
2. On the Ford floor, \(\sum|c|<37\), \(a+b\le9\), \(2\pi<7\), \(|T|\le2\), and \(|xT_x+tT_t|\le63\). Consequently the radial interpolation contributes less than \(2.405\times10^{-29}<3\times10^{-29}\) to each exact nodal value of \(J\), and at most \(4(128)^2(3\cdot10^{-29})<2\cdot10^{-24}\) to the DCT coefficient \(\ell^1\) sum.
3. For the true analytic \(J(x,t)\), a tensor polyellipse \(\rho_x=\rho_t=3/2\) gives \(\sup|J|<20,000,000\) in the complexified rectangle. The tensor degree-127 Gauss–Chebyshev interpolation error is bounded by \(144(2\cdot10^7)(2/3)^{128}<8.313\cdot10^{-14}\). Rational comparisons are in `audit_exact_theorem.py`. This bound uses the real-\(u\) integral representation of \(K_{ir}\), the principal analytic square root for \(y_0\), and classical Chebyshev aliasing estimates.

The final exact-rational audit gives

\[
 \sup_D|J|< 3.05402\cdot10^{-12}+2\cdot10^{-24}
 +144(2\cdot10^7)(2/3)^{128}
 <\boxed{3.14\cdot10^{-12}}.
\]

This is below the \(4\cdot10^{-12}\) theorem gate and meets the tighter \(8\cdot10^{-6}\) inclusion threshold.

## Independent replay

In a Linux environment with Python 3 (`mpmath`, `scipy` needed only to *regenerate* root brackets) and the runtime libraries `libmpfr.so.6`, `libgmp.so`:

```bash
python audit_exact_theorem.py
python audit_seed_slab_simpson_bounds.py
cc -O2 -std=c11 -o certify_radial_192 certify_radial_192.c -l:libmpfr.so.6 -lgmp -lm
cc -O2 -std=c11 -o certify_boundary_dct_mpfr certify_boundary_dct_mpfr.c -l:libmpfr.so.6 -lgmp -lm
./certify_radial_192 0 24 regenerated_nodes.csv regenerated_coeffs.csv
./certify_boundary_dct_mpfr regenerated_coeffs.csv terms.txt 128
```

For faster replay, the included `radial_coeffs_192.csv` and `gl_n192.txt` are frozen MPFR-generated results and the exact C code can be run directly on that CSV. All files' hashes are recorded in `MANIFEST_SHA256.json`. Re-run the 128 node transform after modifying source, precision or data.

## Scope and remaining checks

The *numerical gate* and analytic rational error-budget implications **are closed** in this checkpoint. The proof uses standard but nontrivial external geometric facts about the Picard fundamental domain and normalizer of the quarter-turn, the trace/weak-Green argument, the analytic Gauss and Chebyshev error theorems, and the validity of the signed MPFR interval implementation. They are written explicitly enough for independent audit, but **no external referee has yet endorsed them**. Do not call this a unique eigenvalue, a certificate of the lowest eigenvalue of the full group, or a solution to the \(\Gamma_0(2+i)\) / d=67 spectral exclusion problems.
