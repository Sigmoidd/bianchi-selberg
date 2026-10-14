# Independent MPFR and mathematical review: one-profile Picard eigenvalue existence

**Date:** 2026-10-14  
**Target:** \(\Gamma=\mathrm{PSL}_2(\mathbb Z[i])\), level one, odd quarter-turn sector, \(0<\lambda<45\).  
**Status:** **Computer-assisted existence theorem conditional only on the explicitly recorded and verified MPFR/analytic inputs**. This is a new independent interval replay and proof review, not a proof of the first *full* Picard eigenvalue, a width-0.1 spectral enclosure, or any \(\Gamma_0(2+i)\) certificate. No existing `rung4_certified` flag is promoted.

## The mathematical proof

Use upper-half-space coordinates \((x,t,y)\), positive hyperbolic Laplacian
\(\Delta=-y^2(\partial_x^2+\partial_t^2+\partial_y^2)+y\partial_y\), and the classical Picard Ford cell

\[
\mathcal F=\{(x,t,y):-1/2\le x\le1/2,\quad0\le t\le1/2,\quad x^2+t^2+y^2\ge1\}.
\]

This is a full volume fundamental cell, with sides paired by integer translation \(T_1\), \(R:z\mapsto-z\), \(T_iR:z\mapsto i-z\), and \(S:z\mapsto-1/z\) (interpreted in the three-dimensional quaternionic action). The height floor is \(y_0=\sqrt{1-x^2-t^2}\), with \(y_0\ge1/\sqrt2\). Reference checks: J. Voight, *Quaternion Algebras*, chapter 36.6 (Picard modular group), and GAP/HAP tutorial, §13.9 (Picard group).

Set \(r=6633/1000\), \(\lambda_0=1+r^2=44996689/1000000=44.996689\), \(X=2\pi x\), \(T=2\pi t\), and

\[
f(x,t,y)=2yK_{ir}(2\pi y)(\cos X-\cos T).
\]

**Admissibility / exact quotient form domain.** The function is smooth on the closed core and decays exponentially with its first derivatives as \(y\to\infty\). At the \(x=\pm1/2\) faces, translation matches its values. At \(t=0\) the \(R\)-pairing sends \(x\mapsto-x\), and at \(t=1/2\) the \(T_iR\)-pairing sends \(x\mapsto-x\); \(\cos(2\pi x)\) is even. On the hemispherical floor, the full \(S\) action is \((x,t,y)\mapsto(-x,t,y)\) because \(x^2+t^2+y^2=1\), and the value agrees there as well. These traces are enough to define \(f\in H^1(\Gamma\backslash\mathbb H^3)\); normal-derivative equality at the floor is *not* necessary for a variational test. The finite orbifold stabilizers only impose the same trace identities.

The Gaussian quarter-turn \(Q(x,t,y)=(-t,x,y)\) normalizes \(\Gamma\) and induces an involution on the quotient, as \(Q^2=R\in\Gamma\). On either half of the reference cell, applying \(R\) if necessary to return to the same half-strip, the profile satisfies \(f\circ Q=-f\). Consequently the cusp torus average of any odd function vanishes. Fourier Poincaré, for \(y\ge Y\ge1\), yields

\[
\int_{y\ge Y}|u|^2y^{-3}\,dx\,dt\,dy\le\frac1{4\pi^2Y^2}\int_{y\ge Y}(|u_x|^2+|u_t|^2)y^{-1}\,dx\,dt\,dy.
\]

Together with compact-core Rellich this proves compactness of the odd-sector \(H^1\hookrightarrow L^2\) embedding. Thus its Rayleigh infimum is attained by an odd cuspidal eigenfunction (and 0 is excluded because constants are even).

Every term of \(f\) solves \(\Delta f=\lambda_0f\) in the interior. Green's identity, with metric volume \(y^{-3}\,dx\,dt\,dy\), gives

\[
\mathcal E(f)-45\|f\|_2^2=(\lambda_0-45)\|f\|_2^2+B_{\rm floor}.
\]

Vertical fluxes vanish individually: \(f_x=0\) on \(x=\pm1/2\), and \(f_t=0\) on \(t=0,1/2\). The cusp-top term vanishes by exponential decay. The floor has *outward* unit Euclidean normal \((-x,-t,-y)\), area element \(dx\,dt/y\), and weighted boundary contribution

\[
B_{\rm floor}=-\int_{-1/2}^{1/2}\int_0^{1/2}\frac{f(x,t,y_0)\,[xf_x+tf_t+y_0f_y]}{y_0^2}\,dt\,dx.
\]

Setting \(U=\cos(2\pi x)-\cos(2\pi t)\), \(DU=2\pi[-x\sin(2\pi x)+t\sin(2\pi t)]\), \(z=2\pi y_0\), and \(K=K_{ir}(z)\), this simplifies **exactly** to

\[
\boxed{B_{\rm floor}=-4\int_{-1/2}^{1/2}\int_0^{1/2}[K^2(U^2+U\,DU)+zKK' U^2]dt\,dx.}
\]

Since \(\lambda_0-45=-3311/10^6<0\) and \(f\not\equiv0\), a *strictly negative floor flux* proves \(\inf\sigma(\Delta|_{Q\text{-odd}})<45\), and odd compactness makes this a genuine cuspidal eigenvalue.

## Independent directed-rounding computation

The executable `verify_floor_mpfr.c` uses the system `libmpfr.so.6` directly with 256-bit directed rounding. It neither imports `mpmath.iv` nor calls the previous interval script. Independent decimal root proposals (`gl_n*.txt`) are **not trusted**: their Legendre polynomial signs are revalidated at both bracket endpoints with outward MPFR intervals; ordered disjoint roots and strictly positive Gauss weights are rederived from the exact Legendre derivative formula. The sum-of-weights identity is enclosed. Real Bessel integrals \(K_{ir}(z)=\int_0^\infty e^{-z\cosh u}\cos(ru)du\) and \(K'_{ir}(z)=-\int_0^\infty\cosh u\,e^{-z\cosh u}\cos(ru)du\) are then evaluated with independent directed interval elementary functions.

Every Bessel value and derivative is widened by \(10^{-12}\); the final floor Gauss sum is widened by \(10^{-12}\). Both are *analytic, non-sampling error* allowances, as proved below.

| Spatial Gauss / Bessel Gauss | Independent MPFR final floor-flux enclosure | Sign |
|---|---|---|
| 56×56 / 80 | \([-3.467927098879\cdot10^{-11}, -3.267794988903\cdot10^{-11}]\) | Strict negative |
| 64×64 / 88 | \([-3.521057971848\cdot10^{-11}, -3.214662750934\cdot10^{-11}]\) | Strict negative |

The 64/88 evaluation has wider interval dependency than the first run, but both independently enclose a negative real number. An intentional mutation reversing the floor-normal sign gives a strictly positive interval, and the verifier exits nonzero (`STRICT NEGATIVE NO`). This checks that the asserted sign is load-bearing.

## Audited analytic error bounds

**Bessel integral (u variable).** On the Bernstein ellipse \(\rho_u=13/10\) of the affine image of \([-1,1]\) to \([0,6]\), the imaginary part satisfies \(|\Im u|\le\tfrac32(13/10-10/13)<4/5\), and \(|\Re u|<611/100\). For real \(z=2\pi y_0\ge\sqrt2\pi>4\), \(r=6633/1000\), we checked in independently directed MPFR arithmetic

\[
|\cos(ru)|<102,\quad|\cosh u|<226,\quad\Re(z\cosh u)>8/3,\quad e^{-8/3}<1/14.
\]

Thus \(\lvert e^{-z\cosh u}\cos(ru)\rvert<100\) and \(\lvert\cosh u\,e^{-z\cosh u}\cos(ru)\rvert<2000\). The full \([0,6]\) n-point Gauss errors are bounded by the Chebyshev approximation lemma below:
\(\epsilon_{u,n}\le160000(10/13)^{2n-1}\). For \(n\ge80\) this is < \(1.23\times10^{-13}\). Truncation beyond 6 is less than \(10^{-100}\): substitute \(s=\cosh u\), note \(\cosh6>200\), and bound both derivative and value tails by \(1.01\int_{200}^\infty e^{-4s}ds\), which is much less than \(10^{-100}\). Hence widening each K/K' enclosure by \(10^{-12}\) is sound.

**Spatial floor (x,t variables).** Use product Bernstein ellipses with \(\rho_f=7/5\) over \([-1/2,1/2]\) and \([0,1/2]\). The exact rational bounds \(|\Re x|\le37/70\), \(|\Im x|\le6/35\), \(|\Re t|\le18/35\), \(|\Im t|\le3/35\) yield
\(\Re(1-x^2-t^2)\ge2235/4900>1/4\). Therefore the principal branch of \(y_0=\sqrt{1-x^2-t^2}\) is analytic throughout the closed product ellipses, \(\Re(y_0)>1/2\), \(|y_0|<2\), and \(\Re(z)>\pi\).

For real \(u\ge0\), \(\cosh u\ge1+u^2/2\) and \(\cosh u\le e^{u^2/2}\). Hence the integral representation gives the following valid bounds **even for complex z with \(\Re z>\pi\)**:

\[
|K_{ir}(z)|\le e^{-\pi}\sqrt{\frac\pi{2\pi}}=\frac{e^{-\pi}}{\sqrt2}<\frac1{20},\qquad
|K'_{ir}(z)|\le e^{-\pi}\sqrt{\frac\pi{2(\pi-1)}}<\frac1{20}.
\]

The trig bounds \(\cosh(12\pi/35)<33/20\) and \(\cosh(6\pi/35)<6/5\) imply \(|U|<4\), \(|DU|<15\). Also \(|z|<4\pi<13\). Consequently the analytic *full floor integrand* has modulus bounded by
\(4[(76+208)/400]=2.84<100\), much less than the declared `M_floor=100`.

**Chebyshev/Gauss lemma.** If \(g\) is analytic on and inside a Bernstein ellipse E_ρ over a real interval of length L with \(|g|\le M\), its Chebyshev degree-\((2n-1)\) truncation has uniform error \(\le\frac{2M}{\rho-1}\rho^{-(2n-1)}\). Gauss exactness and positivity of the quadrature weights give an integral error \(\le\frac{4LM}{\rho-1}\rho^{-(2n-1)}\). Applying the one-dimensional estimate successively in x and t, using the positive weight sums and \(L_x=1,L_t=1/2\), gives the tensor upper bound \(\le1000(5/7)^{2n-1}\). We deliberately use \(10000(5/7)^{2n-1}\), which the separate MPFR constants audit verifies is <\(10^{-12}\) for \(n\ge56\).

`audit_constants_mpfr.c` independently verifies 15 strict interval inequalities, including all decisive ellipse, Bessel, and quadrature margins. The analytic assumptions \(\cosh u\le e^{u^2/2}\), Chebyshev approximation, MPFR's directed rounding semantics, and the Picard cell theorem are standard mathematical inputs and are laid out explicitly here rather than hidden in software.

## Reproducibility and negative control

On x86_64 Linux with libmpfr.so.6 and libgmp:

```sh
python generate_gl_brackets.py 56 64 80 88
cc -O2 -std=c11 verify_floor_mpfr.c -o verify_floor_mpfr -l:libmpfr.so.6 -lgmp -lm
cc -O2 -std=c11 audit_constants_mpfr.c -o audit_constants_mpfr -l:libmpfr.so.6 -lgmp -lm
./audit_constants_mpfr
./verify_floor_mpfr 56 80 gl_n56.txt gl_n80.txt
./verify_floor_mpfr 64 88 gl_n64.txt gl_n88.txt
```

The standalone MPFR verifier can be compiled against the official `mpfr.h`, if installed. A marked x86_64 LP64 `libmpfr.so.6` structure declaration is supplied as an emergency fallback when distribution runtime libraries exist without development headers. It must *not* be treated as portable across arbitrary MPFR ABIs.

**Independent versus correlated inputs.** The two MPFR grid sizes use the same independent MPFR verifier and the same analytic theorem; the *earlier* Python `mpmath.iv` implementation is a separate arithmetic backend and a separate programming implementation. The root proposal generator imports SciPy and mpmath, but the proof is independent of their answers because MPFR validates each proposed rational bracket and weight.

**Limitations.** This is a self-contained computer-assisted theorem candidate with two arithmetic backends and the explicit analytical argument above. It has not been externally peer-reviewed, and only the *existence below 45 in the level-one Q-odd sector* follows. No claim for level \((2+i)\), newspace, a precise 44.85 eigenvalue, or the first entire-spectrum eigenvalue is authorized. Preserve all other proof flags unchanged.
