# Analytic implementation links

The trace-formula identity itself remains the cited external input. This
note gives the local reductions needed to reproduce the numerical work
without private OCR extracts. Conventions are
h(r)=sinc^(2k)(delta*r) and g(x)=(2*pi)^-1 integral h(r) exp(-irx) dr.

## Polynomial derivative and admissibility

Convolving 2k uniform densities on [-delta,delta] gives

\[
g(x)=\frac{\sum_{j=0}^{2k}(-1)^j{2k\choose j}
 (x+2\delta(k-j))_+^{2k-1}}{(2k-1)!(2\delta)^{2k}}.
\]

For integer k>=2, differentiate twice at zero. Terms j>=k vanish, giving

\[
g''(0)=\frac{1}{8(2k-3)!\delta^3}
 \sum_{j=0}^{k-1}(-1)^j{2k\choose j}(k-j)^{2k-3}.
\]

Thus int h(r)r^2 dr=2*pi*(-g''(0)), and the identity term is
vol/(2*pi)*(-g''(0)). The rational coefficients for k=2,3,4 are
-1/4,-1/8,-1/12. The old four-point stencil was exact for the cubic
k=2 piece only; it overestimates -g''(0) at k=3 by 3/16.
k=1 is rejected because the identity integral diverges.
Arb's entire sinc evaluates the removable singularity at zero.
Support is [-2k*delta,2k*delta]; the runtime compares its upper endpoint
against the lower endpoint of the proved systole before omitting loxodromics.

## Field and scattering constants

For a negative fundamental discriminant D,
zeta_K(s)=zeta(s)L(s,chi_D), and
L(s,chi)=q^-s sum chi(a) zeta(s,a/q), q=|D|.
Cancellation of the Hurwitz residues gives
L(1,chi)=-(1/q) sum chi(a) psi(a/q).
L'(1,chi) is enclosed by Cauchy integration around s=1, using the exact
Arb angular endpoint 2*pi. Humbert's volume is
|D|^(3/2) zeta_K(2)/(4*pi^2).

For class number one, the lattice Epstein series is w*zeta_K(s), with
covolume V=sqrt(|D|)/2 and w the number of units. Its constant Laurent
coefficient is c0=w*(gamma*L(1)+L'(1)). Partial summation gives the lattice
Euler constant eta=(V/pi)*c0. The cusp rotation index is GG=w/2.
For d=2,7,11,19, units are +/-1, so GG=1 and the cuspidal-elliptic sector
vanishes. Non-cuspidal elliptics still need a complete inventory.

With the single-cusp determinant
phi_K(s)=(2*pi/sqrt(|D|))*zeta_K(s-1)/((s-1)*zeta_K(s)), define
Xi_K(s)=s(s-1)|D|^(s/2)(2*pi)^-s Gamma(s)zeta_K(s) and Q=Xi'/Xi.
The functional equation gives, on s=1+ir,
phi'/phi=-Q(1-ir)-Q(1+ir)+(1-ir)^-1+(1+ir)^-1.
Evenness of h then expresses PHIINT as
(2*pi)^-1 integral h(r)/(1+r^2) dr minus
(2*pi)^-1 integral h(r)Q(1+ir) dr.
Shift the latter to Re(s)=2 through the zero-free strip. There are no
crossed zeros; horizontal sides vanish for sinc^(2k), k>=2.
On the new line Q=zeta_K'/zeta_K+B_K, where
B_K(t)=(2+it)^-1+(1+it)^-1+C_K+psi(2+it),
C_K=log(|D|)/2-log(2*pi).

The absolutely convergent Euler product yields the remaining term

\[
\sum_{\mathfrak p}\sum_{m\ge1}
 \frac{\log N\mathfrak p}{(N\mathfrak p)^m}
 g(m\log N\mathfrak p).
\]

Only norms within exp(support) survive. Ramified, split and inert rational
primes give respectively (norm,multiplicity)=(p,1),(p,2),(p^2,1).
The numerator for a power is log(Np), not log(Np^m).
This reproduces the Picard log(2)/2*g(log(2)) and Eisenstein zero term.
The full contour argument and cusp formula are also archived in
`old_RIGOR_GAPS.md`; its elliptic normalization is subject to the issue note.

## Elementary digamma and tail bounds

For r>=0 the convergent series gives
Re psi(1+ir)=-gamma+sum r^2/[n(n^2+r^2)] and
Im psi(1+ir)=sum r/(n^2+r^2), n>=1.
The real summand decreases with n, so its sum is at most
1+(1/2)log(1+r^2); the imaginary sum is at most the integral from zero
to infinity, pi/2. The real part is bounded below by -gamma. Consequently
|psi(1+ir)|<=log(1+r)+1+pi/2.
This establishes the needed majorant without a complex asymptotic
remainder convention.

On the real line |h(r)|<=(delta*r)^(-2k). Integration by parts gives,
for p=2k, the half-line bound

\[
\int_R^\infty r^{-p}(\log(1+r)+c)\,dr
\le R^{1-p}\left(\frac{\log(1+R)+c}{p-1}
 +\frac1{(p-1)^2}\right).
\]

The elementary h/(1+r^2) tail is bounded by
delta^-p R^(-p-1)/(p+1), and is nonnegative.
For the shifted integrand,
|h(t-i)|<=(cosh(delta)/(delta*t))^p.
The digamma recurrence gives
|B_K(t)|<=log(1+t)+1+pi/2+3/t+|C_K|.
For R>=7 this is at most log(2+t)+3+|C_K|. The implementation takes
an outward bound at least max(|C_K|,log(2*pi)), preserving the historical
bound and making its discriminant dependence explicit.
Each main integral covers [-R,R], so each of these tails has a leading
factor of two. Both PSI and the shifted tail use symmetric interval hulls;
the nonnegative elementary tail uses [0,tail].

Arb quadrature uses ball endpoints, and the cuspidal integral is split at
every B-spline knot with a fixed polynomial integrand on each piece.
No new-field spectral conclusion follows from these analytic reductions
until the missing elliptic data are proved.
