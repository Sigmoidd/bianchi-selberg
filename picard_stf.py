"""
Certified Selberg trace formula engine for the Picard group PSL(2, Z[i]).

Goal: rigorously bound the exceptional spectrum (Laplace eigenvalues
0 < lambda < 1) of the Picard orbifold  M = PSL(2,Z[i]) \\ H^3.

Method
------
Selberg trace formula (Elstrodt-Grunewald-Mennicke Thm 6.5.1 as corrected and
restated by J. Friedman, arXiv:math/0612807 Thm 4.1.1, and by
Balkanova-Chatzakos-Cherubini-Frolenkov-Laaksonen, arXiv:1712.00880 Thm 2.2),
specialized to Gamma = PSL(2, Z[i]), trivial character, one cusp:

  sum_j h(r_j)  =  I + NCE + LOX + SCATT + PHIINT + CE + PAR

with lambda_j = 1 + r_j^2 (r_j in R or r_j in i(0,1]; lambda_0 = 0 <-> r_0 = i)
and g(x) = (1/2pi) int h(r) e^{-irx} dr:

  I      = vol(M)/(4 pi^2) * int_R h(r) r^2 dr,          vol(M) = 2 zeta_K(2)/pi^2
  NCE    = g(0) * log N(T_0) / (4 |E(R)| sin^2(pi k/m)) summed over non-cuspidal
           elliptic classes.  For Z[i]: one class, m = |E| = 3, sin^2 = 3/4,
           N(T_0) = 7 + 4 sqrt(3)   =>   NCE = g(0) * log(7+4sqrt3) / 9
  LOX    = sum over loxodromic classes of g(log N(T)) * (weights)  == 0 here,
           because supp g in [-L, L], L <= l_0 = log((3+sqrt5)/2) = min length
  SCATT  = - tr S(0) h(0)/4 = + h(0)/4          (phi(1) = -1 for Z[i])
  PHIINT = (1/4pi) int_R h(r) (phi'/phi)(1+ir) dr,
           phi(s) = pi * zeta_K(s-1) / ((s-1) zeta_K(s))
  CE     = (5/16) log2 * g(0) + (1/4) int_0^inf g(x) tanh(x/2) dx
           [4 cuspidal-elliptic classes S_w, w in {0,1,i,1+i}: |C| = 4,
            |1-eps^2|^2 = 4, |c_i| = 1, 2, 2, sqrt2]
  PAR    = h(0)/8 + (g(0)/2)(eta/2 - gammaE)
           - (1/4pi) int_R h(r) psi(1+ir) dr
           eta = 2 gammaE + 2 log2 + 3 log pi - 4 log Gamma(1/4)
           (lattice Euler constant of Z[i], via beta'(1) closed form)

PHIINT is evaluated by shifting contours off the critical lines (valid since
xi~(s) = s(s-1) pi^{-s} Gamma(s) zeta_K(s) is entire and zero-free in
-1 <= Re <= 0 and 1 <= Re <= 2), giving the exactly equivalent expression

  PHIINT = (1/2pi) int_R h(r)/(1+r^2) dr
         + (log2 / 2) * g(log 2)
         - (1/2pi) int_R h(t-i) [ 1/(2+it) + 1/(1+it) - log pi + psi(2+it) ] dt

where the prime sum over n >= 3 vanishes because g(log n) = 0 (log 3 > L).

Test functions:  h(r) = (sin(delta r)/(delta r))^{2k}, so that
  g = (2k)-fold convolution of (1/(2 delta)) 1_[-delta,delta]  (cardinal B-spline),
  supp g = [-2k delta, 2k delta],  h >= 0 on R,  h(i sigma) > 1 for sigma in (0,1].

Criterion: with B := RHS - h(i)  (subtracting the constant eigenfunction),
  every exceptional eigenvalue lambda = 1 - sigma^2 contributes
  h(i sigma) = (sinh(delta sigma)/(delta sigma))^{2k} > 1,
  all other discrete eigenvalues contribute >= 0.  Hence:
    B < 1                =>  NO exceptional eigenvalues: lambda_1 >= 1.
    else                 =>  lambda_1 >= 1 - sigma*^2, h(i sigma*) = B.

All arithmetic uses Arb certified balls via python-flint; integrals use
acb_calc adaptive certified quadrature; all infinite tails have explicit
closed-form ball bounds.
"""

from flint import arb, acb, ctx
import math

# ----------------------------------------------------------------------------
# precision
# ----------------------------------------------------------------------------
def set_prec(bits=192):
    ctx.prec = bits

TWO = arb(2)

# ----------------------------------------------------------------------------
# basic special values
# ----------------------------------------------------------------------------
def const_pi():     return arb.pi()
def const_gamma():  return arb.const_euler()

def zeta_hurwitz(s, a):
    """Certified Hurwitz zeta zeta(s, a) for acb s, arb/acb a."""
    return acb(s).zeta(acb(a))

def L_chi4(s):
    """Dirichlet L(s, chi_{-4}) = 4^{-s} (zeta(s,1/4) - zeta(s,3/4)). Certified."""
    s = acb(s)
    return acb(4)**(-s) * (s.zeta(acb(1)/4) - s.zeta(acb(3)/4))

def zeta_K(s):
    """Dedekind zeta of Q(i): zeta(s) * L(s, chi_{-4}). Certified."""
    s = acb(s)
    return s.zeta() * L_chi4(s)

def volume_picard():
    """vol(M) = 2 zeta_K(2) / pi^2 (= Catalan/3 * 2 ... = 0.30532...)."""
    return (2 * zeta_K(2) / const_pi()**2).real

def eta_lattice_Zi():
    """
    Lattice Euler constant eta for Lambda = Z[i]:
      sum_{0<|mu|^2<=x} 1/|mu|^2 = pi (log x + eta) + O(x^{-1/2})
    Closed form: eta = 2 gammaE + 2 log 2 + 3 log pi - 4 log Gamma(1/4),
    derived from sum r_2(n)/n^s = 4 zeta(s) L(s,chi4) and
    beta'(1) = (pi/4)(gammaE + 2 log 2 + 3 log pi - 4 log Gamma(1/4)).
    """
    g = const_gamma()
    return 2*g + 2*arb(2).log() + 3*const_pi().log() - 4*(arb(1)/4).gamma().log()

# ----------------------------------------------------------------------------
# test function h(r) = sinc(delta r)^(2k) and its B-spline transform g
# ----------------------------------------------------------------------------
def _sinc_acb(z):
    """Certified sinc(z) = sin(z)/z, holomorphic, safe for balls containing 0."""
    z = acb(z)
    # magnitude test: use midpoint estimate; both branches are enclosures of the
    # same entire function, so branch choice only affects tightness.
    m = abs(complex(z))
    if m > 0.5:
        return z.sin() / z
    # Taylor series sum_{m=0}^{M-1} (-1)^m z^(2m)/(2m+1)! with rigorous tail:
    # tail <= |z|^(2M)/(2M+1)! * 1/(1 - |z|^2/((2M+2)(2M+3)))  for |z|<=1
    M = 10
    w = z*z
    total = acb(0)
    term = acb(1)   # z^(2m)/(2m+1)! at m=0
    fact = 1
    for mm in range(M):
        if mm > 0:
            term = term * w / ((2*mm)*(2*mm+1))
        total += term if mm % 2 == 0 else -term
    zab = abs(z)
    # |z|^{2M} / (2M+1)! * geometric guard
    tail = zab**(2*M)
    f = arb(1)
    for j in range(2, 2*M+2):
        f *= j
    tail = tail / f
    tail = tail / (1 - arb(1)/((2*M+2)*(2*M+3)))
    rad = arb(0).union(tail).union(-tail)
    return total + acb(rad, rad)  # pad real & imag by tail

def make_h(k, delta):
    """Return h(z) certified for acb z."""
    d = arb(delta)
    def h(z):
        return _sinc_acb(acb(z) * d) ** (2*k)
    return h

def h_at_i_sigma(k, delta, sigma):
    """h(i sigma) = (sinh(delta sigma)/(delta sigma))^{2k}, arb certified."""
    x = arb(delta) * arb(sigma)
    return (x.sinh() / x) ** (2*k)

def binomial(n, j):
    return math.comb(n, j)

def gfun_bspline(k, delta, x):
    """
    Exact g(x) = (1/(2 delta)^{2k}) B_{2k}(x), B-spline = 2k-fold convolution
    of 1_{[-delta,delta]}:
      g(x) = 1/((2d)^{2k} (2k-1)!) sum_{j=0}^{2k} (-1)^j C(2k,j) (x+(k-j)2d)_+^{2k-1}
    x: arb (may be a ball). Returns certified arb.
    """
    d = arb(delta)
    x = arb(x)
    n = 2*k
    fact = 1
    for j in range(2, n):
        fact *= j            # (2k-1)!
    total = arb(0)
    for j in range(n+1):
        y = x + (k - j) * 2 * d
        # (y)_+^{n-1} with ball-safe positive part
        lo, hi = y.lower(), y.upper()
        if hi <= 0:
            continue
        if lo >= 0:
            p = y ** (n-1)
        else:
            # ball straddles 0: hull of [0, hi^(n-1)]
            hb = arb(hi) ** (n-1)
            p = (hb/2).union(arb(0))
            p = arb(0).union(hb)
        t = binomial(n, j) * p
        total = total + t if j % 2 == 0 else total - t
    return total / (fact * (2*d)**n)

def g0_exact(k, delta):
    return gfun_bspline(k, delta, arb(0))

# ----------------------------------------------------------------------------
# certified integration helpers
# ----------------------------------------------------------------------------
def integral(f, a, b, panel=4.0, abs_tol=1e-20):
    """Certified integral of holomorphic f over [a, b] via acb_calc,
    split into panels for tighter adaptive behavior on oscillatory tails."""
    a = float(a); b = float(b)
    n = max(1, int((b - a) / panel + 0.5))
    step = (b - a) / n
    total = acb(0)
    for i in range(n):
        total += acb.integral(lambda z, _: f(z), a + i*step, a + (i+1)*step,
                              abs_tol=abs_tol, eval_limit=10**7)
    return total

# ----------------------------------------------------------------------------
# the individual trace formula terms   (all return arb enclosures)
# ----------------------------------------------------------------------------
def I_term(k, delta, R):
    """(vol/4pi^2) int_R h(r) r^2 dr  with certified tail beyond R.
       tail: |h| <= (delta r)^{-2k} => int_R^inf r^2 (delta r)^{-2k} dr
             = delta^{-2k} R^{3-2k} / (2k-3)."""
    h = make_h(k, delta)
    vol = volume_picard()
    pi = const_pi()
    main = 2 * integral(lambda z: h(z) * z*z, 0, R).real
    d = arb(delta)
    tail = 2 * d**(-2*k) * arb(R)**(3-2*k) / (2*k-3)
    val = (vol / (4*pi*pi)) * (main + arb(0).union(tail))
    return val

def NCE_term(k, delta):
    """One order-3 class: g(0) * log(7+4 sqrt 3) / 9."""
    N0 = 7 + 4*arb(3).sqrt()
    return g0_exact(k, delta) * N0.log() / 9

def CE_terms(k, delta):
    """(5/16) log2 g(0)  +  (1/4) int_0^L g(x) tanh(x/2) dx  (piecewise-analytic)."""
    g0 = g0_exact(k, delta)
    part1 = arb(5)/16 * arb(2).log() * g0
    # integral over B-spline pieces [2 d m, 2 d (m+1)], m = 0..k-1;
    # on each piece g is a fixed polynomial: active terms j <= k+m.
    d = arb(delta)
    n = 2*k
    fact = 1
    for j in range(2, n):
        fact *= j
    total = arb(0)
    for m in range(k):
        a_, b_ = 2*delta*m, 2*delta*(m+1)
        def piece(z, m=m):
            z = acb(z)
            s = acb(0)
            for j in range(0, k+m+1):
                y = z + (k - j)*2*d
                t = binomial(n, j) * y**(n-1)
                s = s + t if j % 2 == 0 else s - t
            gpoly = s / (fact * (2*d)**n)
            return gpoly * (z/2).tanh()
        total += integral(piece, a_, b_).real
    return part1 + total/4

def SCATT_PAR_h0():
    """- tr S(0) h(0)/4 + h(0)/8 = h(0)/4 + h(0)/8 = 3/8   (h(0)=1)."""
    return arb(3)/8

def PAR_g0(k, delta):
    """(g(0)/2) (eta/2 - gammaE)."""
    return g0_exact(k, delta)/2 * (eta_lattice_Zi()/2 - const_gamma())

def PSI_term(k, delta, R):
    """- (1/4pi) int_R h(r) psi(1+ir) dr  = -(1/2pi) Re int_0^R h psi(1+ir) dr - tail.
       tail bound: |psi(1+ir)| <= log(1+r) + 1 + pi/2 for r real."""
    h = make_h(k, delta)
    pi = const_pi()
    I = integral(lambda z: h(z) * (1 + acb(0,1)*z).digamma(), 0, R)
    main = -(1/(2*pi)) * I.real
    d = arb(delta)
    Ra = arb(R)
    # int_R^inf (delta r)^{-2k} (log(1+r) + c) dr,  c = 1 + pi/2
    c = 1 + pi/2
    tail_int = d**(-2*k) * ( ((1+Ra).log() + c) * Ra**(1-2*k)/(2*k-1)
                             + Ra**(1-2*k)/((2*k-1)**2) )
    tail = (1/(2*pi)) * tail_int
    return main + arb(0).union(tail).union(-tail)

def PHI_term(k, delta, R):
    """
    PHIINT = (1/2pi) int_R h/(1+r^2) dr  + (log2/2) g(log2)
             - (1/2pi) int_R h(t-i) [1/(2+it) + 1/(1+it) - log pi + psi(2+it)] dt
    Each doubled from [0,R] by symmetry (second integrand: conjugate symmetry).
    Tails certified via |h(t-i)| <= (cosh(delta)/(delta t))^{2k} etc.
    """
    h = make_h(k, delta)
    pi = const_pi()
    d = arb(delta)
    Ra = arb(R)

    # elementary part
    e_main = 2 * integral(lambda z: h(z) / (1 + z*z), 0, R).real
    e_tail = 2 * d**(-2*k) * Ra**(-1-2*k) / (2*k+1)
    part_elem = (1/(2*pi)) * (e_main + arb(0).union(e_tail))

    # prime part: only n = 2 survives (log 3 > supp g)
    part_prime = arb(2).log()/2 * gfun_bspline(k, delta, arb(2).log())

    # shifted line integral
    logpi = const_pi().log()
    def F(t):
        t = acb(t)
        it = acb(0,1)*t
        return h(t - acb(0,1)) * ( 1/(2+it) + 1/(1+it) - logpi + (2+it).digamma() )
    s_main = 2 * integral(F, 0, R).real
    # tail: |h(t-i)| <= (cosh d / (d t))^{2k};  |bracket| <= log(1+t) + c2
    c2 = arb(3)/2 + logpi + 2 + pi/2   # 1/2 + 1 <= 3/2; |psi(2+it)|<=log(1+t)+2+pi/2... consolidated
    ch = d.cosh()
    s_tail = 2 * (ch/d)**(2*k) * ( ((1+Ra).log() + c2) * Ra**(1-2*k)/(2*k-1)
                                   + Ra**(1-2*k)/((2*k-1)**2) )
    part_shift = -(1/(2*pi)) * (s_main + arb(0).union(s_tail).union(-s_tail))

    return part_elem + part_prime + part_shift

# ----------------------------------------------------------------------------
# assembly
# ----------------------------------------------------------------------------
def compute_B(k, delta, R=80, prec=192, verbose=True):
    """
    Compute certified enclosure of
      B = RHS(trace formula) - h(i)
    which upper-bounds sum over all discrete eigenvalues except lambda_0 = 0.
    Requires 2 k delta <= l_0 = log((3+sqrt5)/2) (checked).
    """
    set_prec(prec)
    l0 = ((3 + arb(5).sqrt())/2).log()
    L = 2*k*arb(delta)
    assert L.upper() <= l0.lower(), f"support 2k*delta={float(L)} exceeds systole {float(l0)}"

    I   = I_term(k, delta, R)
    NCE = NCE_term(k, delta)
    CE  = CE_terms(k, delta)
    S38 = SCATT_PAR_h0()
    Pg  = PAR_g0(k, delta)
    PSI = PSI_term(k, delta, R)
    PHI = PHI_term(k, delta, R)
    hi  = h_at_i_sigma(k, delta, 1)

    rhs = I + NCE + CE + S38 + Pg + PSI + PHI
    B = rhs - hi
    if verbose:
        for name, v in [("I", I), ("NCE", NCE), ("CE", CE), ("SCATT+PARh0", S38),
                        ("PARg0", Pg), ("PSI", PSI), ("PHI", PHI),
                        ("h(i) [const eigenfn]", hi), ("RHS total", rhs), ("B", B)]:
            print(f"  {name:22s} = {v}")
    return B

def sigma_star(k, delta, B_upper):
    """Smallest sigma with h(i sigma) >= B_upper: exceptional eigenvalues with
       sigma > sigma* are excluded; returns certified upper bound for sigma*."""
    if B_upper < 1:
        return 0.0
    lo, hi_ = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi_)/2
        v = h_at_i_sigma(k, delta, mid)
        if v.lower() > B_upper:
            hi_ = mid
        else:
            lo = mid
    return hi_

if __name__ == "__main__":
    import sys
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    delta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.24
    R = float(sys.argv[3]) if len(sys.argv) > 3 else 80
    print(f"PSL(2,Z[i]) level 1, h = sinc^{2*k}(delta r), delta={delta}, "
          f"supp g = [{-2*k*delta:.4f},{2*k*delta:.4f}], R={R}")
    B = compute_B(k, delta, R)
    bu = B.upper()
    print()
    if bu < 1:
        print(f"  B < 1 CERTIFIED  =>  NO exceptional eigenvalues:")
        print(f"  lambda_1(PSL(2,Z[i])) >= 1   [Selberg conjecture, level 1]")
    else:
        ss = sigma_star(k, delta, bu)
        lam = 1 - ss*ss
        print(f"  B >= 1: certified exclusion only for sigma > {ss:.6f}")
        print(f"  lambda_1 >= {lam:.6f}   (Blomer-Brumley general bound: 975/1024 = {975/1024:.6f})")
