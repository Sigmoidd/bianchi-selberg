"""
Exploratory feasibility screen (mpmath; NOT a spectral certificate).

For K = Q(sqrt(-d)), class number 1, compute every term of
    B = I + NCE + CE + Ch0 + PARg0 + PSI + PHIINT - h(i)
that does NOT depend on the elliptic inventory, using the same contour-shift
formulas as bianchi_omega_arb.py generalized to arbitrary fundamental
discriminant D and cusp index GG = |O_K^*|/2.

  B_mech := I + Ch0 + PARg0 + PSI + PHIINT - h(i)          (NCE, CE dropped)

Since NCE >= 0 and CE >= 0 (g >= 0, positive kernels), B >= B_mech.  Hence
B_mech >= 1  =>  this test function CANNOT certify lambda_1 >= 1, whatever the
elliptic inventory is.  Validated on d=1, 3 by adding their known NCE+CE.
"""
import sys
import mpmath as mp
from math import comb
mp.mp.dps = 25

# ---------------------------------------------------------------- field data
FIELDS = {1: -4, 2: -8, 3: -3, 7: -7, 11: -11, 19: -19, 43: -43, 67: -67, 163: -163}

def chi(D, n):
    """Kronecker symbol (D/n) for the fundamental discriminants used here."""
    n = int(n)
    if D == -4:
        return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
    if D == -8:
        r = n % 8
        return {1: 1, 3: 1, 5: -1, 7: -1}.get(r, 0)
    p = -D                                   # odd prime p = 3 mod 4
    r = n % p
    if r == 0:
        return 0
    return 1 if pow(r, (p - 1) // 2, p) == 1 else -1

def L(s, D):
    q = -D
    return mp.mpf(q) ** (-s) * sum(chi(D, a) * mp.zeta(s, mp.mpf(a) / q)
                                   for a in range(1, q + 1) if chi(D, a))

def zetaK(s, D):
    return mp.zeta(s) * L(s, D)

def logder_zetaK(s, D):
    return mp.diff(lambda x: mp.log(zetaK(x, D)), s)

def field(d):
    D = FIELDS[d]
    w = {1: 4, 3: 6}.get(d, 2)
    V = mp.sqrt(-D) / 2
    vol = (-D) ** mp.mpf(1.5) * zetaK(2, D) / (4 * mp.pi ** 2)
    L1 = -(mp.mpf(1) / -D) * sum(chi(D, a) * mp.digamma(mp.mpf(a) / -D)
                                 for a in range(1, -D + 1) if chi(D, a))
    Lp1 = mp.diff(lambda x: L(x, D), 1)
    eta = (V / mp.pi) * w * (mp.euler * L1 + Lp1)
    return dict(d=d, D=D, w=w, GG=w // 2, V=V, vol=vol, eta=eta)

# ---------------------------------------------------------------- systole
def systole(d, search=8):
    """min over tau in O_K, tau not in [-2,2], of 2 log|lambda|, lambda+1/lambda=tau.
    Every tau is a trace in SL2(O_K) via (0 -1; 1 tau), so this is the exact systole."""
    D = FIELDS[d]
    best, arg = None, None
    # O_K = Z[g], g = sqrt(-d) (D=-4d) or (1+sqrt(-d))/2 (D=-d)
    for a in range(-search, search + 1):
        for b in range(-search, search + 1):
            if D % 4 == 0:
                tau = mp.mpc(a, b * mp.sqrt(d))
            else:
                tau = mp.mpc(a + mp.mpf(b) / 2, b * mp.sqrt(d) / 2)
            if abs(tau.imag) < 1e-20 and abs(tau.real) <= 2 + 1e-20:
                continue
            disc = mp.sqrt(tau * tau - 4)
            lam = max(abs((tau + disc) / 2), abs((tau - disc) / 2))
            ell = 2 * mp.log(lam)
            if best is None or ell < best:
                best, arg = ell, tau
    return best, arg

# ---------------------------------------------------------------- B-spline machinery
def gfun(x, k, dl):
    n = 2 * k
    fact = mp.factorial(n - 1)
    tot = mp.mpf(0)
    for j in range(n + 1):
        y = x + (k - j) * 2 * dl
        if y > 0:
            tot += (-1) ** j * comb(n, j) * y ** (n - 1)
    return tot / (fact * (2 * dl) ** n)

def gpp0(k, dl):
    """exact g''(0+) from the polynomial piece (valid for k>=2)."""
    n = 2 * k
    fact = mp.factorial(n - 1)
    tot = mp.mpf(0)
    for j in range(n + 1):
        y = (k - j) * 2 * dl
        if y > 0:
            tot += (-1) ** j * comb(n, j) * (n - 1) * (n - 2) * y ** (n - 3)
    return tot / (fact * (2 * dl) ** n)

def h(r, k, dl):
    x = dl * r
    return mp.mpf(1) if abs(x) < 1e-15 else (mp.sin(x) / x) ** (2 * k)

def prime_term(F, k, dl):
    """sum over ideals a, N(a)<=e^supp:  Lambda_K(a)/N(a) * g(log N(a))."""
    D, supp = F["D"], 2 * k * dl
    tot = mp.mpf(0)
    maxN = int(mp.floor(mp.e ** supp)) + 1
    from math import isqrt
    def primerange(start, stop):
        return (p for p in range(start, stop) if all(p % q for q in range(2, isqrt(p)+1)))
    for p in primerange(2, maxN + 1):
        # splitting type of p in K
        if D % p == 0:
            kind = "ram"
        else:
            c = chi(D, p) if p != 2 else (1 if D % 8 == 1 else -1)
            if p == 2:
                # D = 1 mod 8 splits, 5 mod 8 inert (D odd, D=-7 -> split)
                c = 1 if D % 8 == 1 else -1
            kind = "split" if c == 1 else "inert"
        if kind == "ram":
            primes = [(p, 1)]                       # (norm, multiplicity of such primes)
        elif kind == "split":
            primes = [(p, 2)]
        else:
            primes = [(p * p, 1)]
        for (Np, mult) in primes:
            m = 1
            while Np ** m <= mp.e ** supp:
                N = Np ** m
                tot += mult * mp.log(Np) / N * gfun(mp.log(N), k, dl)
                m += 1
    return tot

# ---------------------------------------------------------------- the terms
def mech_terms(F, k, frac=0.999, R=60):
    ell0, _ = systole(F["d"])
    dl = ell0 * frac / (2 * k)
    GG = F["GG"]
    g0 = gfun(mp.mpf(0), k, dl)
    I = F["vol"] / (2 * mp.pi) * (-gpp0(k, dl))
    Ch0 = mp.mpf(1) / 4 + mp.mpf(1) / (4 * GG)
    PARg0 = g0 * (F["eta"] / 2 - mp.euler) / GG
    PSI = -(mp.mpf(1) / GG) / (2 * mp.pi) * 2 * mp.quad(
        lambda r: h(r, k, dl) * mp.re(mp.digamma(1 + 1j * r)), mp.linspace(0, R, 31))
    CK = mp.log(-F["D"]) / 2 - mp.log(2 * mp.pi)
    part_elem = (1 / (2 * mp.pi)) * 2 * mp.quad(lambda r: h(r, k, dl) / (1 + r * r),
                                                 mp.linspace(0, R, 31))
    part_prime = prime_term(F, k, dl)
    def shift(t):
        it = 1j * t
        z = t - 1j
        x = dl * z
        hz = (mp.sin(x) / x) ** (2 * k)
        return hz * (1 / (2 + it) + 1 / (1 + it) + CK + mp.digamma(2 + it))
    part_shift = -(1 / (2 * mp.pi)) * 2 * mp.re(mp.quad(shift, mp.linspace(0, R, 31)))
    PHIINT = part_elem + part_prime + part_shift
    hi = (mp.sinh(dl) / dl) ** (2 * k)
    return dict(ell0=ell0, delta=dl, supp=2 * k * dl, g0=g0, I=I, Ch0=Ch0, PARg0=PARg0,
                PSI=PSI, PHIINT=PHIINT, hi=hi,
                Bmech=I + Ch0 + PARg0 + PSI + PHIINT - hi)

# known elliptic data for the two validated fields
KNOWN = {
    1: lambda g0, k, dl: g0 * mp.log(7 + 4 * mp.sqrt(3)) / 9
        + mp.mpf(5) / 16 * mp.log(2) * g0
        + mp.mpf(1) / 4 * mp.quad(lambda x: gfun(x, k, dl) * mp.sinh(x) / (mp.cosh(x) + 1),
                                  [0, 2 * k * dl]),
    3: lambda g0, k, dl: g0 * mp.log(7 + 4 * mp.sqrt(3)) / 8
        + mp.mpf(2) / 9 * mp.log(3) * g0
        + mp.mpf(1) / 3 * mp.quad(lambda x: gfun(x, k, dl) * mp.sinh(x) / (mp.cosh(x) + mp.mpf(1) / 2),
                                  [0, 2 * k * dl]),
}

if __name__ == "__main__":
    ds = [int(a) for a in sys.argv[1:]] or [1, 3, 2, 7, 11, 19]
    for d in ds:
        F = field(d)
        print(f"\n=== d={d}  D={F['D']}  w={F['w']}  GG={F['GG']}  vol={mp.nstr(F['vol'],8)}  eta={mp.nstr(F['eta'],8)}")
        ell0, tau = systole(d)
        print(f"    systole l0={mp.nstr(ell0,10)} (tau={mp.nstr(tau,6)})")
        for k in (2, 3):
            T = mech_terms(F, k)
            line = (f"  k={k} delta={mp.nstr(T['delta'],6)}  I={mp.nstr(T['I'],6)} Ch0={mp.nstr(T['Ch0'],4)} "
                    f"PARg0={mp.nstr(T['PARg0'],6)} PSI={mp.nstr(T['PSI'],6)} PHIINT={mp.nstr(T['PHIINT'],6)} "
                    f"h(i)={mp.nstr(T['hi'],6)}")
            print(line)
            msg = f"        B_mech = {mp.nstr(T['Bmech'],8)}"
            if d in KNOWN:
                ell = KNOWN[d](T['g0'], k, T['delta'])
                msg += f"   + known NCE+CE {mp.nstr(ell,6)}  => B = {mp.nstr(T['Bmech']+ell,8)}"
            else:
                msg += "   (+ NCE >= 0 still to add)"
            print(msg)
