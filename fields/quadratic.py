"""Exact characters for the supported class-number-one fields.

No elliptic inventory is inferred from the field discriminant.
"""
from dataclasses import dataclass
from math import isqrt

from flint import arb, acb, acb_series

DISCRIMINANTS = {1: -4, 2: -8, 3: -3, 7: -7, 11: -11, 19: -19,
                 43: -43, 67: -67, 163: -163}


@dataclass(frozen=True)
class QuadraticField:
    d: int

    def __post_init__(self):
        if self.d not in DISCRIMINANTS:
            raise ValueError(f"unsupported field d={self.d}")

    @property
    def D(self):
        return DISCRIMINANTS[self.d]

    @property
    def units(self):
        return {1: 4, 3: 6}.get(self.d, 2)

    def character(self, n):
        if self.D == -4:
            return {1: 1, 3: -1}.get(n % 4, 0)
        if self.D == -8:
            return {1: 1, 3: 1, 5: -1, 7: -1}.get(n % 8, 0)
        p = -self.D
        r = n % p
        return 0 if r == 0 else (1 if pow(r, (p-1)//2, p) == 1 else -1)

    @property
    def character_data(self):
        q = -self.D
        return q, {a: self.character(a) for a in range(1, q) if self.character(a)}

    def covolume(self):
        return arb(-self.D).sqrt()/2

    def prime_norms(self, p):
        """(prime ideal norm, multiplicity) above a rational prime p."""
        if self.D % p == 0:
            return p, 1
        return (p, 2) if self.character(p) == 1 else (p*p, 1)


def get_field(kind):
    if isinstance(kind, QuadraticField):
        return kind
    return QuadraticField({"i": 1, "omega": 3}.get(kind, kind))


def Lval1(kind):
    q, chi = get_field(kind).character_data
    return (-(arb(1)/q)*sum(c*(acb(a)/q).digamma() for a, c in chi.items())).real


def Lprime1(kind, rho=0.4):
    q, chi = get_field(kind).character_data
    rho = arb(rho)
    if not (rho > 0 and rho < 1):
        raise ValueError("Cauchy radius must lie strictly between 0 and 1")
    def Lf(z):
        return acb(q)**(-z)*sum(c*z.zeta(acb(a)/q) for a, c in chi.items())
    def integrand(th, _):
        e = (acb(0, 1)*th).exp()
        z = 1 + rho*e
        return Lf(z)/(z-1)**2*rho*e*acb(0, 1)
    value = acb.integral(integrand, 0, 2*arb.pi())
    return (value/(acb(0, 1)*2*arb.pi())).real


def zetaK2(kind):
    q, chi = get_field(kind).character_data
    L2 = acb(q)**(-acb(2))*sum(c*acb(2).zeta(acb(a)/q) for a, c in chi.items())
    return (acb(2).zeta()*L2).real


def volume(kind):
    F = get_field(kind)
    return arb(-F.D)**arb("1.5")*zetaK2(F)/(4*arb.pi()**2)


def eta(kind):
    F = get_field(kind)
    c0 = F.units*(arb.const_euler()*Lval1(F)+Lprime1(F))
    return F.covolume()/arb.pi()*c0


def _Lseries(x, q, chi):
    logq = arb(q).log()
    qser = acb_series([-x[0]*logq, -logq]).exp()
    hur = acb_series([0, 0])
    for a, c in chi.items():
        hur += c*x.zeta(acb(a)/q)
    return qser*hur


def zetaK_logder(s, kind):
    q, chi = get_field(kind).character_data
    x = acb_series([s, 1])
    zk = x.zeta()*_Lseries(x, q, chi)
    return zk[1]/zk[0]


def FpF(s, kind):
    q, chi = get_field(kind).character_data
    x = acb_series([s, 1])
    zk = x.zeta()*_Lseries(x, q, chi)
    value = acb_series([s-1, 1])*zk
    return value[1]/value[0]


def rational_primes(limit):
    for p in range(2, limit+1):
        if all(p % q for q in range(2, isqrt(p)+1)):
            yield p
