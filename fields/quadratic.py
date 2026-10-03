"""Exact characters for the supported class-number-one fields.

No elliptic inventory is inferred from the field discriminant.
"""
from math import isqrt

from flint import arb, acb, acb_series

from fields.arithmetic import DISCRIMINANTS, QuadraticField, get_field


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
