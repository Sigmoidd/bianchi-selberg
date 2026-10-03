"""Field-dependent Euler-product term in the shifted scattering integral."""
from math import floor
from flint import arb
from fields.quadratic import get_field, rational_primes
from core.bspline import gbspline


def prime_term(kind, k, delta):
    F = get_field(kind)
    support = 2*k*arb(delta)
    # Round outwards and include boundary candidates; gbspline safely encloses
    # the truncated power even when log(N) straddles the support endpoint.
    limit = floor(support.exp().upper())
    total = arb(0)
    for p in rational_primes(limit):
        norm, multiplicity = F.prime_norms(p)
        N = norm
        while N <= limit:
            total += multiplicity*arb(norm).log()/N*gbspline(arb(N).log(), k, delta)
            N *= norm
    return total
