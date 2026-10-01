"""Fourier transform of sinc^(2k), with its exact second derivative."""
from fractions import Fraction
from math import comb, factorial

from flint import arb


def validate_k(k):
    # k=1 has a triangular Fourier transform, and int h(r) r^2 dr diverges.
    if isinstance(k, bool) or not isinstance(k, int) or k < 2:
        raise ValueError("the identity term requires an integer k >= 2")


def gbspline(x, k, d):
    validate_k(k)
    x, d = arb(x), arb(d)
    if not d > 0:
        raise ValueError("delta must be provably positive")
    n = 2*k
    total = arb(0)
    for j in range(n+1):
        y = x + (k-j)*2*d
        if y.lower() > 0:
            p = y**(n-1)
        elif y.upper() > 0:
            p = arb(0).union(arb(y.upper())**(n-1))
        else:
            continue
        total += (-1)**j * comb(n, j) * p
    return total / (factorial(n-1)*(2*d)**n)


def second_derivative_coefficient(k):
    """Return the rational c_k with g''(0) = c_k / delta^3.

    Differentiate the positive truncated powers twice at zero. Since n>=4,
    the zero-base term contributes zero, leaving exactly j=0,...,k-1.
    """
    validate_k(k)
    n = 2*k
    numerator = sum((-1)**j*comb(n, j)*(k-j)**(n-3) for j in range(k))
    return Fraction(numerator, 8*factorial(n-3))


def gpp0(k, d):
    d = arb(d)
    if not d > 0:
        raise ValueError("delta must be provably positive")
    c = second_derivative_coefficient(k)
    return arb(c.numerator) / c.denominator / d**3


def g0(k, d):
    return gbspline(0, k, d)
