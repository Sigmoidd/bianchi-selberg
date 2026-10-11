#!/usr/bin/env python3
r"""Exact Q(i) arithmetic and exact Poincare extension for the global-symmetry work.

Everything here is Fraction arithmetic.  No floating point is used.

A Gaussian rational is a pair ``(re, im)`` of ``Fraction``.  A 2x2 matrix is a
4-tuple ``(a, b, c, d)`` of Gaussian rationals.  The upper-half-space action
uses the Poincare extension

    (z, y) -> ( ((a z + b) conj(c z + d) + a conj(c) y^2) / D ,  |det| y / D ),
    D = |c z + d|^2 + |c|^2 y^2 ,

which is the action of M / sqrt(det M) in SL2(C).  The height formula is exact
only when |det M| is rational; every matrix used with a point action in this
project has |det| = 1.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Iterable

G = tuple[Fraction, Fraction]
Mat = tuple[G, G, G, G]
Point = tuple[G, Fraction]


def g(re: int | Fraction = 0, im: int | Fraction = 0) -> G:
    return (Fraction(re), Fraction(im))


ZERO = g(0, 0)
ONE = g(1, 0)
I = g(0, 1)


def gadd(a: G, b: G) -> G:
    return (a[0] + b[0], a[1] + b[1])


def gsub(a: G, b: G) -> G:
    return (a[0] - b[0], a[1] - b[1])


def gneg(a: G) -> G:
    return (-a[0], -a[1])


def gmul(a: G, b: G) -> G:
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gconj(a: G) -> G:
    return (a[0], -a[1])


def gnorm(a: G) -> Fraction:
    return a[0] * a[0] + a[1] * a[1]


def ginv(a: G) -> G:
    n = gnorm(a)
    if n == 0:
        raise ZeroDivisionError("inverse of zero Gaussian rational")
    return (a[0] / n, -a[1] / n)


def gdiv(a: G, b: G) -> G:
    return gmul(a, ginv(b))


def is_gaussian_integer(a: G) -> bool:
    return a[0].denominator == 1 and a[1].denominator == 1


def mat(a: G, b: G, c: G, d: G) -> Mat:
    return (a, b, c, d)


def mmul(m: Mat, n: Mat) -> Mat:
    return (
        gadd(gmul(m[0], n[0]), gmul(m[1], n[2])),
        gadd(gmul(m[0], n[1]), gmul(m[1], n[3])),
        gadd(gmul(m[2], n[0]), gmul(m[3], n[2])),
        gadd(gmul(m[2], n[1]), gmul(m[3], n[3])),
    )


def mdet(m: Mat) -> G:
    return gsub(gmul(m[0], m[3]), gmul(m[1], m[2]))


def minv(m: Mat) -> Mat:
    """Inverse in GL2 (adjugate over determinant)."""
    det_inv = ginv(mdet(m))
    return (
        gmul(m[3], det_inv),
        gmul(gneg(m[1]), det_inv),
        gmul(gneg(m[2]), det_inv),
        gmul(m[0], det_inv),
    )


def mconj(m: Mat) -> Mat:
    """Entrywise Gaussian conjugate."""
    return tuple(gconj(x) for x in m)  # type: ignore[return-value]


def mneg(m: Mat) -> Mat:
    return tuple(gneg(x) for x in m)  # type: ignore[return-value]


def proj_equal(m: Mat, n: Mat) -> bool:
    """Equality in PSL2 (up to the scalar -1) for SL2 matrices."""
    return m == n or m == mneg(n)


def proj_equal_gl(m: Mat, n: Mat) -> bool:
    """Equality in PGL2 for GL2 matrices: m = lambda n for a scalar lambda."""
    # Cross-multiplication test on all 2x2 minors of the stacked 2x4 matrix.
    cross = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for i, j in cross:
        if gsub(gmul(m[i], n[j]), gmul(m[j], n[i])) != ZERO:
            return False
    return any(x != ZERO for x in m) and any(x != ZERO for x in n)


def h3_action(m: Mat, p: Point) -> Point:
    """Exact Poincare-extension action of the GL2 matrix ``m`` on H^3."""
    a, b, c, d = m
    z, y = p
    q = gadd(gmul(c, z), d)
    den = gnorm(q) + gnorm(c) * y * y
    if den <= 0:
        raise ArithmeticError("non-positive Moebius denominator")
    num = gadd(gmul(gadd(gmul(a, z), b), gconj(q)), gmul(gmul(a, gconj(c)), g(y * y, 0)))
    det_abs = gnorm(mdet(m))
    # |det| must be rational for an exact height; all callers use |det| = 1.
    root = _exact_sqrt(det_abs)
    return ((num[0] / den, num[1] / den), root * y / den)


def _exact_sqrt(x: Fraction) -> Fraction:
    import math
    if x < 0:
        raise ValueError("negative")
    n, d = x.numerator, x.denominator
    rn, rd = math.isqrt(n), math.isqrt(d)
    if rn * rn != n or rd * rd != d:
        raise ArithmeticError(f"|det|={x} is not a rational square; exact height unavailable")
    return Fraction(rn, rd)


# ---------------------------------------------------------------- generators

T1: Mat = (ONE, ONE, ZERO, ONE)
TI: Mat = (ONE, I, ZERO, ONE)
R: Mat = (I, ZERO, ZERO, g(0, -1))
S: Mat = (ZERO, g(-1), ONE, ZERO)
TIR: Mat = mmul(TI, R)
GENERATORS: dict[str, Mat] = {"T1": T1, "R": R, "TiR": TIR, "S": S}

QMAT: Mat = (I, ZERO, ZERO, ONE)          # diag(i, 1): z -> i z, y -> y
QSQ_REP: Mat = (I, ZERO, ZERO, g(0, -1))  # diag(i,-i): z -> -z (equals R)

PI: G = g(2, 1)
PI_BAR: G = g(2, -1)


def rotate_point(p: Point, turns: int) -> Point:
    """Apply z -> i^turns z (an isometry fixing y)."""
    z, y = p
    for _ in range(turns % 4):
        z = gmul(I, z)
    return (z, y)


def conj_point(p: Point) -> Point:
    """C(z, y) = (conj z, y)."""
    return (gconj(p[0]), p[1])


def random_points(count: int, seed: int) -> list[Point]:
    import random
    rng = random.Random(seed)
    out = []
    for _ in range(count):
        z = g(Fraction(rng.randint(-9, 9), rng.randint(1, 7)),
              Fraction(rng.randint(-9, 9), rng.randint(1, 7)))
        y = Fraction(rng.randint(1, 12), rng.randint(1, 7))
        out.append((z, y))
    return out


def random_word(rng, length: int, letters: Iterable[Mat] | None = None) -> Mat:
    """Random product of the four generators and their inverses (an element of PSL2(Z[i]))."""
    pool = list(letters) if letters is not None else [T1, TI, R, S]
    m: Mat = (ONE, ZERO, ZERO, ONE)
    for _ in range(length):
        h = rng.choice(pool)
        if rng.random() < 0.5:
            h = minv(h)
        m = mmul(m, h)
    return m


def divisible_by(a: G, prime: G) -> bool:
    """True iff the Gaussian integer ``a`` lies in the ideal (prime)."""
    if not is_gaussian_integer(a):
        return False
    return is_gaussian_integer(gdiv(a, prime))
