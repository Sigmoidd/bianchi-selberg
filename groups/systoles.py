"""Finite, exhaustive verification of the systoles in docs/SYSTOLES.md."""
from fractions import Fraction
from math import gcd, isqrt
from flint import arb
from groups.data import get_group
from fields.quadratic import DISCRIMINANTS


def trace_pair(group, a, b):
    d = group.field.d
    x = Fraction(a) if group.field.D % 4 == 0 else Fraction(2*a+b, 2)
    y2 = Fraction(d*b*b) if group.field.D % 4 == 0 else Fraction(d*b*b, 4)
    A = x*x+y2
    return A, (A-4)**2+16*y2


def cosh_length(pair):
    A, rad = pair
    return (arb(A.numerator)/A.denominator
            + (arb(rad.numerator)/rad.denominator).sqrt())/4


def verify_systole(kind):
    group = get_group(kind)
    group.require_level_one()
    target = trace_pair(group, *group.systole_trace)
    v0 = cosh_length(target)
    # |tau| > 3 => cosh(length) > 7/2. All chosen targets are <= 7/2.
    if target != (Fraction(9), Fraction(25)) and not v0 < arb(7)/2:
        raise ArithmeticError("systole witness exceeds the analytic tail cutoff")
    checked = 0
    ties = []
    # |tau| <= 3 implies |b| <= 6 and |a| <= 6 in either integral basis.
    for a in range(-6, 7):
        for b in range(-6, 7):
            pair = trace_pair(group, a, b)
            if pair[0] > 9 or (b == 0 and abs(a) <= 2):
                continue
            checked += 1
            if pair == target:
                ties.append((a, b))
            elif not cosh_length(pair) > v0:
                # An overlapping interval never counts as a successful check.
                raise ArithmeticError(f"trace ({a},{b}) not proved longer than witness")
    if group.systole_trace not in ties:
        raise ArithmeticError("witness is absent from exhaustive trace list")
    return dict(d=group.field.d, trace=group.systole_trace, checked=checked,
                ties=ties, systole=group.systole())


def reduced_forms(D):
    """All primitive reduced positive definite forms of discriminant D<0.

    Reduction inequalities imply 3*a^2 <= |D|; boundary sign conventions
    remove duplicate forms. This provides a separate class-number check.
    """
    forms = []
    for a in range(1, isqrt((-D)//3)+1):
        for b in range(-a, a+1):
            if (b*b-D) % (4*a):
                continue
            c = (b*b-D)//(4*a)
            if a > c or gcd(gcd(a, b), c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            forms.append((a, b, c))
    return forms


def verify_all():
    rows = []
    for d, D in DISCRIMINANTS.items():
        if len(reduced_forms(D)) != 1:
            raise ArithmeticError(f"class-number-one assumption failed for D={D}")
        rows.append(verify_systole(d))
    return rows


if __name__ == "__main__":
    for row in verify_all():
        print(f"d={row['d']:3d}: {row['checked']:2d} loxodromic trace candidates; "
              f"witness={row['trace']}; systole={row['systole']}")
