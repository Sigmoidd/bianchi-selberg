"""Exact comparison of real quadratic radicals, shared by inventory proofs."""
from dataclasses import dataclass
from math import isqrt


@dataclass(frozen=True)
class Radical:
    a: int
    b: int
    d: int
    denominator: int = 1

    def __post_init__(self):
        if (any(not isinstance(x, int) or isinstance(x, bool)
                for x in (self.a, self.b, self.d, self.denominator)) or self.denominator <= 0):
            raise ValueError("radical requires integer coefficients and positive denominator")

    def sign(self):
        if self.d <= 0 or isqrt(self.d)**2 == self.d:
            raise ValueError("radicand must be a positive nonsquare integer")
        a, b, d = self.a, self.b, self.d
        if b == 0:
            return (a > 0)-(a < 0)
        if a >= 0 and b > 0:
            return 1
        if a <= 0 and b < 0:
            return -1
        comparison = a*a-d*b*b
        if comparison == 0:
            raise ArithmeticError("a nontrivial rational square equals a nonsquare")
        return ((comparison > 0)-(comparison < 0))*(1 if a > 0 else -1)

    def compare(self, other):
        if isinstance(other, int):
            other = Radical(other, 0, self.d)
        if self.d != other.d:
            raise ValueError("radicands differ")
        return Radical(self.a*other.denominator-other.a*self.denominator,
                       self.b*other.denominator-other.b*self.denominator, self.d).sign()

    def arb(self):
        from flint import arb
        return (arb(self.a)+self.b*arb(self.d).sqrt())/self.denominator

    def square(self):
        from math import gcd
        a, b, denominator = self.a**2+self.d*self.b**2, 2*self.a*self.b, self.denominator**2
        common = gcd(gcd(abs(a), abs(b)), denominator)
        return Radical(a//common, b//common, self.d, denominator//common)

    def unit_reciprocal(self):
        if self.a**2-self.d*self.b**2 != self.denominator**2:
            raise ValueError("radical does not have conjugate product one")
        return Radical(self.a, -self.b, self.d, self.denominator)
