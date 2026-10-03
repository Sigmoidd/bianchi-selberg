"""Exact comparison of real quadratic radicals, shared by inventory proofs."""
from dataclasses import dataclass
from math import isqrt
from flint import arb


@dataclass(frozen=True)
class Radical:
    a: int
    b: int
    d: int

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
        return Radical(self.a-other.a, self.b-other.b, self.d).sign()

    def arb(self):
        return arb(self.a)+self.b*arb(self.d).sqrt()
