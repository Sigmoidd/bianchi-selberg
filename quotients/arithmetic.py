"""Exact O/I and SL2(O/I)/{+I,-I}; no trace-formula or Arb imports.

Over composite rings we quotient only by the global signs, not by all
scalar solutions of u^2=1. This is the image of the chosen PSL convention.
"""
from dataclasses import dataclass
from groups.arithmetic import QuadraticRing
from groups.identity import LevelIdeal
from groups.matrix import MatrixOps, IDENTITY


@dataclass(frozen=True)
class ResidueRing:
    ideal: LevelIdeal

    def reduce(self, value):
        x, y = value
        if any(not isinstance(v, int) or isinstance(v, bool) for v in (x, y)):
            raise ValueError("residue coordinates must be integers")
        a, b, c = self.ideal.hnf
        q, r = divmod(y, c)
        return (x-b*q) % a, r

    @property
    def elements(self):
        a, _, c = self.ideal.hnf
        return tuple((x, y) for x in range(a) for y in range(c))

    def add(self, x, y):
        return self.reduce((x[0]+y[0], x[1]+y[1]))

    def sub(self, x, y):
        return self.add(x, self.neg(y))

    def neg(self, x):
        return self.reduce((-x[0], -x[1]))

    def mul(self, x, y):
        return self.reduce(QuadraticRing(self.ideal.field_d).mul(x, y))


@dataclass(frozen=True)
class ProjectiveMatrices:
    ring: ResidueRing

    @property
    def ops(self):
        return MatrixOps(self.ring)

    @property
    def identity(self):
        return self.reduce(IDENTITY)

    def reduce(self, matrix):
        MatrixOps.validate(matrix)
        M = tuple(self.ring.reduce(x) for x in matrix)
        if self.ops.det(M) != self.ring.reduce((1, 0)):
            raise ValueError("matrix determinant is not one modulo the ideal")
        return self.ops.canon(M)

    def mul(self, A, B):
        return self.reduce(self.ops.mul(A, B))

    def inverse(self, A):
        a, b, c, d = A
        return self.reduce((d, self.ring.neg(b), self.ring.neg(c), a))
