"""Rank-four relative quadratic orders over an integral quadratic ring.

The basis is (1,g,x,g*x) with x^2=trace*x-norm. These operations prove
identities only; maximality, lattice completeness and unit reduction are
supplied by each group's arithmetic backend.
"""
from flint import fmpz_mat
from groups.arithmetic import QuadraticRing

BASIS = tuple(tuple(int(i == j) for i in range(4)) for j in range(4))
IDENTITY = BASIS[0]
ROOT = BASIS[2]


class QuadraticOrder:
    def __init__(self, ring: QuadraticRing, trace, norm):
        self.ring, self.t, self.n = ring, trace, norm

    @staticmethod
    def pair(z):
        return z[:2], z[2:]

    def mul(self, z, w):
        u, v = self.pair(z)
        U, V = self.pair(w)
        constant = self.ring.sub(self.ring.mul(u, U), self.ring.mul(self.n, self.ring.mul(v, V)))
        linear = self.ring.add(self.ring.add(self.ring.mul(u, V), self.ring.mul(v, U)),
                          self.ring.mul(self.t, self.ring.mul(v, V)))
        return constant+linear

    def norm(self, z):
        u, v = self.pair(z)
        return self.ring.add(self.ring.add(self.ring.mul(u, u), self.ring.mul(self.t, self.ring.mul(u, v))),
                        self.ring.mul(self.n, self.ring.mul(v, v)))

    def sigma(self, z):
        u, v = self.pair(z)
        return self.ring.add(u, self.ring.mul(self.t, v))+self.ring.neg(v)

    def matrix(self, z):
        u, v = self.pair(z)
        return u, self.ring.neg(self.ring.mul(self.n, v)), v, self.ring.add(u, self.ring.mul(self.t, v))

    def multiplication_matrix(self, z):
        columns = [self.mul(z, e) for e in BASIS]
        return fmpz_mat([[columns[j][i] for j in range(4)] for i in range(4)])

    def absolute_norm(self, z):
        return int(self.multiplication_matrix(z).det())

    def discriminant(self):
        return int(fmpz_mat([[sum(self.multiplication_matrix(self.mul(e, f))[i, i]
                                  for i in range(4)) for f in BASIS] for e in BASIS]).det())


def power(order, z, exponent):
    if not isinstance(exponent, int) or isinstance(exponent, bool) or exponent < 0:
        raise ValueError("order exponent must be a nonnegative integer")
    result = IDENTITY
    for _ in range(exponent):
        result = order.mul(result, z)
    return result
