"""Shared exact 2x2 matrix operations in an integral quadratic ring."""
from dataclasses import dataclass
from groups.arithmetic import QuadraticRing

Entry = tuple[int, int]
Matrix = tuple[Entry, Entry, Entry, Entry]
IDENTITY: Matrix = ((1, 0), (0, 0), (0, 0), (1, 0))


@dataclass(frozen=True)
class MatrixOps:
    ring: QuadraticRing

    @staticmethod
    def validate(matrix):
        if not isinstance(matrix, tuple) or len(matrix) != 4 or any(
            not isinstance(e, tuple) or len(e) != 2 or any(
                not isinstance(c, int) or isinstance(c, bool) for c in e) for e in matrix):
            raise ValueError("matrix requires four integer coordinate pairs in row-major order")

    def mul(self, A, B):
        a, b, c, d = A
        e, f, g, h = B
        r = self.ring
        return (r.add(r.mul(a, e), r.mul(b, g)), r.add(r.mul(a, f), r.mul(b, h)),
                r.add(r.mul(c, e), r.mul(d, g)), r.add(r.mul(c, f), r.mul(d, h)))

    def neg(self, A):
        return tuple(self.ring.neg(x) for x in A)

    def det(self, A):
        return self.ring.sub(self.ring.mul(A[0], A[3]), self.ring.mul(A[1], A[2]))

    def canon(self, A):
        return min(A, self.neg(A))

    def power(self, A, n):
        if not isinstance(n, int) or isinstance(n, bool) or n < 0:
            raise ValueError("matrix exponent must be a nonnegative integer")
        result = IDENTITY
        while n:
            if n & 1:
                result = self.mul(result, A)
            A = self.mul(A, A)
            n //= 2
        return result
