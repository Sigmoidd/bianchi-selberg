"""Exact identity and membership for full and congruence Bianchi groups.

An ideal uses column Hermite form ((a,0),(b,c)), a,c>0 and 0<=b<a,
relative to the standard integral basis in groups.arithmetic (integer d).
"""
from dataclasses import dataclass
from fields.arithmetic import get_field
from groups.arithmetic import QuadraticRing


def _integer(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _bezout(a, b):
    old_r, r, old_x, x, old_y, y = a, b, 1, 0, 0, 1
    while r:
        q = old_r//r
        old_r, r = r, old_r-q*r
        old_x, x = x, old_x-q*x
        old_y, y = y, old_y-q*y
    sign = 1 if old_r >= 0 else -1
    return sign*old_r, sign*old_x, sign*old_y


@dataclass(frozen=True)
class LevelIdeal:
    field_d: int
    hnf: tuple[int, int, int] = (1, 0, 1)

    def __post_init__(self):
        if not _integer(self.field_d):
            raise ValueError("field d must be an integer")
        get_field(self.field_d)
        if not isinstance(self.hnf, tuple) or len(self.hnf) != 3 or not all(map(_integer, self.hnf)):
            raise ValueError("ideal HNF must be a triple of integers")
        a, b, c = self.hnf
        if a <= 0 or c <= 0 or not 0 <= b < a:
            raise ValueError("ideal HNF requires a,c>0 and 0<=b<a")
        ring = QuadraticRing(self.field_d)
        if not all(self.contains(ring.mul(v, (0, 1))) for v in self.basis):
            raise ValueError("the HNF lattice is not an integral ideal")

    @property
    def basis(self):
        a, b, c = self.hnf
        return (a, 0), (b, c)

    @property
    def norm(self):
        a, b, c = self.hnf
        return a*c

    def contains(self, value):
        x, y = value
        if not _integer(x) or not _integer(y):
            raise ValueError("ideal membership requires integral-basis integer coordinates")
        a, b, c = self.hnf
        return y % c == 0 and (x-b*(y//c)) % a == 0

    @classmethod
    def principal(cls, field_d, generator):
        ring = QuadraticRing(field_d)
        if len(generator) != 2 or not all(map(_integer, generator)):
            raise ValueError("ideal generator requires two integer coordinates")
        u, v = generator, ring.mul(generator, (0, 1))
        determinant = u[0]*v[1]-v[0]*u[1]
        if determinant == 0:
            raise ValueError("level ideal must be nonzero")
        c, x, y = _bezout(u[1], v[1])
        a = abs(determinant)//c
        b = (x*u[0]+y*v[0]) % a
        return cls(field_d, (a, b, c))

    @classmethod
    def rational(cls, field_d, n):
        if not _integer(n) or n <= 0:
            raise ValueError("rational level must be a positive integer")
        return cls(field_d, (n, 0, n))

    def payload(self):
        return dict(hnf=list(self.hnf), norm=self.norm, basis="standard-integral")


@dataclass(frozen=True)
class GroupKey:
    field_d: int
    subgroup: str = "full"
    level_hnf: tuple[int, int, int] = (1, 0, 1)

    def __post_init__(self):
        ideal = LevelIdeal(self.field_d, self.level_hnf)
        if self.subgroup not in ("full", "principal", "gamma0", "gamma1"):
            raise ValueError("subgroup must be full, principal, gamma0, or gamma1")
        if self.subgroup == "full" and ideal.norm != 1:
            raise ValueError("the full group requires the unit level ideal")
        if ideal.norm == 1:
            object.__setattr__(self, "subgroup", "full")

    @classmethod
    def congruence(cls, subgroup, ideal):
        return cls(ideal.field_d, subgroup, ideal.hnf)

    @property
    def level(self):
        return LevelIdeal(self.field_d, self.level_hnf)

    @property
    def is_full(self):
        return self.subgroup == "full"

    def contains(self, matrix):
        from groups.matrix import MatrixOps
        ring = QuadraticRing(self.field_d)
        ops = MatrixOps(ring)
        ops.validate(matrix)
        if ops.det(matrix) != (1, 0):
            return False
        if self.is_full:
            return True
        a, b, c, d = matrix
        ideal = self.level
        if not ideal.contains(c):
            return False
        if self.subgroup == "gamma0":
            return True
        for sign in (1, -1):
            diagonal = ideal.contains(ring.sub(a, (sign, 0))) and ideal.contains(ring.sub(d, (sign, 0)))
            if diagonal and (self.subgroup == "gamma1" or ideal.contains(b)):
                return True
        return False

    def payload(self):
        return dict(d=self.field_d, subgroup=self.subgroup, level=self.level.payload())
