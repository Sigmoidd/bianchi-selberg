"""Integral quadratic-ring arithmetic without a floating point embedding."""
from dataclasses import dataclass
from fields.arithmetic import get_field


@dataclass(frozen=True)
class QuadraticRing:
    kind: object

    def __post_init__(self):
        get_field(self.kind)

    @property
    def relation(self):
        # g^2 = t*g - n; omega uses the legacy basis exp(2*pi*i/3).
        if self.kind == "omega":
            return -1, 1
        F = get_field(self.kind)
        return (0, F.d) if F.D % 4 == 0 else (1, (F.d+1)//4)

    def mul(self, x, y):
        a, b = x
        c, d = y
        t, n = self.relation
        return a*c-n*b*d, a*d+b*c+t*b*d

    def add(self, x, y):
        return x[0]+y[0], x[1]+y[1]

    def sub(self, x, y):
        return x[0]-y[0], x[1]-y[1]

    def neg(self, x):
        return -x[0], -x[1]

    def norm(self, x):
        a, b = x
        t, n = self.relation
        return a*a+t*a*b+n*b*b
