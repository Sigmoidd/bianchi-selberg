"""Exact field identity, characters and splitting; no analytic dependencies."""
from dataclasses import dataclass

DISCRIMINANTS = {1: -4, 2: -8, 3: -3, 7: -7, 11: -11, 19: -19,
                 43: -43, 67: -67, 163: -163}


@dataclass(frozen=True)
class QuadraticField:
    d: int

    def __post_init__(self):
        if not isinstance(self.d, int) or isinstance(self.d, bool) or self.d not in DISCRIMINANTS:
            raise ValueError(f"unsupported field d={self.d}")

    @property
    def D(self):
        return DISCRIMINANTS[self.d]

    @property
    def units(self):
        return {1: 4, 3: 6}.get(self.d, 2)

    def character(self, n):
        if self.D == -4:
            return {1: 1, 3: -1}.get(n % 4, 0)
        if self.D == -8:
            return {1: 1, 3: 1, 5: -1, 7: -1}.get(n % 8, 0)
        p = -self.D
        r = n % p
        return 0 if r == 0 else (1 if pow(r, (p-1)//2, p) == 1 else -1)

    @property
    def character_data(self):
        q = -self.D
        return q, {a: self.character(a) for a in range(1, q) if self.character(a)}

    def covolume(self):
        from flint import arb
        return arb(-self.D).sqrt()/2

    def prime_norms(self, p):
        """(prime ideal norm, multiplicity) above a rational prime p."""
        if self.D % p == 0:
            return p, 1
        return (p, 2) if self.character(p) == 1 else (p*p, 1)


def get_field(kind):
    if isinstance(kind, QuadraticField):
        return kind
    return QuadraticField({"i": 1, "omega": 3}.get(kind, kind))


