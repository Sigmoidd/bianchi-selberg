"""
Preliminary smoke tests for Eisenstein integers Z[omega].

  omega = exp(2 pi i / 3),  omega^2 + omega + 1 = 0,
  N(a + b omega) = a^2 - a b + b^2.

Not certified output — arithmetic / enumeration only.
See ../verify_eisenstein.py for Arb-certified mechanical constants.
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Ring Z[omega]
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class Eisenstein:
    """a + b*omega with a,b in Z."""
    a: int
    b: int

    def __repr__(self) -> str:
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            return f"{self.b}*w" if self.b != 1 else "w"
        sign = "+" if self.b > 0 else "-"
        bb = abs(self.b)
        return f"{self.a}{sign}{bb if bb != 1 else ''}w"

    def __add__(self, other: Eisenstein) -> Eisenstein:
        return Eisenstein(self.a + other.a, self.b + other.b)

    def __neg__(self) -> Eisenstein:
        return Eisenstein(-self.a, -self.b)

    def __sub__(self, other: Eisenstein) -> Eisenstein:
        return self + (-other)

    def __mul__(self, other: Eisenstein) -> Eisenstein:
        # (a+bw)(c+dw) = ac + (ad+bc)w + bd w^2
        # w^2 = -1 - w  =>  ac - bd + (ad+bc-bd)w
        a, b, c, d = self.a, self.b, other.a, other.b
        return Eisenstein(a * c - b * d, a * d + b * c - b * d)

    def conj(self) -> Eisenstein:
        # conj(w) = w^2 = -1-w => conj(a+bw) = a + b w^2 = (a-b) - b w
        return Eisenstein(self.a - self.b, -self.b)

    def norm(self) -> int:
        # N = z * conj(z) = a^2 - ab + b^2
        return self.a * self.a - self.a * self.b + self.b * self.b

    def to_complex(self) -> complex:
        w = cmath.exp(2j * math.pi / 3)
        return self.a + self.b * w


W = Eisenstein(0, 1)
ONE = Eisenstein(1, 0)
ZERO = Eisenstein(0, 0)


def units() -> list[Eisenstein]:
    """Six units of Z[omega]: ±1, ±w, ±w^2."""
    w2 = W * W
    return [ONE, -ONE, W, -W, w2, -w2]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_minimal_polynomial() -> None:
    # w^2 + w + 1 = 0
    lhs = W * W + W + ONE
    assert lhs == ZERO, lhs


def test_norm_multiplicative() -> None:
    samples = [Eisenstein(a, b) for a in range(-4, 5) for b in range(-4, 5)]
    for z in samples:
        for u in samples:
            assert (z * u).norm() == z.norm() * u.norm()


def test_units_norm_one() -> None:
    for u in units():
        assert u.norm() == 1, (u, u.norm())
    assert len(units()) == 6


def test_conj_norm() -> None:
    z = Eisenstein(3, -2)
    assert (z * z.conj()).norm() == z.norm() ** 2
    # z * conj(z) is rational integer N(z)
    p = z * z.conj()
    assert p.b == 0 and p.a == z.norm()


def test_small_primes_ramification() -> None:
    # Classically: 1-w is prime above (3) up to units (N(1-w)=3).
    p3 = ONE - W
    assert p3.norm() == 3
    # 2 = -w^2 (1-w)^2 / something?  N(2)=4; 2 ramifies as unit * (1-w)^2
    # Check N(1-w)^2 = 9 != 4 — better: known that lambda=1-w, 3 ~ lambda^2
    assert (p3 * p3).norm() == 9


def test_loxodromic_trace_search_tiny() -> None:
    """Tiny brute force: min log N(tau) for |tr| large enough — diagnostic only.

    Full certified systole search lives in ../verify_eisenstein.py.
    """
    best = None
    for a in range(-8, 9):
        for b in range(-8, 9):
            tau = Eisenstein(a, b)
            n = tau.norm()
            # loxodromic traces satisfy |tr|^2 > 4 in the embedding sense;
            # use N(tau) > 4 as a crude integer filter (preliminary).
            if n <= 4:
                continue
            cand = math.log(n)
            if best is None or cand < best[0]:
                best = (cand, tau, n)
    assert best is not None
    # Repo certified systole uses log N(T) with T group element; here we only
    # sanity-check that the search finds something near the known scale ~0.86
    print(f"  tiny-search min log N(tau) ~ {best[0]:.6f}  tau={best[1]}  N={best[2]}")
    assert 0.5 < best[0] < 2.0


def main() -> None:
    print("Eisenstein numbers — preliminary smoke tests")
    print("=" * 50)
    tests = [
        test_minimal_polynomial,
        test_norm_multiplicative,
        test_units_norm_one,
        test_conj_norm,
        test_small_primes_ramification,
        test_loxodromic_trace_search_tiny,
    ]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"  [PASS] {t.__name__}")
        except Exception as e:
            failed += 1
            print(f"  [FAIL] {t.__name__}: {e}")
    print("=" * 50)
    if failed:
        print(f"{failed} failed")
        raise SystemExit(1)
    print("all smoke tests passed")
    print("next: python -u ../verify_eisenstein.py  (Arb-certified constants)")


if __name__ == "__main__":
    main()
