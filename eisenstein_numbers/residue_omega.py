"""Stage 8 foundation — residue rings and P¹ for Γ₀(𝔫) ⊂ PSL(2, ℤ[ω]).

Parallel to independent_exclusion/congruence_prototype.py FieldResidueRing
and P1_ACTION.md, for O = ℤ[ω].

Norm: N(a + b ω) = a² − a b + b².
Units: ±1, ±ω, ±ω² (order 6).

Prime behaviour in ℤ[ω]:
  • 3 ramifies:  (3) = (λ)² up to units,  λ = 1−ω,  N(λ)=3
  • p ≡ 2 (mod 3) inert: remains prime, N(p)=p²
  • p ≡ 1 (mod 3) splits: p = π π̄,  N(π)=p

First ladder candidates (prime ideal 𝔭, index = N(𝔭)+1 for PSL Γ₀):

  | 𝔭 (gen)     | N(𝔭) | index | R = O/𝔭        |
  |-------------|------|-------|----------------|
  | λ=1−ω       | 3    | 4     | 𝔽₃             |
  | 2 (inert)   | 4    | 5     | 𝔽₄ ≅ 𝔽₂[ω]/(2)|
  | π|7         | 7    | 8     | 𝔽₇             |
  | π|13        | 13   | 14    | 𝔽₁₃            |

This module: ring ops, field residues for split/ramified prime norms,
P¹(R) transversals and right action.  Not a certificate.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

# ---------------------------------------------------------------------------
# Ring O = ℤ[ω]
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Eisenstein:
    """a + b·ω with a,b ∈ ℤ,  ω² + ω + 1 = 0."""
    a: int
    b: int

    def __repr__(self) -> str:
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            return f"{self.b}*w" if abs(self.b) != 1 else ("w" if self.b == 1 else "-w")
        sign = "+" if self.b > 0 else "-"
        bb = abs(self.b)
        return f"{self.a}{sign}{'' if bb == 1 else bb}w"

    def __add__(self, o: Eisenstein) -> Eisenstein:
        return Eisenstein(self.a + o.a, self.b + o.b)

    def __neg__(self) -> Eisenstein:
        return Eisenstein(-self.a, -self.b)

    def __sub__(self, o: Eisenstein) -> Eisenstein:
        return self + (-o)

    def __mul__(self, o: Eisenstein) -> Eisenstein:
        a, b, c, d = self.a, self.b, o.a, o.b
        # w^2 = -1 - w
        return Eisenstein(a * c - b * d, a * d + b * c - b * d)

    def conj(self) -> Eisenstein:
        # conj(w)=w^2=-1-w ⇒ conj(a+bw)=(a-b)-b w
        return Eisenstein(self.a - self.b, -self.b)

    def norm(self) -> int:
        return self.a * self.a - self.a * self.b + self.b * self.b

    def mod_int(self, n: int) -> Eisenstein:
        return Eisenstein(self.a % n, self.b % n)


W = Eisenstein(0, 1)
ONE = Eisenstein(1, 0)
ZERO = Eisenstein(0, 0)
LAMBDA = Eisenstein(1, -1)  # 1−ω, N=3


def units() -> list[Eisenstein]:
    w2 = W * W
    return [ONE, -ONE, W, -W, w2, -w2]


def associates(z: Eisenstein) -> list[Eisenstein]:
    return [u * z for u in units()]


# ---------------------------------------------------------------------------
# Residue field R = O/π  for N(π)=q prime (split or ramified)
# Represent as 𝔽_q via a fixed O → 𝔽_q with ω ↦ w_img
# ---------------------------------------------------------------------------


def find_omega_image(q: int) -> int | None:
    """Find w in 0..q-1 with w²+w+1 ≡ 0 (mod q). None if no root."""
    for w in range(q):
        if (w * w + w + 1) % q == 0:
            return w
    return None


@dataclass
class OmegaFieldRing:
    """R ≅ 𝔽_q via a+bω ↦ a + b·w_img (mod q), q = N(𝔭) prime.

    Works when O/𝔭 is a prime field (N(𝔭)=q prime), i.e. split or
    ramified primes — not inert p≡2 mod 3 (those give 𝔽_{p²}).
    """
    q: int
    w_img: int
    label: str = ""

    def __post_init__(self):
        self.norm = self.q
        self.index = self.q + 1  # |P¹(𝔽_q)| for PSL cosets of Γ₀
        if not self.label:
            self.label = f"N={self.q}"

    def embed(self, z: Eisenstein) -> int:
        return (z.a + z.b * self.w_img) % self.q

    def add(self, a: int, b: int) -> int:
        return (a + b) % self.q

    def neg(self, a: int) -> int:
        return (-a) % self.q

    def mul(self, a: int, b: int) -> int:
        return (a * b) % self.q

    def inv(self, a: int) -> int:
        a %= self.q
        if a == 0:
            raise ZeroDivisionError
        # Fermat
        return pow(a, self.q - 2, self.q)

    def is_unit(self, a: int) -> bool:
        return a % self.q != 0

    def elements(self) -> range:
        return range(self.q)

    def p1_normalize(self, c: int, d: int) -> tuple:
        """Projective class (c:d) as (0,1) or (1, x)."""
        c, d = c % self.q, d % self.q
        if c == 0:
            if d == 0:
                raise ValueError("non-unimodular (0:0)")
            return (0, 1)
        return (1, self.mul(d, self.inv(c)))

    def p1_points(self) -> list[tuple]:
        """Transversal of ℙ¹(𝔽_q): (0:1) and (1:x) for x∈𝔽_q."""
        return [(0, 1)] + [(1, x) for x in range(self.q)]

    def p1_act(self, pt: tuple, M: tuple) -> tuple:
        """Right action (c:d)·[[a,b],[c',d']] on ℙ¹."""
        c, d = pt
        a, b = M[0]
        cp, dp = M[1]
        c2 = self.add(self.mul(c, a), self.mul(d, cp))
        d2 = self.add(self.mul(c, b), self.mul(d, dp))
        return self.p1_normalize(c2, d2)

    def mat_embed(self, A: tuple[tuple[Eisenstein, Eisenstein],
                                  tuple[Eisenstein, Eisenstein]]):
        return (
            (self.embed(A[0][0]), self.embed(A[0][1])),
            (self.embed(A[1][0]), self.embed(A[1][1])),
        )


def make_ring_norm_q(q: int, label: str = "") -> OmegaFieldRing:
    """Build OmegaFieldRing for prime norm q (q prime, ω root exists)."""
    if q < 2:
        raise ValueError(q)
    w = find_omega_image(q)
    if w is None:
        raise ValueError(
            f"no ω-image mod {q}: inert or no split (need q=3 or q≡1 mod 3)"
        )
    return OmegaFieldRing(q=q, w_img=w, label=label or f"N={q}")


# Standard first-rung rings
def ring_lambda3() -> OmegaFieldRing:
    """𝔭 = (1−ω), N=3, R=𝔽₃, ω↦1 (since 1−ω≡0 ⇒ ω≡1; check 1+1+1=3≡0)."""
    return make_ring_norm_q(3, label="(1-w)")


def ring_split7() -> OmegaFieldRing:
    """Some π|7, N=7.  ω image: w²+w+1≡0 mod 7 → w=2 (4+2+1=7) or w=4."""
    return make_ring_norm_q(7, label="pi|7")


def ring_split13() -> OmegaFieldRing:
    return make_ring_norm_q(13, label="pi|13")


# Generators of PSL(2,O) embedded in residue (for gluing permutations)
# Same abstract generators as geometry_fund: T1, Tw, U, S
def gen_matrices() -> dict[str, tuple]:
    """Eisenstein matrix lifts (as Eisenstein entries)."""
    return {
        "T1": ((ONE, ONE), (ZERO, ONE)),
        "Tw": ((ONE, W), (ZERO, ONE)),
        "U": ((W * W, ZERO), (ZERO, W)),  # diag(ω², ω) → z↦ωz
        "S": ((ZERO, -ONE), (ONE, ZERO)),
    }


def gluing_perms(ring: OmegaFieldRing) -> dict[str, list[int]]:
    """For each generator δ, permutation of P¹ indices: i ↦ j where
    copy i's δ-source meets copy j (j = i · δ⁻¹ in P¹ action convention
    matching CONGRUENCE.md: j = p · δ⁻¹).

    Here we store the right action by δ itself on the transversal
    (float prototype of the permutation; cert uses the same combinatorics).
    """
    pts = ring.p1_points()
    gens = gen_matrices()
    out = {}
    for name, M_e in gens.items():
        M = ring.mat_embed(M_e)
        # inverse in PSL for gluing map p ↦ p·δ⁻¹ — compute δ⁻¹ via adjugate
        # For SL: inv = [[d,-b],[-c,a]]
        a, b = M[0]
        c, d = M[1]
        Minv = ((d, ring.neg(b)), (ring.neg(c), a))
        perm = []
        for p in pts:
            p2 = ring.p1_act(p, Minv)
            perm.append(pts.index(p2))
        out[name] = perm
    return out


# ---------------------------------------------------------------------------
# Ladder inventory
# ---------------------------------------------------------------------------

LADDER_RUNGS = [
    dict(name="(1-w)", N=3, index=4, note="ramified; first rung"),
    dict(name="pi|7", N=7, index=8, note="split p≡1 mod 3"),
    dict(name="pi|13", N=13, index=14, note="split; parallel to Gaussian N=13"),
    dict(name="(2) inert", N=4, index=5, note="𝔽₄ — needs quad residue ring"),
]


def inventory() -> None:
    print("Eisenstein congruence ladder — residue inventory")
    print("=" * 56)
    print(f"  λ=1−ω = {LAMBDA},  N(λ)={LAMBDA.norm()}")
    print(f"  units ({len(units())}): {units()}")
    print()
    for row in LADDER_RUNGS:
        print(f"  {row['name']:12}  N={row['N']:3}  index={row['index']:3}  "
              f"#  {row['note']}")
    print()
    for factory, tag in (
        (ring_lambda3, "N=3"),
        (ring_split7, "N=7"),
        (ring_split13, "N=13"),
    ):
        R = factory()
        pts = R.p1_points()
        assert len(pts) == R.index
        perms = gluing_perms(R)
        print(f"  [{tag}] w_img={R.w_img}  |P1|={len(pts)}  "
              f"T1 perm={perms['T1']}  S perm={perms['S']}")
        # S should be an involution on P1
        s = perms["S"]
        assert all(s[s[i]] == i for i in range(len(s))), "S^2 != id"
        print(f"         S^2=id OK;  U order check on P1...")
        u = perms["U"]
        # U^3 = id on P1 (up to proj)
        def comp(p, q):
            return [p[q[i]] for i in range(len(p))]
        u3 = comp(comp(u, u), u)
        assert u3 == list(range(len(u))), (u, u3)
        print(f"         U^3=id on P1 OK")
    print()
    print("  inert N=4 (p=2) not yet: need F_4 arithmetic.")
    print("  Next: multi-cusp data for Γ₀(λ), float CR on 4 copies.")


def main():
    inventory()


if __name__ == "__main__":
    main()
