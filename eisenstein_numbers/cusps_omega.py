"""Cusp data for Γ₀(𝔭) ⊂ PSL(2, ℤ[ω]) — Stage 8.

**Proved** (see T0_AREA.md): for prime 𝔭 with N(𝔭)=q and O/𝔭 a field,

  |T_∞| = √3/6,
  |T_0| = q · √3/6 = N(𝔭) · |T_∞|.

Cosets ℙ¹(O/𝔭): cusp ∞ has 1 point (0:1); cusp 0 has q points (1:x).

  β_∞(s) = (1−s)/(|T_∞| Y²),   β_0(s) = (1−s)/(|T_0| Y²),
  κ_c(s) = 1/((1+s) Y²).

Mode bounds (design check): |μ|_∞ = 2/√3, |μ|_0 = (2/√3)/√q.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from residue_omega import ring_lambda3, ring_split7, ring_split13

try:
    from reference_cell import AREA_T as AREA_T_LEVEL1
except Exception:
    AREA_T_LEVEL1 = math.sqrt(3) / 6.0

Y_DEFAULT = 1.25
SQRT3 = math.sqrt(3.0)


@dataclass(frozen=True)
class CuspTable:
    """Two-cusp data for prime Γ₀(𝔭), N(𝔭)=q."""
    q: int
    label: str
    area_inf: float
    area_0: float
    n_inf: int  # cosets at ∞
    n_0: int    # cosets at 0

    @property
    def index(self) -> int:
        return self.n_inf + self.n_0

    def beta_inf(self, s: float, Y: float = Y_DEFAULT) -> float:
        return (1.0 - s) / (self.area_inf * Y * Y)

    def beta_0(self, s: float, Y: float = Y_DEFAULT) -> float:
        return (1.0 - s) / (self.area_0 * Y * Y)

    def mode_bound_inf(self, Y: float = Y_DEFAULT) -> float:
        mu = 2.0 / SQRT3
        return 4 * math.pi ** 2 * Y * Y / (mu * mu)

    def mode_bound_0(self, Y: float = Y_DEFAULT) -> float:
        mu = (2.0 / SQRT3) / math.sqrt(self.q)
        return 4 * math.pi ** 2 * Y * Y / (mu * mu)

    def ok_modes(self, Y: float = Y_DEFAULT) -> bool:
        return self.mode_bound_inf(Y) > 1.0 and self.mode_bound_0(Y) > 1.0


def cusp_table_prime(q: int, label: str = "") -> CuspTable:
    """Build two-cusp table for prime norm q."""
    a_inf = AREA_T_LEVEL1  # √3/6
    a_0 = q * AREA_T_LEVEL1
    return CuspTable(
        q=q,
        label=label or f"N={q}",
        area_inf=a_inf,
        area_0=a_0,
        n_inf=1,
        n_0=q,
    )


def tables() -> dict[str, CuspTable]:
    return {
        "(1-w)": cusp_table_prime(3, "(1-w)"),
        "pi|7": cusp_table_prime(7, "pi|7"),
        "pi|13": cusp_table_prime(13, "pi|13"),
    }


def verify_cusp_classes():
    """Check P¹ partitions match (1, q) for field levels."""
    from collections import Counter
    rows = []
    for factory, name, q in (
        (ring_lambda3, "(1-w)", 3),
        (ring_split7, "pi|7", 7),
        (ring_split13, "pi|13", 13),
    ):
        R = factory()
        pts = R.p1_points()
        # ∞ iff first coord 0
        cls = [0 if p[0] == 0 else 1 for p in pts]
        ctr = Counter(cls)
        assert ctr[0] == 1 and ctr[1] == q, (name, ctr)
        rows.append((name, q, dict(ctr)))
    return rows


def main():
    print("Cusp tables for Γ₀(𝔭) ⊂ PSL(2,ℤ[ω])")
    print("=" * 56)
    print(f"  level-1 |T| = √3/6 = {AREA_T_LEVEL1:.12f}")
    print()
    verify_cusp_classes()
    print("  P¹ cusp partition (∞ : 0) = (1 : N) for field primes: OK")
    print()
    Y = Y_DEFAULT
    for name, ct in tables().items():
        print(f"  {name:8}  N={ct.q:3}  index={ct.index:3}")
        print(f"    |T_∞|={ct.area_inf:.6f}  |T_0|={ct.area_0:.6f}")
        print(f"    mode ∞: {ct.mode_bound_inf(Y):.1f}  "
              f"mode 0: {ct.mode_bound_0(Y):.1f}  "
              f"OK={ct.ok_modes(Y)}")
        print(f"    β_∞(0)={ct.beta_inf(0, Y):.4f}  "
              f"β_0(0)={ct.beta_0(0, Y):.4f}")
    print()
    print("  |T_0|=N·|T_∞| proved for field primes: T0_AREA.md")


if __name__ == "__main__":
    main()
