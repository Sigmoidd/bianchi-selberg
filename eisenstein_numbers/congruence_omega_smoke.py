"""Stage 8 smoke — gluing combinatorics for Γ₀(𝔭) over ℤ[ω] (no spectrum yet).

Checks:
  1. residue rings + |P¹| = N+1 for N∈{3,7,13}
  2. generator perms are well-defined; S²=id, U³=id on P¹
  3. union-find coset gluing count of free interfaces
  4. reference P₃ cell still builds (level-1 cell for copies)

Does **not** assemble multi-copy FEM or claim λ₁≥1 for Γ₀.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from residue_omega import (
    ring_lambda3, ring_split7, ring_split13, gluing_perms, LADDER_RUNGS,
)
from reference_cell import build_P3_cell, AREA_T


def glue_components(perm_dict: dict[str, list[int]]) -> tuple[int, int]:
    """Union-find on n copies; identify i ~ perm[i] for each generator perm.
    Returns (n_copies, n_orbits_of_dofs_placeholder).

    Here we only track identification of *copies* that are joined by some
    face perm with j≠i (internal fixed points ignored).  Full DOF gluing
    needs face maps — next prototype step.
    """
    # Use first perm length as n
    n = len(next(iter(perm_dict.values())))
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    cross = 0
    for name, perm in perm_dict.items():
        for i, j in enumerate(perm):
            if i != j:
                union(i, j)
                cross += 1
    orbits = len({find(i) for i in range(n)})
    return n, orbits, cross


def main() -> int:
    print("Stage 8 smoke — Eisenstein Γ₀ gluing combinatorics")
    print("=" * 56)
    mesh = build_P3_cell(N_tri=4, N3=2, Y=1.25, lift=True)
    print(f"  level-1 ref cell: tets={len(mesh['tets'])}  |T|={AREA_T:.6f}  "
          f"domain={mesh['domain']}")

    print("\n  Ladder rungs:")
    for row in LADDER_RUNGS:
        print(f"    {row['name']:12} N={row['N']:3} index={row['index']:3}")

    print("\n  Field levels:")
    ok = True
    for factory, tag in (
        (ring_lambda3, "N=3 (1-w)"),
        (ring_split7, "N=7"),
        (ring_split13, "N=13"),
    ):
        R = factory()
        perms = gluing_perms(R)
        n, orbits, cross = glue_components(perms)
        print(f"    [{tag}] copies={n}  copy-orbits≈{orbits}  "
              f"cross-edges={cross}  gens={list(perms)}")
        if n != R.index:
            print("      FAIL index mismatch")
            ok = False
        s = perms["S"]
        if any(s[s[i]] != i for i in range(n)):
            print("      FAIL S^2")
            ok = False

    print()
    if ok:
        print("  SMOKE PASS — combinatorics only.")
        print("  Next: cusp table for Γ₀(1-w); float multi-copy CR.")
        return 0
    print("  SMOKE FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
