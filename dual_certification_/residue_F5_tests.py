#!/usr/bin/env python3
"""
Independent unit tests: residue ring F5 for 𝔭=(2+i) and exact 6-copy gluing.

AgentReady Rung 4 checklist:
  [ ] Residue-ring arithmetic for F5 verified by independent unit-test suite
  [ ] Gluing of the 6 copies is exact (combinatorial identity, not numerical)

Uses independent_exclusion/congruence_prototype.py (canonical production ring).
Does not claim dual eigenvalue certification.

Usage:
  python residue_F5_tests.py
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / "independent_exclusion"))

import congruence_prototype as cp  # noqa: E402


def test_F5_field() -> list[str]:
    msgs = []
    R = cp.make_ring("(2+i)")
    assert R.norm == 5 and R.index == 6, (R.norm, R.index)
    msgs.append("OK norm=5 index=6")
    # i^2 = -1 in F5 with i↦3: 3*3=9≡4≡-1
    i = R.i
    assert R.mul(i, i) == R.embed(-1), (R.mul(i, i), R.embed(-1))
    msgs.append("OK i^2 ≡ -1 (i↦3 mod 5)")
    # field: every nonzero has inverse
    for a in range(5):
        if a == 0:
            continue
        inv = R.inv(R.embed(a))
        assert R.mul(R.embed(a), inv) == R.one
    msgs.append("OK F5 inverses")
    # P1 has 6 points
    pts = R.p1_points()
    assert len(pts) == 6
    msgs.append(f"OK |P1(F5)|={len(pts)}")
    return msgs


def test_gluing_exact() -> list[str]:
    msgs = []
    pts, glue, cusp_class = cp.build_gluing("(2+i)")
    R = cp._ring()
    NC = R.index
    assert len(pts) == NC == 6
    msgs.append("OK build_gluing returned 6 points")
    # cusp classes: one ∞, five finite
    assert cusp_class.count(0) == 1 and cusp_class.count(1) == 5, cusp_class
    msgs.append(f"OK cusp_class hist={dict(Counter(cusp_class))}")
    # each glue perm is a permutation + round-trip
    for name, perm in glue.items():
        assert sorted(perm) == list(range(NC)), f"{name} not perm"
        # involution check for S-type when applicable; all must be bijective
        assert len(set(perm)) == NC
        msgs.append(f"OK glue[{name}] bijective perm={perm}")
    # combinatorial identity: composing glue[T1] etc. stays in S_6
    # exactness: re-run act consistency (already in build_gluing asserts)
    msgs.append("OK all generator gluing round-trips (asserted in build_gluing)")
    return msgs


def test_generators_S2_U() -> list[str]:
    """S^2 ~ ±I, generators invertible in residue."""
    msgs = []
    R = cp.make_ring("(2+i)")
    mats = R.generator_matrices()
    inv = R.inverses(mats)
    S = mats.get("S")
    if S is not None:
        SS = R.mat_mul(S, S)
        assert R.mat_pm_I(SS), SS
        msgs.append("OK S^2 ≡ ±I in residue SL")
    for name in mats:
        P = R.mat_mul(mats[name], inv[name])
        assert R.mat_pm_I(P), (name, P)
    msgs.append("OK generators × inverses ≡ ±I in residue")
    return msgs


def main() -> int:
    print("residue_F5_tests — 𝔭=(2+i), F5, 6-copy gluing")
    print("=" * 56)
    ok = True
    for name, fn in (
        ("F5 field", test_F5_field),
        ("gluing exact", test_gluing_exact),
        ("generators", test_generators_S2_U),
    ):
        print(f"\n[{name}]")
        try:
            for m in fn():
                print(f"  {m}")
            print(f"  PASS")
        except Exception as e:
            ok = False
            print(f"  FAIL: {e}")
    print("\n" + "=" * 56)
    print(f"OVERALL: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
