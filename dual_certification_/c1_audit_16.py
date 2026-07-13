#!/usr/bin/env python3
"""Measure C1 factor from dropping (1+√λ) in A_bdry (audit §1.6 only)."""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lemma_K import C1_C2_constants, ball_mid  # noqa: E402


def main() -> int:
    eta = 5.04e-4
    print("=== C1 audit §1.6 (bdry_sqrt_lam) ===")
    print("Production default keeps bdry_sqrt_lam=True.\n")
    rows = []
    for label, kw in [
        ("default sharp", dict(sharp_geom=True, r=6.0, Y0=1.5, bdry_sqrt_lam=True)),
        ("audit no (1+√λ)", dict(sharp_geom=True, r=6.0, Y0=1.5, bdry_sqrt_lam=False)),
        ("legacy y=1/2 D_Y=3", dict(sharp_geom=True, r=6.0, Y0=1.5, y_min=0.5, D_Y=3.0)),
    ]:
        d = C1_C2_constants(field="i", **kw)
        C1 = ball_mid(d["C1"])
        e0 = ball_mid(d["eta0"])
        orders = math.log10(eta / e0) if e0 > 0 else float("inf")
        rows.append((label, C1, e0, orders))
        print(f"{label}")
        print(f"  C1={C1:.6e}  eta0={e0:.6e}  orders@eta={orders:.2f}")
        print(f"  A_bdry={ball_mid(d['A_bdry']):.6e}  A_ell={ball_mid(d['A_ell']):.6e}")
        print(f"  bdry_sqrt_lam={d.get('bdry_sqrt_lam')}")
        print()
    c_def, c_aud = rows[0][1], rows[1][1]
    print(f"factor default/audit C1 = {c_def/c_aud:.3f}×")
    print(f"1+√37 = {1+math.sqrt(37):.6f}  (expected ≈ factor)")
    print("See C1_AUDIT_1_6.md — do not flip hard map / production eta0.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
