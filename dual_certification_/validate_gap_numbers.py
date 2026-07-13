#!/usr/bin/env python3
"""Independent validation of Rung 4 quantitative gap numbers. Non-certifying."""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lemma_K import C1_C2_constants, ball_mid, lemma_K_tail  # noqa: E402


def main() -> int:
    eta_scan = 5.04e-4
    print("=== Independent C1 / eta0 / gap verification ===")
    print(f"scan residual eta = {eta_scan:.6e}\n")
    rows = []
    for label, kw in [
        ("cons r=6 Y0=0.8", dict(sharp_geom=False, r=6.0, Y0=0.8)),
        # Pre-1.1/1.2 sharp freeze (y_min=1/2, D_Y=3) for quote continuity
        (
            "sharp_legacy r=6 Y0=0.8",
            dict(sharp_geom=True, r=6.0, Y0=0.8, y_min=0.5, D_Y=3.0),
        ),
        # Current default: y_min=1/√2 + box D_Y (morningbrief 1.1–1.2)
        ("sharp r=6 Y0=0.8", dict(sharp_geom=True, r=6.0, Y0=0.8)),
        ("sharp r=6 Y0=1.5", dict(sharp_geom=True, r=6.0, Y0=1.5)),
        ("sharp U=2 Y0=1.5", dict(sharp_geom=True, r=6.0, Y0=1.5, U_norm=2.0)),
        ("sharp r=6.622 Y0=0.8", dict(sharp_geom=True, r=6.62212, Y0=0.8)),
    ]:
        d = C1_C2_constants(field="i", **kw)
        C1 = ball_mid(d["C1"])
        C2 = ball_mid(d["C2"])
        e0 = ball_mid(d["eta0"])
        y0 = kw["Y0"]
        collar = C2 * math.exp(-2 * math.pi * y0)
        dlam = C1 * eta_scan + collar
        gap = eta_scan / e0
        orders = math.log10(gap)
        eta_w = max(0.0, (0.1 - collar) / C1) if C1 > 0 else float("nan")
        rows.append((label, C1, C2, e0, gap, orders, dlam, collar, eta_w, d))
        print(label)
        print(f"  C1={C1:.6e}  C2={C2:.6e}  eta0={e0:.6e}")
        print(
            f"  C_tr={ball_mid(d['C_trace']):.4f}  D_Y={ball_mid(d['D_Y']):.4f}  "
            f"y_min={ball_mid(d['y_min']):.4f}"
        )
        print(f"  eta/eta0={gap:.3e}  orders={orders:.2f}")
        print(f"  dlam@eta_scan={dlam:.6e}  collar={collar:.6e}")
        print(f"  eta for width<0.1 ≈ {eta_w:.3e}")
        print()

    print("=== Cross-check user quoted numbers ===")
    # User: before sharp C1 1.02e6 eta0 2.4e-13; after C1 7.3e4 eta0 4.6e-11
    c0, e0c = rows[0][1], rows[0][3]
    c_leg, e_leg = rows[1][1], rows[1][3]
    c1, e01 = rows[2][1], rows[2][3]
    print(f"conservative C1={c0:.4e} (user ~1.02e6) rel={c0/1.02e6:.3f}")
    print(f"conservative eta0={e0c:.4e} (user ~2.4e-13) rel={e0c/2.4e-13:.3f}")
    print(f"sharp_legacy C1={c_leg:.4e} (user ~7.3e4) rel={c_leg/7.3e4:.3f}")
    print(f"sharp_legacy eta0={e_leg:.4e} (user ~4.6e-11) rel={e_leg/4.6e-11:.3f}")
    print(f"sharp_default (1.1+1.2) C1={c1:.4e}  eta0={e01:.4e}")
    print(
        f"  factor vs legacy C1 = {c_leg/c1:.2f}×  "
        f"orders gap now={rows[2][5]:.2f} (was {rows[1][5]:.2f})"
    )
    print(f"gap cons orders={rows[0][5]:.2f} (user ~9)")
    print(f"gap legacy sharp orders={rows[1][5]:.2f} (user ~7)")

    print("\n=== kappa summary file ===")
    p = Path(__file__).resolve().parent / "hejhal_kappa_summary.txt"
    print(p.read_text(encoding="utf-8"))

    print("=== Lemma K tail M=400 r=6.62212 Y0=0.8 ===")
    t = lemma_K_tail(400, 0.8, 6.62212, 0.5)
    print(t)
    print("\nVALIDATION OK if relatives near 1.0 for cons + sharp_legacy quotes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
