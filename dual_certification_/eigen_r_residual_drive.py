#!/usr/bin/env python3
"""
Drive residual: coarse r-grid + golden-section refine on two-cusp η_proxy.

Uses production_hejhal_residual.residual_at_r (two-cusp near-kernel +
PAIRINGS δ_aut, optional periodize). Reports best r and factor vs baseline.

Language: engineering residual drive, not certified Maass / dual GREEN.

Usage:
  python eigen_r_residual_drive.py
  python eigen_r_residual_drive.py --M 32 --r-min 5.0 --r-max 8.0 --n-grid 9
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from production_hejhal_residual import residual_at_r  # noqa: E402


def drive(
    M: int = 32,
    r_min: float = 5.0,
    r_max: float = 8.0,
    n_grid: int = 9,
    Y0: float = 0.8,
    periodize: bool = True,
    refine_iters: int = 10,
) -> Dict[str, Any]:
    t0 = time.time()
    # coarse grid on relative residual of two-cusp operator (fast proxy)
    # then refine on full η = δ_aut + τ_proxy
    print(f"=== eigen-r residual drive  M={M} Y0={Y0} periodize={periodize} ===")
    print(f"coarse grid r∈[{r_min},{r_max}] n={n_grid}")

    rs = [r_min + (r_max - r_min) * i / max(n_grid - 1, 1) for i in range(n_grid)]
    rows: List[Dict[str, Any]] = []
    best: Optional[Dict[str, Any]] = None
    for r in rs:
        print(f"  grid r={r:.5f} ...", flush=True)
        rec = residual_at_r(M, float(r), Y0=Y0, periodize=periodize)
        row = dict(
            r=float(r),
            delta_aut=float(rec["cusp_inf"]["delta_aut"]),
            delta_0=float(rec["cusp_0"]["delta_aut"]),
            rel=float(rec["two_cusp"]["rel"]),
            tau_proxy=float(rec["tau_proxy"]),
            eta_proxy=float(rec["eta_proxy"]),
            kappa_eq=float(rec["two_cusp"]["kappa_eq"]),
            seconds=float(rec["seconds"]),
            phase="grid",
        )
        rows.append(row)
        print(
            f"    δ={row['delta_aut']:.4e}  rel={row['rel']:.3e}  "
            f"η={row['eta_proxy']:.4e}",
            flush=True,
        )
        if best is None or row["eta_proxy"] < best["eta_proxy"]:
            best = dict(row)
            best["full"] = rec

    assert best is not None
    # golden-section refine around best grid point
    span = (r_max - r_min) / max(n_grid - 1, 1)
    a = max(r_min, best["r"] - span)
    b = min(r_max, best["r"] + span)
    phi = (math.sqrt(5) - 1) / 2
    print(f"\ngolden refine on [{a:.5f},{b:.5f}] iters={refine_iters}")
    for it in range(refine_iters):
        r1 = b - phi * (b - a)
        r2 = a + phi * (b - a)
        rec1 = residual_at_r(M, r1, Y0=Y0, periodize=periodize)
        rec2 = residual_at_r(M, r2, Y0=Y0, periodize=periodize)
        e1, e2 = rec1["eta_proxy"], rec2["eta_proxy"]
        print(
            f"  it={it}  r1={r1:.5f} η1={e1:.4e}  r2={r2:.5f} η2={e2:.4e}",
            flush=True,
        )
        for r_try, rec in ((r1, rec1), (r2, rec2)):
            row = dict(
                r=float(r_try),
                delta_aut=float(rec["cusp_inf"]["delta_aut"]),
                delta_0=float(rec["cusp_0"]["delta_aut"]),
                rel=float(rec["two_cusp"]["rel"]),
                tau_proxy=float(rec["tau_proxy"]),
                eta_proxy=float(rec["eta_proxy"]),
                kappa_eq=float(rec["two_cusp"]["kappa_eq"]),
                seconds=float(rec["seconds"]),
                phase=f"refine_{it}",
            )
            rows.append(row)
            if row["eta_proxy"] < best["eta_proxy"]:
                best = dict(row)
                best["full"] = rec
        if e1 < e2:
            b = r2
        else:
            a = r1

    baseline = rows[0]
    fac = baseline["eta_proxy"] / best["eta_proxy"] if best["eta_proxy"] > 0 else float("inf")
    out = dict(
        M=M,
        Y0=Y0,
        periodize=periodize,
        r_min=r_min,
        r_max=r_max,
        n_grid=n_grid,
        refine_iters=refine_iters,
        rows=[{k: v for k, v in row.items() if k != "full"} for row in rows],
        best=dict(
            r=best["r"],
            delta_aut=best["delta_aut"],
            delta_0=best["delta_0"],
            rel=best["rel"],
            tau_proxy=best["tau_proxy"],
            eta_proxy=best["eta_proxy"],
            kappa_eq=best["kappa_eq"],
            phase=best.get("phase", ""),
            per_generator=best.get("full", {})
            .get("cusp_inf", {})
            .get("per_generator", {}),
        ),
        baseline_eta=baseline["eta_proxy"],
        factor_vs_first_grid=fac,
        seconds=time.time() - t0,
        language=(
            "Eigen-r residual drive on two-cusp model + PAIRINGS δ_aut. "
            "Not certified Maass form. Hard map unchanged."
        ),
    )
    print(
        f"\nBEST r≈{out['best']['r']:.6f}  δ={out['best']['delta_aut']:.4e}  "
        f"η={out['best']['eta_proxy']:.4e}  "
        f"factor vs first grid={fac:.2f}×  ({out['seconds']:.1f}s)"
    )
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--M", type=int, default=32)
    p.add_argument("--r-min", type=float, default=5.0)
    p.add_argument("--r-max", type=float, default=8.0)
    p.add_argument("--n-grid", type=int, default=7)
    p.add_argument("--Y0", type=float, default=0.8)
    p.add_argument("--refine-iters", type=int, default=8)
    p.add_argument("--no-periodize", action="store_true")
    p.add_argument("--json-out", type=str, default="")
    args = p.parse_args(argv)

    out = drive(
        M=args.M,
        r_min=args.r_min,
        r_max=args.r_max,
        n_grid=args.n_grid,
        Y0=args.Y0,
        periodize=not args.no_periodize,
        refine_iters=args.refine_iters,
    )
    path = Path(args.json_out) if args.json_out else _HERE / "eigen_r_drive_result.json"
    path.write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
