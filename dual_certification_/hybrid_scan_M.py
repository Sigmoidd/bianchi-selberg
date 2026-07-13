#!/usr/bin/env python3
"""
Larger-M hybrid / periodize δ_aut scan (morningbrief priority after §1.6).

For each M in the list, measure δ_aut under:
  collocation | multi (hybrid) | periodize

Writes hybrid_scan_M.csv and hybrid_scan_M.md.

Usage:
  python hybrid_scan_M.py
  python hybrid_scan_M.py --M 64,100,200,400 --jump-weight 2
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from delta_aut_pairing import delta_aut_on_pairings  # noqa: E402


def run_scan(
    Ms: Sequence[int],
    r: float = 6.0,
    Y0: float = 0.8,
    n_face: int = 16,
    jump_weight: float = 2.0,
    modes: Sequence[str] = ("collocation", "multi", "periodize"),
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    t_all = time.time()
    for M in Ms:
        print(f"\n=== M={M} ===", flush=True)
        base_delta = None
        for mode in modes:
            print(f"  mode={mode} ...", flush=True)
            out = delta_aut_on_pairings(
                M=M,
                r=r,
                Y0=Y0,
                n_face=n_face,
                mode=mode,
                jump_weight=jump_weight,
            )
            d = float(out["delta_aut"])
            if mode == "collocation":
                base_delta = d
            fac = (base_delta / d) if (base_delta and d > 0) else float("nan")
            rec = dict(
                M=M,
                mode=mode,
                r=r,
                Y0=Y0,
                jump_weight=jump_weight,
                delta_aut=d,
                tau_proxy=float(out["tau_proxy"]),
                eta_proxy=float(out["eta_production_proxy"]),
                rel=float(out["collocation"]["rel"]),
                kappa_eq=float(out["collocation"]["kappa_eq"]),
                n_modes=int(out["collocation"]["n_modes"]),
                n_pts=int(out["collocation"]["n_pts"]),
                factor_vs_col=fac,
                seconds=float(out["seconds"]),
            )
            # per-gen maxima
            for gname, grec in out["per_generator"].items():
                rec[f"d_{gname}"] = float(grec["delta_max"])
            rows.append(rec)
            print(
                f"    δ={d:.4e}  τ={rec['tau_proxy']:.3e}  "
                f"η={rec['eta_proxy']:.3e}  κ={rec['kappa_eq']:.1f}  "
                f"fac={fac:.2f}×  ({rec['seconds']:.1f}s)",
                flush=True,
            )
    print(f"\nTotal wall time: {time.time() - t_all:.1f}s")
    return rows


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--M", type=str, default="64,100,200,400")
    p.add_argument("--r", type=float, default=6.0)
    p.add_argument("--Y0", type=float, default=0.8)
    p.add_argument("--n-face", type=int, default=16)
    p.add_argument("--jump-weight", type=float, default=2.0)
    p.add_argument(
        "--modes",
        type=str,
        default="collocation,multi,periodize",
        help="comma list of modes",
    )
    args = p.parse_args(argv)

    Ms = [int(x.strip()) for x in args.M.split(",") if x.strip()]
    modes = [x.strip() for x in args.modes.split(",") if x.strip()]
    print("=== Hybrid / periodize δ_aut scan vs M ===")
    print(f"M={Ms}  r={args.r}  Y0={args.Y0}  jump_w={args.jump_weight}")
    rows = run_scan(
        Ms,
        r=args.r,
        Y0=args.Y0,
        n_face=args.n_face,
        jump_weight=args.jump_weight,
        modes=modes,
    )

    csv_path = _HERE / "hybrid_scan_M.csv"
    fields = list(rows[0].keys()) if rows else []
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for rec in rows:
            w.writerow(rec)
    print(f"wrote {csv_path}")

    md = _HERE / "hybrid_scan_M.md"
    lines = [
        "# Hybrid / periodize δ_aut vs M",
        "",
        f"r={args.r}, Y0={args.Y0}, jump_weight={args.jump_weight}, n_face={args.n_face}",
        "",
        "| M | mode | δ_aut | τ_proxy | η_proxy | κ_eq | fac vs col | s |",
        "|--:|------|------:|--------:|--------:|-----:|-----------:|--:|",
    ]
    for rec in rows:
        lines.append(
            f"| {rec['M']} | {rec['mode']} | {rec['delta_aut']:.3e} | "
            f"{rec['tau_proxy']:.3e} | {rec['eta_proxy']:.3e} | "
            f"{rec['kappa_eq']:.1f} | {rec['factor_vs_col']:.2f} | "
            f"{rec['seconds']:.1f} |"
        )
    lines += [
        "",
        "Non-certifying. Hard map unchanged.",
        "fac vs col = δ_collocation / δ_mode at same M.",
    ]
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
