#!/usr/bin/env python3
"""Recompute selected floor cells with Taylor arithmetic attribution enabled."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from flint import arb, ctx

from continuum_defect_arb import parse_trial
from track_b_direct_weighted_arb import DifferentialTrial, floor_cell_taylor_model
from track_b_projected_mass_arb import projected_coefficients


ROOT = Path(__file__).resolve().parent


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trial", type=Path,
                        default=ROOT / "six_copy_hejhal_balanced_coeffs.json")
    parser.add_argument("--grid", default="8,4,8")
    parser.add_argument("--cell", required=True)
    parser.add_argument("--degree", type=int, default=0)
    parser.add_argument("--tensor-degree", default="")
    parser.add_argument("--width", default="0.30")
    parser.add_argument("--bits", type=int, default=192)
    parser.add_argument("--arithmetic-audit", action="store_true")
    parser.add_argument("--json-out", type=Path, default=None)
    ns = parser.parse_args()
    ctx.prec = max(128, ns.bits)
    grid = tuple(int(q) for q in ns.grid.split(","))
    cell = tuple(int(q) for q in ns.cell.split(","))
    if len(grid) != 3 or len(cell) != 3:
        parser.error("grid and cell require three comma-separated integers")
    nx, nt, n_s = grid
    i, j, k = cell
    width = arb(ns.width)
    bounds = (
        (-arb(1) / 2 + arb(i) / nx, -arb(1) / 2 + arb(i + 1) / nx),
        (arb(j) / (2 * nt), arb(j + 1) / (2 * nt)),
        (width * k / n_s, width * (k + 1) / n_s),
    )
    data, row = parse_trial(ns.trial)
    ev = DifferentialTrial(data["parameters"]["M"], str(row["r"]), row["coefficients"])
    odd = projected_coefficients(ev)["odd"]
    degree: int | tuple[int, int, int]
    if ns.tensor_degree.strip():
        degree = tuple(int(q) for q in ns.tensor_degree.split(","))
        if len(degree) != 3:
            parser.error("tensor degree must be px,py,ps")
    elif ns.degree:
        degree = ns.degree
    else:
        parser.error("provide --degree or --tensor-degree")
    result = floor_cell_taylor_model(
        ev, odd, *bounds, width, degree,
        collect_arithmetic_audit=ns.arithmetic_audit,
    )
    output = {
        "grid": list(grid),
        "cell": list(cell),
        "degree": list(degree) if isinstance(degree, tuple) else degree,
        "bounds": [[str(lo), str(hi)] for lo, hi in bounds],
        "polynomial_l2_upper": str(result["polynomial_l2_upper_ball"]),
        "remainder_l2_upper": str(result["remainder_l2_upper_ball"]),
        "combined_cell_l2_upper": str(result["combined_cell_l2_upper_ball"]),
        "bessel_direct_count": result["bessel_direct_count"],
        "bessel_majorant_count": result["bessel_majorant_count"],
        "bessel_fallback_count": result["bessel_fallback_count"],
        "bessel_real_order_counts": result["bessel_real_order_counts"],
        "arithmetic_audit": result["arithmetic_audit"],
    }
    text = json.dumps(output, indent=2) + "\n"
    if ns.json_out is not None:
        ns.json_out.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
