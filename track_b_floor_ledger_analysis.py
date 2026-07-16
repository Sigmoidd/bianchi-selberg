#!/usr/bin/env python3
"""Independent concentration audit for validated Track-B floor cell ledgers.

The script performs no midpoint recomputation.  Every ranking weight is the
square of the rigorous ``remainder_l2_upper`` already written by the Arb run.
Geometry labels are exact incidences in the normalized coordinates
``(x1,x2,s=log rho)`` for the closed width-0.30 Humbert collar.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from flint import arb, ctx


SPECIAL = {
    "v_00": (arb(0), arb(0), arb(0)),
    "v_0h": (arb(0), arb(1) / 2, arb(0)),
    "v_mh0": (-arb(1) / 2, arb(0), arb(0)),
    "v_ph0": (arb(1) / 2, arb(0), arb(0)),
    "v_mhh": (-arb(1) / 2, arb(1) / 2, arb(0)),
    "v_phh": (arb(1) / 2, arb(1) / 2, arb(0)),
}


def _contains(box: arb, point: arb) -> bool:
    return bool(box.contains(point))


def _record_bounds(
    record: dict[str, Any], grid: tuple[int, int, int], width: arb
) -> tuple[tuple[arb, arb], tuple[arb, arb], tuple[arb, arb]]:
    names = ("x1_bounds", "x2_bounds", "s_bounds")
    if all(name in record for name in names):
        return tuple(tuple(arb(q) for q in record[name]) for name in names)  # type: ignore[return-value]
    i, j, k = (int(q) for q in record["cell_index"])
    nx, nt, ns = grid
    return (
        (-arb(1) / 2 + arb(i) / nx, -arb(1) / 2 + arb(i + 1) / nx),
        (arb(j) / (2 * nt), arb(j + 1) / (2 * nt)),
        (width * k / ns, width * (k + 1) / ns),
    )


def _cell_geometry(
    record: dict[str, Any], width: arb, grid: tuple[int, int, int]
) -> dict[str, Any]:
    (xa, xb), (ta, tb), (sa, sb) = _record_bounds(record, grid, width)
    floor = bool(sa.contains(0))
    outer = bool(sb.contains(width))
    def hits(lo: arb, hi: arb, point: arb) -> bool:
        return bool(lo <= point and point <= hi)
    edges = []
    if floor:
        if hits(xa, xb, -arb(1) / 2):
            edges.append("floor_x1m")
        if hits(xa, xb, arb(1) / 2):
            edges.append("floor_x1p")
        if hits(ta, tb, arb(0)):
            edges.append("floor_x2m")
        if hits(ta, tb, arb(1) / 2):
            edges.append("floor_x2p")
    vertices = [
        name for name, (vx, vt, vs) in SPECIAL.items()
        if hits(xa, xb, vx) and hits(ta, tb, vt) and hits(sa, sb, vs)
    ]
    return {
        "exact_bounds": {
            "x1": [str(xa), str(xb)],
            "x2": [str(ta), str(tb)],
            "s": [str(sa), str(sb)],
        },
        "touches_collar_boundary": bool(floor or outer),
        "collar_boundary_components": [
            name for name, hit in (("rho=1", floor), ("log(rho)=0.30", outer))
            if hit
        ],
        "touches_inversion_sphere": floor,
        "elliptic_edges": edges,
        "elliptic_vertices": vertices,
    }


def analyze(path: Path, width: arb, top: int) -> dict[str, Any]:
    records = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
               if line.strip()]
    weighted = []
    total_sq = arb(0)
    for record in records:
        remainder = arb(record["remainder_l2_upper"]).upper()
        contribution = (remainder * remainder).upper()
        total_sq += contribution
        weighted.append((contribution, remainder, record))
    weighted.sort(key=lambda row: float(row[0]), reverse=True)
    total_sq = total_sq.upper()
    grid = tuple(max(int(q[2]["cell_index"][axis]) for q in weighted) + 1
                 for axis in range(3))

    def count_for(fraction: str) -> int:
        target = arb(fraction) * total_sq
        partial = arb(0)
        for index, (contribution, _remainder, _record) in enumerate(weighted, 1):
            partial += contribution
            if bool(partial >= target):
                return index
        return len(weighted)

    rows = []
    cumulative = arb(0)
    for rank, (contribution, remainder, record) in enumerate(weighted[:top], 1):
        cumulative += contribution
        row = {
            "rank": rank,
            "cell_index": record["cell_index"],
            "remainder_l2_upper": record["remainder_l2_upper"],
            "squared_remainder_fraction_upper": str(
                (contribution / total_sq.lower()).upper()
            ),
            "cumulative_squared_fraction_upper": str(
                (cumulative / total_sq.lower()).upper()
            ),
            **_cell_geometry(record, width, grid),
            "dominant_omitted_multi_indices": record.get(
                "dominant_omitted_multi_indices", "not recorded by baseline run"
            ),
        }
        rows.append(row)
    direct = sum(int(q.get("bessel_direct_count", 0)) for q in records)
    fallback = sum(int(q.get("bessel_fallback_count", 0)) for q in records)
    majorant = sum(int(q.get("bessel_majorant_count", 0)) for q in records)
    return {
        "schema": "track-b-floor-ledger-analysis/v1",
        "source": str(path.resolve()),
        "cell_count": len(records),
        "grid_reconstructed_from_cell_indices": list(grid),
        "rigorous_remainder_l2_reconstructed": str(total_sq.sqrt().upper()),
        "cells_for_50_percent_squared_remainder": count_for("0.50"),
        "cells_for_90_percent_squared_remainder": count_for("0.90"),
        "top_cells": rows,
        "bessel_breakdown": {
            "direct_complex_order_evaluations": direct,
            "real_order_majorant_evaluations": majorant,
            "fallback_evaluations": fallback,
            "fourier_tail_upper": "0",
        },
        "geometry_basis": (
            "exact floor-collar ledger: inversion sphere s=0; elliptic edges are "
            "s=0 intersected with x1=+-1/2 or x2=0,1/2"
        ),
        "legacy_interval_field_warning": (
            "baseline *_interval strings used arb.union and are not endpoint "
            "serializations; exact_bounds are reconstructed from the complete grid indices"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--width", default="0.30")
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--json-out", type=Path, default=None)
    parser.add_argument("--bits", type=int, default=192)
    ns = parser.parse_args()
    ctx.prec = max(128, ns.bits)
    result = analyze(ns.ledger, arb(ns.width), ns.top)
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if ns.json_out is not None:
        ns.json_out.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
