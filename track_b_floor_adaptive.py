#!/usr/bin/env python3
"""Rigorous adaptive replacement of validated Track-B floor cells.

The input ledger is a complete uniform Arb audit.  Selected leaf cells are
removed and replaced by a disjoint Cartesian partition into validated child
cells.  The output contains a complete parent-child event ledger and a leaf
ledger; aggregation therefore counts every point of the original grid once.

This utility is intentionally local and fail-closed.  It never promotes the
global Track-B or rung-4 theorem gates.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import math
from pathlib import Path
import time
from typing import Any

from flint import arb, ctx

from continuum_defect_arb import parse_trial
from track_b_direct_weighted_arb import (
    DifferentialTrial,
    _floor_cell_worker,
    _initialize_floor_worker,
    _serialized_floor_cell,
    certified_allowed_floor_budget,
    floor_cell_taylor_model,
)
from track_b_projected_mass_arb import projected_coefficients


ROOT = Path(__file__).resolve().parent
AXIS_FIELDS = ("x1_interval", "x2_interval", "s_interval")
AXIS_NAMES = ("x1", "x2", "s")


def _load_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()]


def _bounds(record: dict[str, Any], axis: int) -> tuple[arb, arb]:
    if "_adaptive_bounds" in record:
        lo, hi = record["_adaptive_bounds"][axis]
        return arb(lo), arb(hi)
    names = ("x1_bounds", "x2_bounds", "s_bounds")
    if names[axis] in record:
        lo, hi = record[names[axis]]
        return arb(lo), arb(hi)
    raise ValueError("cell endpoints absent; reconstruct the base grid before refinement")


def _leaf_id(record: dict[str, Any]) -> str:
    if "adaptive_id" in record:
        return str(record["adaptive_id"])
    return "base:" + ",".join(str(q) for q in record["cell_index"])


def _summary(records: list[dict[str, Any]], budget: arb) -> dict[str, Any]:
    polynomial_sq = arb(0)
    remainder_sq = arb(0)
    combined_sq = arb(0)
    for record in sorted(records, key=_leaf_id):
        polynomial_sq += arb(record["polynomial_l2_squared_upper"]).upper()
        rem = arb(record["remainder_l2_upper"]).upper()
        combined = arb(record["combined_cell_l2_upper"]).upper()
        remainder_sq += rem * rem
        combined_sq += combined * combined
    polynomial = polynomial_sq.sqrt().upper()
    remainder = remainder_sq.sqrt().upper()
    combined = combined_sq.sqrt().upper()
    return {
        "polynomial_l2_upper": str(polynomial),
        "remainder_l2_upper": str(remainder),
        "floor_l2_upper": str(combined),
        "allowed_budget_lower": str(budget.lower()),
        "certified_margin_lower": str((budget.lower() - combined).lower()),
        "certified_ratio_upper": str((combined / budget.lower()).upper()),
        "leaf_count": len(records),
    }


def _choose(
    leaves: list[dict[str, Any]],
    top_fraction: float,
    contribution_fraction: float,
    threshold: arb | None,
) -> list[dict[str, Any]]:
    ranked = sorted(
        leaves,
        key=lambda q: float(arb(q["remainder_l2_upper"]).upper() ** 2),
        reverse=True,
    )
    selected: dict[str, dict[str, Any]] = {}
    if top_fraction > 0:
        count = max(1, math.ceil(top_fraction * len(ranked)))
        for record in ranked[:count]:
            selected[_leaf_id(record)] = record
    if contribution_fraction > 0:
        total = sum(
            (arb(q["remainder_l2_upper"]).upper() ** 2 for q in ranked), arb(0)
        ).upper()
        target = arb(str(contribution_fraction)) * total
        partial = arb(0)
        for record in ranked:
            selected[_leaf_id(record)] = record
            partial += arb(record["remainder_l2_upper"]).upper() ** 2
            if bool(partial >= target):
                break
    if threshold is not None:
        for record in ranked:
            if bool(arb(record["remainder_l2_upper"]).upper() >= threshold):
                selected[_leaf_id(record)] = record
    return sorted(selected.values(), key=_leaf_id)


def _child_specs(parent: dict[str, Any], axes: tuple[int, ...]) -> list[dict[str, Any]]:
    pieces: list[list[tuple[arb, arb, int]]] = []
    for axis in range(3):
        lo, hi = _bounds(parent, axis)
        if axis in axes:
            mid = (lo + hi) / 2
            pieces.append([(lo, mid, 0), (mid, hi, 1)])
        else:
            pieces.append([(lo, hi, 0)])
    out = []
    parent_id = _leaf_id(parent)
    for x in pieces[0]:
        for t in pieces[1]:
            for s in pieces[2]:
                suffix = "".join(
                    f"{AXIS_NAMES[axis]}{piece[2]}" for axis, piece in enumerate((x, t, s))
                    if axis in axes
                )
                out.append({
                    "parent_id": parent_id,
                    "adaptive_id": parent_id + "/" + suffix,
                    "cell_index": parent["cell_index"],
                    "bounds": ((x[0], x[1]), (t[0], t[1]), (s[0], s[1])),
                })
    return out


def _compute_children(
    specs: list[dict[str, Any]], trial: Path, width: arb, degree: int,
    bits: int, workers: int, ev: DifferentialTrial, odd: list[Any],
) -> list[dict[str, Any]]:
    payloads = []
    for job, spec in enumerate(specs):
        xb, tb, sb = spec["bounds"]
        payloads.append((
            (job, 0, 0),
            (str(xb[0]), str(xb[1])),
            (str(tb[0]), str(tb[1])),
            (str(sb[0]), str(sb[1])),
            str(width), degree,
        ))
    if workers > 1:
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=_initialize_floor_worker,
            initargs=(str(trial.resolve()), bits),
        ) as pool:
            children = list(pool.map(_floor_cell_worker, payloads))
    else:
        children = []
        for payload in payloads:
            index, xb, tb, sb, width_text, model_degree = payload
            xbounds = (arb(xb[0]), arb(xb[1]))
            tbounds = (arb(tb[0]), arb(tb[1]))
            sbounds = (arb(sb[0]), arb(sb[1]))
            model = floor_cell_taylor_model(
                ev, odd, xbounds, tbounds, sbounds, arb(width_text), model_degree
            )
            children.append(_serialized_floor_cell(
                index, xbounds, tbounds, sbounds, model_degree, bits, model
            ))
    for child, spec in zip(children, specs):
        child["cell_index"] = spec["cell_index"]
        child["adaptive_id"] = spec["adaptive_id"]
        child["parent_id"] = spec["parent_id"]
        child["_adaptive_bounds"] = [
            [str(lo), str(hi)] for lo, hi in spec["bounds"]
        ]
    return children


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adaptive-floor", action="store_true")
    parser.add_argument("--base-audit", type=Path, required=True)
    parser.add_argument("--trial", type=Path,
                        default=ROOT / "six_copy_hejhal_balanced_coeffs.json")
    parser.add_argument("--mass-certificate", type=Path,
                        default=ROOT / "track_b_projected_mass_arb_result.json")
    parser.add_argument("--width", default="0.30")
    parser.add_argument("--spectral-width-tol", default="0.1")
    parser.add_argument("--degree", type=int, default=8)
    parser.add_argument("--base-grid", default="8,4,8")
    parser.add_argument("--bits", type=int, default=192)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--adaptive-top-fraction", type=float, default=0.0)
    parser.add_argument("--adaptive-contribution-fraction", type=float, default=0.0)
    parser.add_argument("--adaptive-max-depth", type=int, default=1)
    parser.add_argument("--adaptive-axes", default="x1,x2,s")
    parser.add_argument("--adaptive-remainder-threshold", default="")
    parser.add_argument("--adaptive-target-remainder", default="0.004")
    parser.add_argument("--audit-jsonl", type=Path,
                        default=ROOT / "track_b_floor_adaptive_leaves.jsonl")
    parser.add_argument("--parent-child-jsonl", type=Path,
                        default=ROOT / "track_b_floor_adaptive_tree.jsonl")
    parser.add_argument("--json-out", type=Path,
                        default=ROOT / "track_b_floor_adaptive_result.json")
    ns = parser.parse_args()
    if not ns.adaptive_floor:
        parser.error("--adaptive-floor is required (fail-closed explicit mode)")
    if ns.workers <= 0 or ns.adaptive_max_depth < 0:
        parser.error("workers must be positive and max depth nonnegative")
    if not (0 <= ns.adaptive_top_fraction <= 1):
        parser.error("top fraction must lie in [0,1]")
    if not (0 <= ns.adaptive_contribution_fraction <= 1):
        parser.error("contribution fraction must lie in [0,1]")
    axes_by_name = {name: axis for axis, name in enumerate(AXIS_NAMES)}
    try:
        axes = tuple(sorted({axes_by_name[q.strip()] for q in ns.adaptive_axes.split(",")
                             if q.strip()}))
    except KeyError as error:
        parser.error(f"unknown refinement axis {error.args[0]!r}")
    if not axes:
        parser.error("at least one refinement axis is required")

    ctx.prec = max(128, ns.bits)
    width = arb(ns.width)
    mass = json.loads(ns.mass_certificate.read_text(encoding="utf-8"))
    budget = certified_allowed_floor_budget(mass, arb(ns.spectral_width_tol))
    data, row = parse_trial(ns.trial)
    ev = DifferentialTrial(data["parameters"]["M"], str(row["r"]), row["coefficients"])
    odd = projected_coefficients(ev)["odd"]
    leaves = _load_jsonl(ns.base_audit)
    base_grid = tuple(int(q) for q in ns.base_grid.split(","))
    if len(base_grid) != 3 or min(base_grid) <= 0:
        parser.error("base grid must be nx1,nx2,ns")
    for record in leaves:
        record.setdefault("adaptive_id", _leaf_id(record))
        if int(record["polynomial_degree"]) != ns.degree:
            raise ValueError("base ledger degree does not match --degree")
        if not all(name in record for name in ("x1_bounds", "x2_bounds", "s_bounds")):
            i, j, k = (int(q) for q in record["cell_index"])
            nx, nt, n_s = base_grid
            record["_adaptive_bounds"] = [
                [str(-arb(1) / 2 + arb(i) / nx),
                 str(-arb(1) / 2 + arb(i + 1) / nx)],
                [str(arb(j) / (2 * nt)), str(arb(j + 1) / (2 * nt))],
                [str(width * k / n_s), str(width * (k + 1) / n_s)],
            ]
    base_summary = _summary(leaves, budget)
    target = arb(ns.adaptive_target_remainder)
    threshold = (arb(ns.adaptive_remainder_threshold)
                 if ns.adaptive_remainder_threshold.strip() else None)
    events = []
    started = time.perf_counter()
    for depth in range(1, ns.adaptive_max_depth + 1):
        current = _summary(leaves, budget)
        if bool(arb(current["remainder_l2_upper"]).upper() <= target):
            break
        selected = _choose(
            leaves, ns.adaptive_top_fraction,
            ns.adaptive_contribution_fraction, threshold,
        )
        if not selected:
            break
        specs = []
        for parent in selected:
            specs.extend(_child_specs(parent, axes))
        children = _compute_children(
            specs, ns.trial, width, ns.degree, int(ctx.prec), ns.workers, ev, odd
        )
        selected_ids = {_leaf_id(q) for q in selected}
        kept = [q for q in leaves if _leaf_id(q) not in selected_ids]
        if len(kept) != len(leaves) - len(selected):
            raise AssertionError("parent replacement was not one-for-one")
        child_ids = [_leaf_id(q) for q in children]
        if len(child_ids) != len(set(child_ids)):
            raise AssertionError("duplicate child leaf identifier")
        events.extend({
            "depth": depth,
            "parent_id": _leaf_id(parent),
            "parent_remainder_l2_upper": parent["remainder_l2_upper"],
            "children": [_leaf_id(q) for q in children
                         if q["parent_id"] == _leaf_id(parent)],
        } for parent in selected)
        leaves = kept + children
        all_ids = [_leaf_id(q) for q in leaves]
        if len(all_ids) != len(set(all_ids)):
            raise AssertionError("leaf partition contains a duplicate")

    elapsed = time.perf_counter() - started
    final_summary = _summary(leaves, budget)
    ns.audit_jsonl.write_text(
        "".join(json.dumps(q, sort_keys=True) + "\n" for q in sorted(leaves, key=_leaf_id)),
        encoding="utf-8",
    )
    ns.parent_child_jsonl.write_text(
        "".join(json.dumps(q, sort_keys=True) + "\n" for q in events),
        encoding="utf-8",
    )
    result = {
        "schema": "track-b-floor-adaptive/v1",
        "base_audit": str(ns.base_audit.resolve()),
        "polynomial_degree": ns.degree,
        "arb_bits": int(ctx.prec),
        "base": base_summary,
        "adaptive": final_summary,
        "refinement_axes": [AXIS_NAMES[q] for q in axes],
        "maximum_depth": ns.adaptive_max_depth,
        "top_fraction": ns.adaptive_top_fraction,
        "squared_contribution_fraction": ns.adaptive_contribution_fraction,
        "remainder_threshold": ns.adaptive_remainder_threshold or None,
        "target_global_remainder": ns.adaptive_target_remainder,
        "parents_replaced": len(events),
        "parent_child_events": len(events),
        "elapsed_seconds": elapsed,
        "leaf_audit_jsonl": str(ns.audit_jsonl.resolve()),
        "parent_child_audit_jsonl": str(ns.parent_child_jsonl.resolve()),
        "partition_exact_once_asserted": True,
        "floor_residual_certified": False,
        "global_partition_certified": False,
        "global_weight_bounds_certified": False,
        "rung4_integrator_comparison_certified": False,
        "rung4_certified": False,
        "status": "DIAGNOSTIC_STABILITY_OR_GLOBAL_GATES_PENDING",
    }
    ns.json_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
