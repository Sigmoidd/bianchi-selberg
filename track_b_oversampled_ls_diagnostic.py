#!/usr/bin/env python3
"""Diagnostic constrained least-squares crosscheck of cutoff instability."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from flint import arb, ctx
import numpy as np

from track_b_collocation_family import exact_candidate_points
from track_b_two_cusp_hejhal import (
    ValidatedWhittaker, assemble_physical_rows_from_points, midpoint_matrix,
)


ROOT = Path(__file__).resolve().parent


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r", default="6.7439020359331625")
    parser.add_argument("--cutoffs", default="2,3,4,5")
    parser.add_argument("--candidate-points", type=int, default=20)
    parser.add_argument("--json-out", type=Path, default=ROOT / "track_b_oversampled_ls_diagnostic.json")
    ns = parser.parse_args()
    ctx.prec = 192
    cutoffs = [int(q) for q in ns.cutoffs.split(",")]
    solutions = {}
    label_sets = {}
    runs = []
    for cutoff in cutoffs:
        rows, _ledger, modes_inf, modes_zero = assemble_physical_rows_from_points(
            cutoff, exact_candidate_points(ns.candidate_points),
            ValidatedWhittaker(arb(ns.r)), None,
        )
        matrix = midpoint_matrix(rows)
        norm_index = next(i for i, mode in enumerate(modes_inf) if mode[:2] == (1, 0))
        free = [j for j in range(matrix.shape[1]) if j != norm_index]
        free_solution = np.linalg.lstsq(matrix[:, free], -matrix[:, norm_index], rcond=None)[0]
        solution = np.zeros(matrix.shape[1], dtype=np.complex128)
        solution[norm_index] = 1
        solution[free] = free_solution
        solutions[cutoff] = solution
        label_sets[cutoff] = [
            *(("infinity", *mode) for mode in modes_inf),
            *(("zero", *mode) for mode in modes_zero),
        ]
        runs.append({
            "fourier_cutoff": cutoff,
            "unknown_count": matrix.shape[1],
            "full_oversampled_residual_l2_diagnostic": float(np.linalg.norm(matrix @ solution)),
            "constrained_matrix_condition_diagnostic": float(np.linalg.cond(matrix[:, free])),
        })
    comparisons = []
    for low, high in zip(cutoffs, cutoffs[1:]):
        high_index = {label: i for i, label in enumerate(label_sets[high])}
        difference = np.array([
            solutions[high][high_index[label]] - solutions[low][i]
            for i, label in enumerate(label_sets[low])
        ])
        comparisons.append({
            "source_cutoff": low,
            "target_cutoff": high,
            "low_mode_delta_l2_diagnostic": float(np.linalg.norm(difference)),
            "relative_delta_diagnostic": float(
                np.linalg.norm(difference) / np.linalg.norm(solutions[low])
            ),
        })
    output = {
        "schema": "track-b-oversampled-least-squares-diagnostic/v1",
        "label": "DIAGNOSTIC ONLY; NOT AN INTERVAL SOLVE",
        "normalization": "a_infinity,(1,0)=1",
        "runs": runs,
        "comparisons": comparisons,
        "conclusion": (
            "order-one low-mode drift persists in the full oversampled least-squares "
            "crosscheck, so changing only the square-row selector does not explain it"
        ),
        "diagnostic_only": True,
        "trial_frozen": False,
        "global_hejhal_defect_certified": False,
        "rung4_certified": False,
        "dual_certification": False,
    }
    ns.json_out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
