#!/usr/bin/env python3
"""Floating diagnostic for the target's level-one Picard symmetry class.

This reads the fixed six-copy trial and checks the infinity-cusp coefficient
relations of Then's D/G/C/H classes.  It is evidence for choosing a theorem
projector, not an interval certificate and not a rung-4 input.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from six_copy_hejhal import component_mid_rad, gaussian_modes


SIGNS = {
    "D": (+1, +1),  # a_{i beta}=+a_beta, a_{conj beta}=+a_beta
    "G": (+1, -1),
    "C": (-1, +1),
    "H": (-1, -1),
}


def main() -> int:
    src = Path("six_copy_hejhal_balanced_coeffs.json")
    data = json.loads(src.read_text(encoding="utf-8"))
    row = data["best"]
    modes = gaussian_modes(int(data["parameters"]["M"]))
    coeff = np.array(
        [complex(float(q["real"]), float(q["imag"])) for q in row["coefficients"]],
        dtype=np.complex128,
    )[: len(modes)]
    index = {(a, b): k for k, (a, b, _nn) in enumerate(modes)}
    rot = np.array([coeff[index[(-b, a)]] for a, b, _nn in modes])
    conj = np.array([coeff[index[(a, -b)]] for a, b, _nn in modes])
    denom = max(float(np.linalg.norm(coeff)), 1e-300)
    classes = {}
    for name, (s_rot, s_conj) in SIGNS.items():
        dr = float(np.linalg.norm(rot - s_rot * coeff) / denom)
        dc = float(np.linalg.norm(conj - s_conj * coeff) / denom)
        classes[name] = {
            "rotation_relation_relative_l2": dr,
            "reflection_relation_relative_l2": dc,
            "combined_root_sum_square": float((dr * dr + dc * dc) ** 0.5),
        }
    best = min(classes, key=lambda name: classes[name]["combined_root_sum_square"])

    # Small horosphere diagnostic for the actual pointwise oldspace and
    # quarter-turn-odd projector.  This uses Arb K midpoints but floating phases,
    # coefficients, sums, and quotients; it is deliberately not certifying.
    seeds = ((0.21, 0.08), (0.33, -0.11), (0.14, 0.27))
    xy = []
    for x1, x2 in seeds:
        xy.extend(((x1, x2), (-x2, x1), (-x1, -x2), (x2, -x1)))
    points = np.array([(x1, x2, 0.9) for x1, x2 in xy], dtype=float)
    qpoints = points.copy()
    qpoints[:, 0] = -points[:, 1]  # z -> i z
    qpoints[:, 1] = points[:, 0]
    modes0 = gaussian_modes(5 * int(data["parameters"]["M"]))
    all_coeff = np.array(
        [complex(float(q["real"]), float(q["imag"])) for q in row["coefficients"]],
        dtype=np.complex128,
    )
    vals, qvals = [], []
    for copy in range(6):
        a, _ = component_mid_rad(
            copy, points, modes, modes0, float(row["r"]), 160
        )
        aq, _ = component_mid_rad(
            copy, qpoints, modes, modes0, float(row["r"]), 160
        )
        vals.append(a @ all_coeff)
        qvals.append(aq @ all_coeff)
    vals = np.asarray(vals)
    qvals = np.asarray(qvals)
    old = np.mean(vals, axis=0)
    old_q = np.mean(qvals, axis=0)
    target = 0.5 * (old - old_q)
    total_mass = max(float(np.sum(np.abs(vals) ** 2)), 1e-300)
    old_fraction = float(6.0 * np.sum(np.abs(old) ** 2) / total_mass)
    target_fraction = float(6.0 * np.sum(np.abs(target) ** 2) / total_mass)
    out = {
        "status": "floating coefficient diagnostic; not a certificate",
        "source": str(src.resolve()),
        "r": row["r"],
        "M": data["parameters"]["M"],
        "n_infinity_modes": len(modes),
        "class_relations": classes,
        "best_class": best,
        "published_target_class_from_Then_table": "C",
        "agrees_with_published_class": best == "C",
        "quarter_turn_operator": {
            "geometric_order_on_H3": 4,
            "square_on_base": "z -> -z = diag(i,-i).z",
            "diag_i_minus_i_is_in_PSL2_Zi": True,
            "pullback_order_on_level_one_quotient": 2,
            "I_minus_Q_over_2_is_exact_orthogonal_projector": True,
            "eigenvalue_simplicity_required": False,
        },
        "sampled_projector_mass": {
            "status": "floating Q-invariant horosphere diagnostic; not a lower bound",
            "height": 0.9,
            "oldspace_mass_fraction": old_fraction,
            "old_quarter_turn_odd_mass_fraction": target_fraction,
        },
        "level_one_lift_diagnostic": row.get("lifted_level1_structure"),
        "hard_map_changed": False,
    }
    target = Path("six_copy_target_symmetry_result.json")
    target.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    print(target.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
