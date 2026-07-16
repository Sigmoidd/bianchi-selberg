#!/usr/bin/env python3
"""High-precision diagnostic decomposition of the old six-point M=4 nullspace."""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from flint import acb, acb_mat, arb, ctx
import numpy as np

from track_b_two_cusp_data import canonical_hash, cusp_modes
from track_b_two_cusp_hejhal import (
    ValidatedWhittaker, assemble_physical_rows, midpoint_matrix,
)


ROOT = Path(__file__).resolve().parent


def midpoint_complex(value: acb) -> complex:
    return complex(float(value.real.mid()), float(value.imag.mid()))


def unit_orbit(a: int, b: int) -> tuple[tuple[int, int], ...]:
    return tuple(sorted({(a, b), (-b, a), (-a, -b), (b, -a)}))


def energies(labels: list[tuple[str, int, int, int]], vector: np.ndarray) -> dict[str, Any]:
    by_cusp: dict[str, float] = defaultdict(float)
    by_shell: dict[str, float] = defaultdict(float)
    by_orbit: dict[str, float] = defaultdict(float)
    modes = []
    for label, value in zip(labels, vector):
        cusp, a, b, norm = label
        energy = float(abs(value) ** 2)
        by_cusp[cusp] += energy
        shell = f"{cusp}:{norm if cusp == 'infinity' else str(norm) + '/5'}"
        by_shell[shell] += energy
        by_orbit[str((cusp, unit_orbit(a, b)))] += energy
        modes.append({
            "cusp": cusp, "mode": [a, b, norm], "energy": energy,
            "real_energy": float(value.real * value.real),
            "imaginary_energy": float(value.imag * value.imag),
            "symmetry_orbit": [list(q) for q in unit_orbit(a, b)],
            "six_copy_fiber": "shared cusp-channel coefficient",
        })
    modes.sort(key=lambda q: q["energy"], reverse=True)
    shells = sorted(by_shell.items(), key=lambda q: q[1], reverse=True)
    return {
        "by_cusp": dict(sorted(by_cusp.items())),
        "by_radial_shell": [{"shell": k, "energy": v} for k, v in shells],
        "dominant_shell": None if not shells else {"shell": shells[0][0], "energy": shells[0][1]},
        "by_symmetry_orbit": [
            {"orbit": k, "energy": v}
            for k, v in sorted(by_orbit.items(), key=lambda q: q[1], reverse=True)
        ],
        "dominant_modes": modes[:16],
        "real_component_energy": float(sum(value.real * value.real for value in vector)),
        "imaginary_component_energy": float(sum(value.imag * value.imag for value in vector)),
        "fiber_note": "the physical formulation has cusp-channel coefficients shared by all six fibers",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", type=int, default=256)
    parser.add_argument("--r", default="6.7439020359331625")
    parser.add_argument("--json-out", type=Path, default=ROOT / "track_b_M4_nullspace_audit.json")
    ns = parser.parse_args()
    ctx.prec = max(192, ns.bits)
    rows, ledger, modes_inf, modes_zero = assemble_physical_rows(
        4, 6, ValidatedWhittaker(arb(ns.r)), None
    )
    n = len(rows[0])
    norm_index = next(i for i, mode in enumerate(modes_inf) if mode[:2] == (1, 0))
    normal = [acb(0) for _ in range(n)]
    normal[norm_index] = acb(1)
    interval_matrix = acb_mat([*rows, normal])
    gram = interval_matrix.conjugate().transpose() * interval_matrix
    eigenvalues, right = gram.eig(right=True, algorithm="approx")
    order = sorted(range(n), key=lambda i: abs(midpoint_complex(eigenvalues[i])))
    physical_mid = midpoint_matrix(rows)
    labels = [
        *(("infinity", *mode) for mode in modes_inf),
        *(("zero", *mode) for mode in modes_zero),
    ]
    low3_inf, low3_zero = cusp_modes(3)
    old_labels = {
        *(("infinity", *mode) for mode in low3_inf),
        *(("zero", *mode) for mode in low3_zero),
    }
    null_vectors = []
    for vector_number, column in enumerate(order[:7]):
        vector = np.array([midpoint_complex(right[i, column]) for i in range(n)])
        vector /= np.linalg.norm(vector)
        action = physical_mid @ vector
        grouped: dict[str, list[float]] = defaultdict(list)
        for value, record in zip(action, ledger):
            grouped[f"block:{record['matrix_block']}"] .append(float(abs(value)))
            grouped[f"relation:{record['group_element']}"] .append(float(abs(value)))
            grouped[f"fiber:{record['source_copy']}"] .append(float(abs(value)))
        sensitivity = [
            {
                "row_family": key,
                "rms_action": float(np.linalg.norm(values) / max(len(values), 1) ** 0.5),
                "maximum_action": max(values),
            }
            for key, values in grouped.items()
        ]
        sensitivity.sort(key=lambda q: (q["rms_action"], q["row_family"]))
        new_energy = float(sum(
            abs(value) ** 2 for label, value in zip(labels, vector) if label not in old_labels
        ))
        decomposition = energies(labels, vector)
        null_vectors.append({
            "basis_vector": vector_number,
            "gram_eigenvalue_high_precision_diagnostic": str(eigenvalues[column]),
            "physical_action_l2": float(np.linalg.norm(action)),
            "fraction_in_new_M4_modes": new_energy,
            "fraction_in_old_M_le_3_modes": 1.0 - new_energy,
            "row_families_least_sensitive": sensitivity[:12],
            **decomposition,
        })
    singular = np.linalg.svd(np.vstack([
        physical_mid,
        np.eye(1, n, norm_index, dtype=np.complex128),
    ]), compute_uv=False)
    result = {
        "schema": "track-b-M4-six-point-nullspace-audit/v1",
        "label": "DIAGNOSTIC ONLY",
        "arb_midpoint_precision_bits": int(ctx.prec),
        "matrix_shape": [len(rows) + 1, n],
        "observed_numerical_nullity": 7,
        "numpy_smallest_singular_values": [float(q) for q in singular[-12:]],
        "basis_method": "256-bit Arb midpoint Gram matrix; approximate high-precision acb eigensystem",
        "basis_nonuniqueness_note": "the seven nearly degenerate null directions may rotate within their invariant subspace",
        "null_vectors": null_vectors,
        "diagnostic_only": True,
        "cutoff_rank_certified": False,
        "trial_frozen": False,
        "global_hejhal_defect_certified": False,
        "rung4_certified": False,
        "dual_certification": False,
    }
    result["audit_hash"] = canonical_hash(result)
    ns.json_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "observed_numerical_nullity": 7,
        "nullspace_basis_vectors_reported": len(null_vectors),
        "new_mode_energy_range": [
            min(q["fraction_in_new_M4_modes"] for q in null_vectors),
            max(q["fraction_in_new_M4_modes"] for q in null_vectors),
        ],
        "diagnostic_only": True,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
