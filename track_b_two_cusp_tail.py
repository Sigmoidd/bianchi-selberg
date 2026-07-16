#!/usr/bin/env python3
"""Fail-closed Fourier-tail audit for a finite two-cusp Track-B trial."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from track_b_two_cusp_data import canonical_hash


ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_tail(trial_path: Path, theorem_path: Path) -> dict[str, Any]:
    trial = json.loads(trial_path.read_text(encoding="utf-8"))
    theorem = json.loads(theorem_path.read_text(encoding="utf-8"))
    cutoff = int(trial["fourier_cutoff"])
    finite_field_tail = {
        "definition": "the recorded trial is the finite sum over its mode_order",
        "outside_cutoff_coefficients": "exactly zero by definition of this finite trial",
        "value_tail_upper": "0",
        "horizontal_gradient_tail_upper": "0",
        "vertical_gradient_tail_upper": "0",
        "hyperbolic_gradient_tail_upper": "0",
        "laplacian_tail_upper": "0",
        "certified": True,
        "note": "This is not an error bound to an unknown infinite automorphic form.",
    }
    infinite_tail = {
        "definition": "sum of unknown coefficients beyond the recorded cutoff",
        "starting_shells": {"infinity": cutoff + 1, "zero": 5 * cutoff + 1},
        "coefficient_growth_assumption": None,
        "certified_coefficient_bound": None,
        "bessel_exponential_decay_bound": None,
        "summed_tail_interval": None,
        "value_tail_upper": None,
        "horizontal_gradient_tail_upper": None,
        "vertical_gradient_tail_upper": None,
        "hyperbolic_gradient_tail_upper": None,
        "certified": False,
        "blocker": (
            "No coefficients outside the independently solved finite space and no "
            "certified coefficient-growth envelope are present.  Bessel decay alone "
            "cannot bound an arbitrary coefficient sequence."
        ),
    }
    result = {
        "schema": "track-b-two-cusp-fourier-tail/v1",
        "trial_path": str(trial_path.resolve()),
        "trial_sha256": sha256(trial_path),
        "coefficient_vector_hash": trial.get("coefficient_vector_hash"),
        "coefficient_to_field_map_hash": trial.get("coefficient_to_field_map_hash"),
        "theorem_defect_definition_hash": theorem.get("theorem_defect_definition_hash"),
        "fourier_cutoff": cutoff,
        "finite_trial_formal_tail": finite_field_tail,
        "unknown_infinite_expansion_tail": infinite_tail,
        "theorem_interpretation": (
            "D(K) applies directly to a finite Whittaker quasimode and requires its "
            "cusp L2 decay; it does not require approximation to an unknown infinite expansion."
        ),
        "requested_fourier_tail_certified": False,
        "resolved_retained_mode_bessel_fallback_count": trial.get(
            "unresolved_bessel_fallback_count"
        ),
        "tail_majorant_count": 0,
        "global_hejhal_defect_certified": False,
        "rung4_certified": False,
        "status": "RED",
    }
    result["tail_audit_hash"] = canonical_hash(result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trial", type=Path, default=ROOT / "track_b_two_cusp_result.json")
    parser.add_argument("--theorem-defect-definition", type=Path, default=ROOT / "track_b_theorem_defect_definition.json")
    parser.add_argument("--json-out", type=Path, default=ROOT / "track_b_two_cusp_tail_result.json")
    ns = parser.parse_args()
    result = audit_tail(ns.trial, ns.theorem_defect_definition)
    ns.json_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "finite_trial_formal_tail_certified": result["finite_trial_formal_tail"]["certified"],
        "requested_fourier_tail_certified": False,
        "blocker": result["unknown_infinite_expansion_tail"]["blocker"],
    }, indent=2))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
