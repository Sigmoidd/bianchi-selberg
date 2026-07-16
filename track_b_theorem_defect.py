#!/usr/bin/env python3
"""Freeze the exact Track-B D(K) defect interface and its trial dependencies.

This module deliberately does not evaluate a residual.  It records the norm
and threshold-producing formula used by ``theorem_DK_sixcopy.tex`` and checks
whether the available projected-mass witness belongs to the requested trial.
An incompatible mass witness leaves the admissible threshold unset.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from flint import arb, ctx


ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_hash(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _required_theorem_clauses(text: str) -> dict[str, bool]:
    return {
        "working_draft_disclosed": "Working draft" in text,
        "patch_residual_defined": r"R:=\norm{(\D-\lambda)U}_2" in text,
        "direct_weighted_bound_defined": r"\norm{B_0d_0+B_1d_1}_2" in text,
        "supremum_fallback_defined": r"R\le\tau+b_0\delta_0+b_1\delta_1" in text,
        "finite_whittaker_tau_zero": r"\tau=0" in text and "termwise" in text,
        "target_mass_formula_defined": r"\mu_{\rm target}=m_B-\varepsilon_B>0" in text,
        "spectral_error_formula_defined": (
            r"\frac{\tau+b_0\delta_0+b_1\delta_1}{\mu_{\rm target}}" in text
        ),
        "two_cusp_decay_required": "correctly scaled coordinates of both subgroup cusps" in text,
    }


def build_definition(
    theorem_path: Path,
    mass_path: Path,
    trial_path: Path | None,
    full_width_tolerance: str,
    bits: int,
) -> dict[str, Any]:
    ctx.prec = max(128, bits)
    theorem_text = theorem_path.read_text(encoding="utf-8")
    clauses = _required_theorem_clauses(theorem_text)
    norm_frozen = all(clauses.values())

    mass = json.loads(mass_path.read_text(encoding="utf-8"))
    mass_trial_path = Path(str(mass.get("trial", "")))
    mass_parameters = mass.get("parameters", {})
    mu_text = str(
        mass.get("plateau_construction", {}).get("certified_mu_B_lower", "nan")
    )

    trial: dict[str, Any] | None = None
    trial_hash = None
    trial_r = None
    trial_cutoff = None
    if trial_path is not None:
        trial = json.loads(trial_path.read_text(encoding="utf-8"))
        trial_hash = sha256(trial_path)
        trial_r = str(trial.get("spectral_parameter_input", ""))
        trial_cutoff = trial.get("fourier_cutoff")

    mass_trial_sha = mass.get("trial_sha256")
    mass_trial_file_matches = bool(
        mass_trial_path.is_file()
        and mass_trial_sha
        and sha256(mass_trial_path) == mass_trial_sha
    )
    mass_r = str(mass_parameters.get("r", ""))
    requested_r_is_point = bool(trial_r and "," in trial_r and len(set(trial_r.split(","))) == 1)
    requested_r_value = trial_r.split(",")[0] if requested_r_is_point else trial_r
    same_r = trial is not None and requested_r_value == mass_r
    same_coefficient_artifact = bool(
        trial is not None
        and mass_trial_sha
        and trial_hash == mass_trial_sha
    )
    # The plateau proposition is for the fixed recorded projected trial.  An
    # equal r alone cannot transfer its mass lower bound to another vector.
    mass_compatible = bool(
        norm_frozen and mass_trial_file_matches and same_r and same_coefficient_artifact
    )

    threshold: dict[str, Any] = {
        "formula": "R_allowed = mu_target_lower * full_spectral_width_tolerance / 2",
        "justification": "|lambda_j-lambda| <= R/mu_target",
        "full_spectral_width_tolerance": full_width_tolerance,
        "mass_witness_lower": mu_text,
        "mass_witness_trial": str(mass_trial_path),
        "mass_witness_trial_sha256": mass_trial_sha,
        "mass_witness_r": mass_r,
        "mass_compatible_with_requested_trial": mass_compatible,
        "allowed_defect_interval": None,
        "allowed_defect_lower": None,
        "certified_for_requested_trial": False,
    }
    if mass_compatible:
        mu = arb(mu_text).lower()
        tolerance = arb(full_width_tolerance)
        allowed = (mu * tolerance / 2).lower()
        threshold.update({
            "allowed_defect_interval": str(allowed),
            "allowed_defect_lower": str(allowed),
            "certified_for_requested_trial": bool(allowed > 0),
        })

    definition: dict[str, Any] = {
        "schema": "track-b-theorem-defect-definition/v1",
        "label": "TRACK-B D(K) THEOREM DEFECT DEFINITION",
        "source": {
            "path": str(theorem_path.resolve()),
            "sha256": sha256(theorem_path),
            "status": "working draft",
            "references": ["lem:patch", "cor:target", "prop:trackBwitness", "cor:trackBplateau"],
        },
        "defect": {
            "global_field": "U=sum_j psi_j F_j",
            "operator": "Delta-lambda",
            "spectral_parameter_convention": "lambda=1+r^2",
            "norm_type": "L2 graph residual of the covariant Laplacian",
            "norm": "R=||(Delta-lambda)U||_{L2(X;C^6)}",
            "base_space": "X=PSL2(Z[i])\\H^3 with the induced flat C^6 bundle",
            "integration_measure": "hyperbolic volume dx1 dx2 dy / y^3",
            "pointwise_value_norm": "Euclidean C^6 fiber norm",
            "pointwise_gradient_norm": "hyperbolic cotangent norm tensor Euclidean C^6 norm",
            "value_jump": "d0=max_active ||t_kj F_j-F_k||_{C^6}",
            "gradient_jump": "d1=max_active ||nabla(t_kj F_j)-nabla F_k||_{T*X tensor C^6}",
            "value_jump_weight": "B0=sum_j |Delta psi_j|",
            "gradient_jump_weight": "B1=2 sum_j |nabla psi_j|",
            "partition_commutator": "||B0*d0+B1*d1||_L2, with the pointwise sum formed before the norm",
            "interior_term": "tau=||sum_j psi_j(Delta-lambda)F_j||_L2",
            "finite_whittaker_interior_term": "tau=0 termwise",
            "direct_bound": "R<=tau+||B0*d0+B1*d1||_L2",
            "uniform_fallback": "R<=tau+b0*delta0+b1*delta1; b_i=||B_i||_L2",
            "cusp_requirement": "zero constant term and exponential L2 decay in correctly scaled coordinates at both cusps",
            "trace_or_boundary_norm": None,
        },
        "required_clause_checks": clauses,
        "theorem_defect_norm_frozen": norm_frozen,
        "theorem_source_is_final": False,
        "requested_trial": {
            "path": None if trial_path is None else str(trial_path.resolve()),
            "sha256": trial_hash,
            "spectral_parameter_input": trial_r,
            "fourier_cutoff": trial_cutoff,
        },
        "mass_compatibility": {
            "mass_artifact_path": str(mass_path.resolve()),
            "mass_artifact_sha256": sha256(mass_path),
            "recorded_trial_file_hash_matches": mass_trial_file_matches,
            "spectral_parameter_matches": same_r,
            "coefficient_artifact_matches": same_coefficient_artifact,
            "compatible": mass_compatible,
            "blocker": None if mass_compatible else (
                "The certified plateau mass belongs to a different fixed coefficient artifact; "
                "D(K) does not transfer that lower bound by spectral-parameter proximity."
            ),
        },
        "threshold": threshold,
        "legacy_eta0_applicable": False,
        "local_floor_budget_is_global_threshold_only_for_same_fixed_trial": True,
        "global_threshold_certified_for_requested_trial": bool(
            threshold["certified_for_requested_trial"]
        ),
    }
    # The interface hash must be independent of a particular residual result:
    # residual artifacts carry this hash, so including their file hash here
    # would create a circular and non-reproducible dependency.
    hash_payload = {
        "schema": definition["schema"],
        "source": definition["source"],
        "defect": definition["defect"],
        "required_clause_checks": definition["required_clause_checks"],
        "theorem_defect_norm_frozen": definition["theorem_defect_norm_frozen"],
        "theorem_source_is_final": definition["theorem_source_is_final"],
        "threshold_definition": {
            "formula": threshold["formula"],
            "justification": threshold["justification"],
            "full_spectral_width_tolerance": threshold["full_spectral_width_tolerance"],
        },
        "legacy_eta0_applicable": definition["legacy_eta0_applicable"],
    }
    definition["theorem_defect_hash_payload"] = hash_payload
    definition["theorem_defect_definition_hash"] = canonical_hash(hash_payload)
    return definition


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--theorem", type=Path, default=ROOT / "theorem_DK_sixcopy.tex")
    parser.add_argument("--mass", type=Path, default=ROOT / "track_b_projected_mass_arb_result.json")
    parser.add_argument("--trial", type=Path, default=ROOT / "track_b_two_cusp_result.json")
    parser.add_argument("--full-width-tolerance", default="0.1")
    parser.add_argument("--bits", type=int, default=192)
    parser.add_argument("--json-out", type=Path, default=ROOT / "track_b_theorem_defect_definition.json")
    ns = parser.parse_args()
    result = build_definition(ns.theorem, ns.mass, ns.trial, ns.full_width_tolerance, ns.bits)
    ns.json_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "theorem_defect_norm_frozen": result["theorem_defect_norm_frozen"],
        "theorem_defect_definition_hash": result["theorem_defect_definition_hash"],
        "mass_compatible": result["mass_compatibility"]["compatible"],
        "global_threshold_certified_for_requested_trial": result["global_threshold_certified_for_requested_trial"],
    }, indent=2))
    return 0 if result["theorem_defect_norm_frozen"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
