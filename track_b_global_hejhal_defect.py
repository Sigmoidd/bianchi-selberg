#!/usr/bin/env python3
"""Fail-closed global Track-B Hejhal defect channel aggregator.

Only compatible, hash-bound continuum channels may enter the theorem total.
Missing channels are represented by ``None`` and prevent construction of a
numeric global bound; they are never silently replaced by collocation data.
"""
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


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def optional_channel(path: Path | None, channel: str, theorem_hash: str,
                     coefficient_hash: str) -> tuple[dict[str, Any] | None, str | None]:
    if path is None or not path.is_file():
        return None, f"missing {channel} certificate"
    data = load(path)
    if data.get("theorem_defect_definition_hash") != theorem_hash:
        return None, f"{channel} theorem-definition hash mismatch"
    if data.get("coefficient_vector_hash") != coefficient_hash:
        return None, f"{channel} coefficient-vector hash mismatch"
    if not data.get("certified", False):
        return None, f"{channel} certificate flag is false"
    return data, None


def build(ns: argparse.Namespace) -> dict[str, Any]:
    theorem = load(ns.theorem_defect_definition)
    partition = load(ns.global_partition)
    floor = load(ns.floor_certificate)
    trial = load(ns.two_cusp_result)
    tail = load(ns.tail_result)
    ladder = load(ns.cutoff_ladder)
    theorem_hash = str(theorem.get("theorem_defect_definition_hash"))
    coefficient_hash = str(trial.get("coefficient_vector_hash"))

    blockers: list[str] = []
    if not theorem.get("theorem_defect_norm_frozen", False):
        blockers.append("the exact theorem defect norm is not frozen")
    if not theorem.get("global_threshold_certified_for_requested_trial", False):
        blockers.append("no admissible threshold is certified for this coefficient vector")
    if not partition.get("global_partition_certified", False):
        blockers.append("global partition certificate is false")
    if not partition.get("global_weight_bounds_certified", False):
        blockers.append("global weight certificate is false")
    if not trial.get("physical_residual_certified", False):
        blockers.append("full finite physical residual is not certified")
    if trial.get("theorem_defect_definition_hash") != theorem_hash:
        blockers.append("two-cusp result theorem-definition hash mismatch")
    if not ladder.get("all_verified_solves_certified", False):
        blockers.append("cutoff ladder has an uncertified solve")
    if not ladder.get("cutoff_convergence_checked", False):
        blockers.append("cutoff convergence is not certified")
    if not tail.get("requested_fourier_tail_certified", False):
        blockers.append("requested unknown-expansion Fourier tail is uncertified")

    # The closing floor artifact did not preserve trial_sha256 or coefficient
    # map/hash.  Its known generating code defaults to a different fixed trial,
    # so it cannot be imported into this M=1 field.
    floor_compatible = bool(
        floor.get("trial_sha256") == sha256(ns.two_cusp_result)
        and floor.get("coefficient_vector_hash") == coefficient_hash
        and floor.get("theorem_defect_definition_hash") == theorem_hash
    )
    if not floor.get("floor_residual_certified", False):
        blockers.append("floor residual certificate is false")
    if not floor_compatible:
        blockers.append("floor certificate lacks compatible trial/coefficient/theorem provenance")

    cusp_inf, error = optional_channel(
        ns.cusp_infinity, "cusp-infinity continuum", theorem_hash, coefficient_hash
    )
    if error:
        blockers.append(error)
    cusp_zero, error = optional_channel(
        ns.cusp_zero, "cusp-zero continuum", theorem_hash, coefficient_hash
    )
    if error:
        blockers.append(error)
    faces, error = optional_channel(
        ns.nonfloor_faces, "non-floor face", theorem_hash, coefficient_hash
    )
    if error:
        blockers.append(error)
    reprojection, error = optional_channel(
        ns.reprojection, "reprojection", theorem_hash, coefficient_hash
    )
    if error:
        blockers.append(error)

    channels = [
        {
            "id": "floor",
            "kind": "theorem_commutator_L2",
            "upper": floor.get("floor_l2_upper") if floor_compatible else None,
            "certified": bool(floor.get("floor_residual_certified", False) and floor_compatible),
            "source_sha256": sha256(ns.floor_certificate),
            "aggregation_rule": "support-aware rule not yet available globally",
        },
        {
            "id": "cusp_infinity_continuum",
            "kind": "theorem_commutator_L2",
            "upper": None if cusp_inf is None else cusp_inf.get("upper"),
            "certified": cusp_inf is not None,
            "aggregation_rule": None,
        },
        {
            "id": "cusp_zero_continuum",
            "kind": "theorem_commutator_L2",
            "upper": None if cusp_zero is None else cusp_zero.get("upper"),
            "certified": cusp_zero is not None,
            "aggregation_rule": None,
        },
        {
            "id": "nonfloor_face_value",
            "kind": "d0 channel before B0/B1 combination",
            "upper": None if faces is None else faces.get("value_upper"),
            "certified": faces is not None and faces.get("value_upper") is not None,
            "aggregation_rule": None,
        },
        {
            "id": "nonfloor_face_gradient",
            "kind": "d1 channel before B0/B1 combination",
            "upper": None if faces is None else faces.get("gradient_upper"),
            "certified": faces is not None and faces.get("gradient_upper") is not None,
            "aggregation_rule": None,
        },
        {
            "id": "unknown_infinite_fourier_tail",
            "kind": "requested approximation tail",
            "upper": tail.get("unknown_infinite_expansion_tail", {}).get("value_tail_upper"),
            "certified": bool(tail.get("requested_fourier_tail_certified", False)),
            "aggregation_rule": None,
        },
        {
            "id": "reprojection",
            "kind": "coefficient-to-field discrepancy",
            "upper": None if reprojection is None else reprojection.get("upper"),
            "certified": reprojection is not None,
            "aggregation_rule": None,
        },
    ]
    channel_ids = [q["id"] for q in channels]
    exact_once = len(channel_ids) == len(set(channel_ids))
    if not exact_once:
        blockers.append("duplicate channel id")

    threshold = theorem.get("threshold", {})
    result: dict[str, Any] = {
        "schema": "track-b-global-hejhal-defect/v1",
        "label": "GLOBAL HEJHAL DEFECT CERTIFICATE",
        "arb_bits": trial.get("arb_bits"),
        "spectral_parameter_interval": trial.get("spectral_parameter_interval"),
        "fourier_cutoff": trial.get("fourier_cutoff"),
        "theorem_defect_definition_hash": theorem_hash,
        "partition_hash": partition.get("partition_definition_hash"),
        "weight_hash": partition.get("weight_formulas_hash"),
        "floor_certificate_hash": sha256(ns.floor_certificate),
        "two_cusp_assembly_hash": trial.get("assembly_hash"),
        "coefficient_vector_hash": coefficient_hash,
        "coefficient_to_field_map_hash": trial.get("coefficient_to_field_map_hash"),
        "physical_collocation_residual_upper": trial.get("physical_residual_l2_upper"),
        "physical_collocation_is_theorem_channel": False,
        "cusp_infinity_continuum_upper": None,
        "cusp_zero_continuum_upper": None,
        "nonfloor_face_value_upper": None,
        "nonfloor_face_gradient_upper": None,
        "floor_upper": floor.get("floor_l2_upper") if floor_compatible else None,
        "fourier_tail_upper": None,
        "bessel_tail_upper": None,
        "reprojection_error_upper": None,
        "global_defect_upper": None,
        "allowed_defect_interval": threshold.get("allowed_defect_interval"),
        "allowed_defect_lower": threshold.get("allowed_defect_lower"),
        "certified_margin_lower": None,
        "certified_ratio_upper": None,
        "resolved_bessel_fallback_count": trial.get("unresolved_bessel_fallback_count"),
        "tail_majorant_count": tail.get("tail_majorant_count", 0),
        "continuum_cell_count": 0,
        "face_count": 0,
        "channel_ledger": channels,
        "channel_ids_exact_once": exact_once,
        "aggregation_complete": False,
        "aggregation_note": "No numeric total is formed while a required theorem channel is absent.",
        "cutoff_convergence_checked": bool(ladder.get("cutoff_convergence_checked", False)),
        "cutoff_stability_passed": False,
        "independent_global_verification": False,
        "floor_residual_certified": bool(floor.get("floor_residual_certified", False)),
        "floor_certificate_compatible": floor_compatible,
        "global_partition_certified": bool(partition.get("global_partition_certified", False)),
        "global_weight_bounds_certified": bool(partition.get("global_weight_bounds_certified", False)),
        "two_cusp_assembly_certified": bool(trial.get("two_cusp_assembly_certified", False)),
        "full_physical_residual_certified": bool(trial.get("physical_residual_certified", False)),
        "both_cusp_continuum_certified": cusp_inf is not None and cusp_zero is not None,
        "all_nonfloor_faces_certified": faces is not None,
        "fourier_tails_certified": bool(tail.get("requested_fourier_tail_certified", False)),
        "reprojection_certified": reprojection is not None,
        "global_threshold_certified": bool(theorem.get("global_threshold_certified_for_requested_trial", False)),
        "blockers": blockers,
        "global_hejhal_defect_certified": False,
        "rung4_certified": False,
        "hejhal_existence_certified": False,
        "dual_certification": False,
        "status": "RED",
    }
    result["aggregation_ledger_hash"] = canonical_hash(channels)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--theorem-defect-definition", type=Path, default=ROOT / "track_b_theorem_defect_definition.json")
    parser.add_argument("--global-partition", type=Path, default=ROOT / "track_b_global_partition_result.json")
    parser.add_argument("--floor-certificate", type=Path, default=ROOT / "track_b_floor_stability_d10_result.json")
    parser.add_argument("--two-cusp-result", type=Path, default=ROOT / "track_b_two_cusp_result.json")
    parser.add_argument("--tail-result", type=Path, default=ROOT / "track_b_two_cusp_tail_result.json")
    parser.add_argument("--cutoff-ladder", type=Path, default=ROOT / "track_b_cutoff_ladder_result.json")
    parser.add_argument("--cusp-infinity", type=Path, default=None)
    parser.add_argument("--cusp-zero", type=Path, default=None)
    parser.add_argument("--nonfloor-faces", type=Path, default=None)
    parser.add_argument("--reprojection", type=Path, default=None)
    parser.add_argument("--json-out", type=Path, default=ROOT / "track_b_global_hejhal_defect_result.json")
    ns = parser.parse_args()
    result = build(ns)
    ns.json_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "global_hejhal_defect_certified": False,
        "rung4_certified": False,
        "blockers": result["blockers"],
    }, indent=2))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
