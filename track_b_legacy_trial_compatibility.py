#!/usr/bin/env python3
"""Record validity and incompatibility of legacy floor/mass trial artifacts."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from track_b_two_cusp_data import canonical_hash


ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--floor", type=Path, default=ROOT / "track_b_floor_stability_d10_result.json")
    parser.add_argument("--mass", type=Path, default=ROOT / "track_b_projected_mass_arb_result.json")
    parser.add_argument("--current-trial", type=Path, default=ROOT / "track_b_observable_M5_result.json")
    parser.add_argument("--json-out", type=Path, default=ROOT / "track_b_legacy_trial_compatibility.json")
    ns = parser.parse_args()
    floor = json.loads(ns.floor.read_text(encoding="utf-8"))
    mass = json.loads(ns.mass.read_text(encoding="utf-8"))
    current = json.loads(ns.current_trial.read_text(encoding="utf-8"))
    original_trial = Path(str(mass["trial"]))
    original_hash = mass.get("trial_sha256")
    original_hash_verified = bool(
        original_trial.is_file() and original_hash == sha256(original_trial)
    )
    current_hash = current.get("coefficient_hash")
    records = [
        {
            "artifact": "floor",
            "path": str(ns.floor.resolve()),
            "sha256": sha256(ns.floor),
            "valid_for_original_trial": bool(floor.get("floor_residual_certified", False)),
            "original_trial": str(original_trial),
            "original_trial_sha256": original_hash,
            "original_trial_hash_verified_from_mass_provenance": original_hash_verified,
            "current_trial_coefficient_hash": current_hash,
            "compatible_with_current_trial": False,
            "reason": "legacy standalone floor artifact has no matching current coefficient/trial hash",
        },
        {
            "artifact": "projected_mass",
            "path": str(ns.mass.resolve()),
            "sha256": sha256(ns.mass),
            "valid_for_original_trial": bool(
                mass.get("theorem_DK_projected_mass_admissible", False)
                or mass.get("witness_certified", False)
            ),
            "original_trial": str(original_trial),
            "original_trial_sha256": original_hash,
            "original_trial_hash_verified": original_hash_verified,
            "current_trial_coefficient_hash": current_hash,
            "compatible_with_current_trial": False,
            "reason": "projected mass proposition is certified for the recorded fixed legacy trial only",
        },
    ]
    output = {
        "schema": "track-b-legacy-trial-compatibility/v1",
        "current_trial_path": str(ns.current_trial.resolve()),
        "current_trial_coefficient_hash": current_hash,
        "artifacts": records,
        "legacy_artifacts_preserved": True,
        "legacy_artifacts_imported_into_current_trial": False,
        "global_hejhal_defect_certified": False,
        "rung4_certified": False,
        "dual_certification": False,
    }
    output["compatibility_ledger_hash"] = canonical_hash(output)
    ns.json_out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "legacy_artifacts_preserved": True,
        "legacy_artifacts_imported_into_current_trial": False,
        "records": [{q["artifact"]: {
            "valid_for_original_trial": q["valid_for_original_trial"],
            "compatible_with_current_trial": q["compatible_with_current_trial"],
        }} for q in records],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
