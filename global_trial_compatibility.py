#!/usr/bin/env python3
r"""Read-only audit: do the existing Track B certificates refer to the same trial?

Never modifies any artifact.  Reports, for each analytic certificate, the spectral parameter,
Fourier cutoff and trial identity it carries, and whether they are mutually compatible.
Hashes of files are reported both for the bytes on disk (LF) and for the CRLF form, because
the recorded hashes in this repository were produced on Windows (see audit).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def _load(name: str):
    path = ROOT / name
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def sha_variants(path: Path) -> dict[str, str]:
    raw = path.read_bytes()
    lf = raw.replace(b"\r\n", b"\n")
    return {
        "sha256_bytes_on_disk": hashlib.sha256(raw).hexdigest(),
        "sha256_lf": hashlib.sha256(lf).hexdigest(),
        "sha256_crlf": hashlib.sha256(lf.replace(b"\n", b"\r\n")).hexdigest(),
    }


def audit() -> dict:
    physical = _load("track_b_two_cusp_result.json") or {}
    defect = _load("track_b_global_hejhal_defect_result.json") or {}
    mass = _load("track_b_projected_mass_arb_result.json") or {}
    floor = _load("track_b_floor_stability_d10_result.json") or {}
    coeff_path = ROOT / "six_copy_hejhal_balanced_coeffs.json"

    rec = {
        "physical_two_cusp_solve": {
            "r": str(physical.get("spectral_parameter_input", "")).split(",")[0],
            "fourier_cutoff_M": physical.get("fourier_cutoff"),
            "trial_identity": defect.get("coefficient_vector_hash"),
        },
        "projected_mass": {
            "r": (mass.get("parameters") or {}).get("r"),
            "M_infinity": (mass.get("parameters") or {}).get("M_infinity"),
            "M_zero": (mass.get("parameters") or {}).get("M_zero"),
            "trial_identity": mass.get("trial_sha256"),
        },
        "floor": {
            "r": floor.get("r"),
            "trial_identity": floor.get("trial_sha256"),
            "carries_trial_hash": "trial_sha256" in floor,
            "carries_spectral_parameter": any(k in floor for k in ("r", "spectral_parameter", "spectral_parameter_interval")),
        },
    }
    coeff_hashes = sha_variants(coeff_path) if coeff_path.is_file() else {}
    mass_hash = rec["projected_mass"]["trial_identity"]
    hash_note = None
    if coeff_hashes and mass_hash:
        hash_note = {
            "recorded_mass_trial_sha_matches_LF_bytes": mass_hash == coeff_hashes["sha256_lf"],
            "recorded_mass_trial_sha_matches_CRLF_bytes": mass_hash == coeff_hashes["sha256_crlf"],
        }
    r_phys, r_mass = rec["physical_two_cusp_solve"]["r"], rec["projected_mass"]["r"]
    same_r = bool(r_phys) and bool(r_mass) and abs(float(r_phys) - float(r_mass)) < 1e-12
    same_trial = (rec["physical_two_cusp_solve"]["trial_identity"] is not None
                  and rec["physical_two_cusp_solve"]["trial_identity"] == mass_hash)
    blockers = []
    if not same_r:
        blockers.append(f"spectral parameters differ: physical r={r_phys}, mass/floor trial r={r_mass}")
    if not same_trial:
        blockers.append("physical solve and projected-mass certificate refer to different coefficient vectors")
    if not rec["floor"]["carries_trial_hash"]:
        blockers.append("floor certificate carries no trial hash")
    if not rec["floor"]["carries_spectral_parameter"]:
        blockers.append("floor certificate carries no spectral parameter")
    return {
        "schema": "global-trial-compatibility/v1",
        "tier": "AUDIT (read-only)",
        "artifacts": rec,
        "balanced_coeffs_file_hashes": coeff_hashes,
        "hash_note": hash_note,
        "same_spectral_parameter": same_r,
        "same_trial": same_trial,
        "compatible": not blockers,
        "blockers": blockers,
        "defect_flags": {k: defect.get(k) for k in (
            "global_hejhal_defect_certified", "rung4_certified", "floor_certificate_compatible",
            "both_cusp_continuum_certified", "all_nonfloor_faces_certified", "global_threshold_certified")},
    }


def parity_of_balanced_trial() -> dict | None:
    """EXPERIMENTAL: coefficient-norm Q-parity of the mass/floor trial (not an L2 mass fraction)."""
    path = ROOT / "six_copy_hejhal_balanced_coeffs.json"
    if not path.is_file():
        return None
    import global_quarter_turn as gq
    best = json.loads(path.read_text(encoding="utf-8"))["best"]
    c = np.array([complex(float(x["real"]), float(x["imag"])) for x in best["coefficients"]])
    mi, mz = gq.gaussian_modes(8), gq.gaussian_modes(40)
    n = len(c)
    assert n == len(mi) + len(mz)
    perm = np.zeros((n, n))
    for off, modes in ((0, mi), (len(mi), mz)):
        idx = {(m[0], m[1]): k for k, m in enumerate(modes)}
        for k, m in enumerate(modes):
            perm[off + k, off + idx[gq.rot_i((m[0], m[1]))]] = 1.0
    even, odd = 0.5 * (c + perm @ c), 0.5 * (c - perm @ c)
    nc = np.linalg.norm(c)
    return {
        "tier": "EXPERIMENTAL (coefficient-norm fractions; not L2 mass fractions)",
        "r": best["r"], "M_infinity": 8, "M_zero": 40,
        "even_fraction": float(np.linalg.norm(even) / nc),
        "odd_fraction": float(np.linalg.norm(odd) / nc),
    }


def main() -> int:
    out = audit()
    out["parity_of_balanced_trial"] = parity_of_balanced_trial()
    Path("global_trial_compatibility_result.json").write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("compatible", "blockers", "hash_note", "parity_of_balanced_trial")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
