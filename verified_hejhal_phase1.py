#!/usr/bin/env python3
"""Phase-1 repair and scaling audit for the dual-certification Hejhal path.

This is deliberately separate from the certified theorem residual.  It:

* cross-validates Arb K_{ir} against mpmath;
* demonstrates r-dependence of the repaired model matrix;
* records every diagonal in Aeq = L A D^{-1} C^{-1} R;
* maps the SVD vector back to physical coefficients exactly;
* compares the old and repaired model residuals.

It does not flip any certification flag and does not promote sampled face
residuals to Theorem D(K) inputs.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np

from verified_kir import kir_enclosure


DEFAULT_PIPELINE = Path(r"C:\Users\Admin\.grok\worktrees\bearings-bianchi-selberg\np9\dual_certification_")


def equilibrate_with_ledger(A: np.ndarray, n_iter: int = 6):
    """Return (L A R, Ldiag, Rdiag), retaining every diagonal factor."""
    B = np.array(A, dtype=np.complex128, copy=True)
    left = np.ones(B.shape[0], dtype=float)
    right = np.ones(B.shape[1], dtype=float)
    for _ in range(n_iter):
        row = np.maximum(np.max(np.abs(B), axis=1), 1e-300)
        rs = 1.0 / row
        B = rs[:, None] * B
        left *= rs
        col = np.maximum(np.max(np.abs(B), axis=0), 1e-300)
        cs = 1.0 / col
        B = B * cs[None, :]
        right *= cs
    reconstruction = np.linalg.norm(B - left[:, None] * A * right[None, :])
    reconstruction /= max(np.linalg.norm(B), 1e-300)
    return B, left, right, float(reconstruction)


def svd_record(A: np.ndarray) -> dict[str, Any]:
    s = np.linalg.svd(A, compute_uv=False)
    cutoff = 1e-12 * max(float(s[0]), 1e-300)
    rank = int(np.count_nonzero(s > cutoff))
    return {
        "shape": list(A.shape), "sigma_max": float(s[0]), "sigma_min": float(s[-1]),
        "condition_2": float(s[0] / max(s[-1], 1e-300)),
        "rank_rtol_1e-12": rank,
        "nullity_rtol_1e-12": int(A.shape[1] - rank),
    }


def patch_model_backend(tc, pairing, hc, bits: int):
    """Patch midpoint calls only; all values still originate in Arb balls."""
    def abs_k(x: float, r: float) -> float:
        q = kir_enclosure(r, x, bits)
        return abs(q.midpoint)

    def signed_k(x: float, r: float) -> float:
        return kir_enclosure(r, x, bits).midpoint

    def k_rad(x: float, r: float) -> float:
        q = kir_enclosure(r, x, bits)
        return q.radius

    # The legacy model S assumes nonnegative amplitudes.  Preserve that model's
    # algebra for the r-dependence experiment, but use the signed value in the
    # independently evaluated Fourier field.
    tc.k_amp = abs_k
    tc.k_amp_radius = k_rad
    hc.k_bessel_amp = signed_k
    pairing.k_bessel_amp = signed_k


def build_and_solve(tc, M: int, r: float, Y0: float):
    system = tc.build_block_system(M, r, Y0, theta=0.5)
    physical, radius = tc.coupled_mid_rad(system)
    pre, pre_rad, D = tc.precondition_right(
        physical, radius, system.w_inf, system.w_0
    )
    col = np.maximum(np.linalg.norm(pre, axis=0), 1e-300)
    G = pre / col[None, :]
    Aeq, left, right, ledger_error = equilibrate_with_ledger(G, 6)
    _, s, vh = np.linalg.svd(Aeq, full_matrices=False)
    veq = vh[-1].conj()

    # Aeq = L * physical * D^{-1} * C^{-1} * R.
    # Hence physical coefficients are D^{-1} C^{-1} R veq.
    coeff = (right * veq) / (D * col)
    coeff /= max(np.linalg.norm(coeff), 1e-300)
    weighted_residual = np.linalg.norm(Aeq @ veq) / max(np.linalg.norm(veq), 1e-300)
    physical_rel = np.linalg.norm(physical @ coeff)
    physical_rel /= max(np.linalg.norm(physical, 2) * np.linalg.norm(coeff), 1e-300)
    interval_action_bound = np.linalg.norm(physical @ coeff) + np.linalg.norm(
        radius @ np.abs(coeff)
    )
    return {
        "system": system, "physical": physical, "radius": radius, "coeff": coeff,
        "record": {
            "physical": svd_record(physical), "preconditioned": svd_record(pre),
            "column_normalized": svd_record(G), "equilibrated": svd_record(Aeq),
            "equilibration_ledger_relative_error": ledger_error,
            "weighted_svd_residual": float(weighted_residual),
            "physical_relative_residual": float(physical_rel),
            "interval_action_2norm_bound": float(interval_action_bound),
            "left_scale_range": [float(left.min()), float(left.max())],
            "right_scale_range": [float(right.min()), float(right.max())],
            "amplitude_scale_range": [float(D.min()), float(D.max())],
            "column_scale_range": [float(col.min()), float(col.max())],
            "sigma_values_tail": [float(x) for x in s[-5:]],
        },
    }


def mpmath_crosscheck(points: list[tuple[float, float]], bits: int):
    import mpmath as mp
    mp.mp.dps = max(60, int(bits * 0.30103) + 10)
    rows = []
    for r, x in points:
        q = kir_enclosure(r, x, bits)
        ref = float(mp.re(mp.besselk(1j * mp.mpf(repr(r)), mp.mpf(repr(x)))))
        rows.append({
            **q.to_dict(), "mpmath": ref,
            "mpmath_inside_arb_float_endpoints": bool(q.lower <= ref <= q.upper),
            "midpoint_minus_mpmath": q.midpoint - ref,
        })
    return rows


def pairing_measure(pairing, coeff: np.ndarray, M: int, r: float, Y0: float, n_face: int):
    modes = pairing.gaussian_modes(M)
    per, delta, scale = pairing.measure_jumps(
        coeff, modes, r, Y0, 0.5, n_face, 1.0 / math.sqrt(2.0), True, False
    )
    return {"delta_sampled": delta, "scale": scale, "per_pairing": per}


def run(pipeline: Path, M: int, r: float, Y0: float, bits: int, n_face: int):
    sys.path.insert(0, str(pipeline))
    import two_cusp_hejhal_N5 as tc
    import delta_aut_pairing as pairing
    import hejhal_conditioning as hc

    patch_model_backend(tc, pairing, hc, bits)
    base = build_and_solve(tc, M, r, Y0)
    lo = build_and_solve(tc, M, 6.0, Y0)
    hi = build_and_solve(tc, M, 8.0, Y0)
    r_change = np.linalg.norm(lo["physical"] - hi["physical"])
    r_change /= max(np.linalg.norm(lo["physical"]), 1e-300)
    n = base["system"].n
    coeff_inf = base["coeff"][:n]
    coeff_0 = base["coeff"][n:]
    return {
        "status": "diagnostic only; hard map unchanged",
        "parameters": {"M": M, "r": r, "Y0": Y0, "bits": bits, "n_face": n_face},
        "backend": {
            "name": "python-flint acb.bessel_k", "fail_closed": True,
            "mpmath_crosscheck": mpmath_crosscheck([
                (r, 2 * math.pi * Y0),
                (r, 2 * math.pi * math.sqrt(2) * Y0),
                (r, 2 * math.pi * math.sqrt(M) * Y0),
            ], bits),
            "repaired_model_relative_change_r6_to_r8": float(r_change),
        },
        "solve": base["record"],
        "sampled_face_residual": {
            "infinity": pairing_measure(pairing, coeff_inf, M, r, Y0, n_face),
            "zero": pairing_measure(pairing, coeff_0, M, r, Y0, n_face),
        },
        "limitations": [
            "The two-cusp scattering block remains the legacy model S.",
            "The face residual is sampled, not a continuum interval norm.",
            "The SVD is floating point although every K midpoint came from an Arb enclosure.",
            "No value here is admissible as eta in Theorem D(K).",
        ],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline", type=Path, default=DEFAULT_PIPELINE)
    ap.add_argument("--M", type=int, default=28)
    ap.add_argument("--r", type=float, default=6.7439020359331625)
    ap.add_argument("--Y0", type=float, default=0.8)
    ap.add_argument("--bits", type=int, default=160)
    ap.add_argument("--n-face", type=int, default=16)
    ap.add_argument("--json-out", type=Path, default=Path("verified_hejhal_phase1_result.json"))
    ns = ap.parse_args()
    out = run(ns.pipeline, ns.M, ns.r, ns.Y0, ns.bits, ns.n_face)
    ns.json_out.write_text(json.dumps(out, indent=2, allow_nan=False), encoding="utf-8")
    print(ns.json_out.resolve())
    print(json.dumps({
        "r6_to_r8_change": out["backend"]["repaired_model_relative_change_r6_to_r8"],
        "physical_residual": out["solve"]["physical_relative_residual"],
        "delta_inf": out["sampled_face_residual"]["infinity"]["delta_sampled"],
        "delta_0": out["sampled_face_residual"]["zero"]["delta_sampled"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
