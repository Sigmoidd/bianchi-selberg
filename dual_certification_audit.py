#!/usr/bin/env python3
"""Diagnostic-only audit for the external dual_certification_ worktree.

It does not alter the pipeline or claim a certificate.  The script makes the
normalisations in the model Hejhal path explicit and writes a reproducible JSON
record suitable for the accompanying mathematical audit.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


DEFAULT_PIPELINE = Path(r"C:\Users\Admin\.grok\worktrees\bearings-bianchi-selberg\np9\dual_certification_")


def finite(x: Any) -> Any:
    """Convert NumPy scalar values to JSON and preserve nonfinite values visibly."""
    x = float(x)
    return x if math.isfinite(x) else str(x)


def svd_metrics(A: np.ndarray, rtol: float = 1e-12) -> dict[str, Any]:
    s = np.linalg.svd(A, compute_uv=False)
    smax = float(s[0]) if s.size else 0.0
    smin = float(s[-1]) if s.size else 0.0
    rank = int(np.count_nonzero(s > rtol * max(smax, 1e-300)))
    return {
        "shape": list(A.shape), "sigma_max": finite(smax), "sigma_min": finite(smin),
        "condition_2": finite(smax / max(smin, 1e-300)), "rank_rtol_1e-12": rank,
        "nullity_rtol_1e-12": int(A.shape[1] - rank),
    }


def face_mode_decomposition(name, coeffs, modes, r, Y0, theta, n_face, pairing):
    """Per-mode jump amplitudes.  They are diagnostics, not an additive norm split."""
    pts = pairing.face_samples(name, n_face, Y0, y_min=1.0 / math.sqrt(2.0))
    mat = pairing.GEN[name]
    pulled = np.zeros_like(pts)
    for j, (x1, x2, y) in enumerate(pts):
        z, yp = pairing.h3_act(mat, complex(x1, x2), float(y))
        pulled[j] = (z.real, z.imag, max(yp, 1e-12))
    phi0 = pairing._mode_row(pts, modes, r, theta)
    phi1 = pairing._mode_row(pulled, modes, r, theta)
    term = (phi0 - phi1) * coeffs[None, :]
    per_mode = []
    by_norm: dict[int, np.ndarray] = {}
    for k, (a, b, nn) in enumerate(modes):
        v = term[:, k]
        per_mode.append({"index": k, "a": a, "b": b, "norm": nn,
                         "max_abs": finite(np.max(np.abs(v))),
                         "rms_abs": finite(np.sqrt(np.mean(np.abs(v) ** 2)))})
        by_norm.setdefault(nn, np.zeros(pts.shape[0], dtype=np.complex128))
        by_norm[nn] += v
    per_mode.sort(key=lambda q: float(q["max_abs"]), reverse=True)
    grouped = []
    for nn, v in by_norm.items():
        grouped.append({"norm": nn, "n_basis": sum(m[2] == nn for m in modes),
                        "max_abs_coherent": finite(np.max(np.abs(v))),
                        "rms_abs_coherent": finite(np.sqrt(np.mean(np.abs(v) ** 2)))})
    grouped.sort(key=lambda q: float(q["max_abs_coherent"]), reverse=True)
    return {"basis_sorted": per_mode, "fourier_norm_sorted": grouped}


def bessel_crosscheck(r: float) -> dict[str, Any]:
    """Independent mpmath comparison at representative Bessel arguments."""
    try:
        import mpmath as mp
        mp.mp.dps = 60
    except Exception as exc:  # pragma: no cover - environment dependent
        return {"available": False, "reason": repr(exc)}
    import hejhal_conditioning as hc
    rows = []
    for x in (2 * math.pi * 0.8, 2 * math.pi * math.sqrt(2) * 0.8,
              2 * math.pi * math.sqrt(28) * 0.8):
        truth = abs(complex(mp.besselk(1j * r, x)))
        impl = hc.k_bessel_amp(x, r)
        rows.append({"x": x, "mpmath_abs_Kir": truth, "implementation": impl,
                     "relative_difference": abs(impl - truth) / max(truth, 1e-300)})
    return {"available": True, "rows": rows}


def audit(M: int, r: float, Y0: float, n_face: int) -> dict[str, Any]:
    import two_cusp_hejhal_N5 as tc
    import production_hejhal_residual as prod
    import hejhal_iterate as itr
    import delta_aut_pairing as pairing
    from hejhal_conditioning import equilibrate

    theta = 0.5
    system = tc.build_block_system(M, r, Y0, theta=theta)
    mid, rad = tc.coupled_mid_rad(system)
    mid_p, rad_p, weights = tc.precondition_right(mid, rad, system.w_inf, system.w_0)
    col = np.maximum(np.linalg.norm(mid_p, axis=0), 1e-300)
    G = mid_p / col[None, :]
    Geq = equilibrate(G, n_iter=6)

    # Correct mapping for the *unequilibrated* G right singular vector.
    _, _, vh_g = np.linalg.svd(G, full_matrices=False)
    ug = vh_g[-1].conj()
    physical_g = ug / (weights * col)
    physical_g /= np.linalg.norm(physical_g)
    # The currently reported vector uses an SVD after equilibration but loses its
    # diagonal scalings, then omits D^{-1}; reproduce it exactly.
    a_inf_current, a0_current, meta = prod.two_cusp_near_kernel(M, r, Y0, theta=theta)
    current = np.concatenate([a_inf_current, a0_current])

    n = system.n
    modes = tc.gaussian_modes(M)
    current_inf = current[:n]
    reprojected = itr.reproject_periodized(current_inf, modes, r, Y0)
    before_plain = pairing.measure_jumps(current_inf, modes, r, Y0, theta, n_face,
                                         1.0 / math.sqrt(2.0), True, False)
    before_period = pairing.measure_jumps(current_inf, modes, r, Y0, theta, n_face,
                                          1.0 / math.sqrt(2.0), True, True)
    after_plain = pairing.measure_jumps(reprojected, modes, r, Y0, theta, n_face,
                                        1.0 / math.sqrt(2.0), True, False)
    after_period = pairing.measure_jumps(reprojected, modes, r, Y0, theta, n_face,
                                         1.0 / math.sqrt(2.0), True, True)

    # Independent consistency check of the two equivalent S actions used here.
    action_errors = []
    for z, y in ((0.13 + 0.21j, 0.8), (-0.31 + 0.17j, 0.71)):
        z1, y1 = pairing.h3_act(pairing.GEN["S"], z, y)
        z2, y2 = tc.sigma0_plane(z, y)
        action_errors.append({"z_error": abs(z1 - z2), "y_error": abs(y1 - y2)})

    def measured(rec):
        per_gen, delta, scale = rec
        return {"delta_aut": delta, "scale": scale, "per_pairing": per_gen}

    decompositions = {name: face_mode_decomposition(name, current_inf, modes, r, Y0,
                                                      theta, n_face, pairing)
                      for name in pairing.GEN}
    jump_matrix, _, _ = pairing.build_multi_pairing_V(M, r, Y0, theta=theta,
                                                        n_face=n_face)
    # Actual unscaled residuals expose the effect of reporting only conditioned SVD data.
    def rel_res(A, x):
        return float(np.linalg.norm(A @ x) / max(np.linalg.norm(A, ord=2) * np.linalg.norm(x), 1e-300))
    report = {
        "kind": "non-certifying numerical diagnosis; no hard-map flags changed",
        "parameters": {"M": M, "r": r, "Y0": Y0, "n_face": n_face},
        "assembly_model": system.meta,
        "conditioning": {
            "mid": svd_metrics(mid), "mid_preconditioned": svd_metrics(mid_p),
            "column_normalized": svd_metrics(G), "equilibrated": svd_metrics(Geq),
            "jump_operator": svd_metrics(jump_matrix),
            "amplitude_dynamic_range": finite(np.max(weights) / np.min(weights)),
            "interval_radius_to_mid_fro": finite(np.linalg.norm(rad, "fro") / max(np.linalg.norm(mid, "fro"), 1e-300)),
        },
        "coefficient_map": {
            "reported_meta": meta,
            "reported_physical_mid_relative_residual": rel_res(mid, current),
            "reported_G_relative_residual": rel_res(G, current * col),
            "correct_G_backmap_mid_relative_residual": rel_res(mid, physical_g),
            "reported_vs_correct_G_coefficient_overlap": finite(abs(np.vdot(current, physical_g))),
            "note": "The reported path does not retain equilibrate's diagonal factors and applies a_full=u/col, whereas G=mid*D^{-1}*diag(1/col).",
        },
        "residual_by_cusp": {
            "infinity": measured(before_period),
            "zero": measured(pairing.measure_jumps(current[n:], modes, r, Y0, theta, n_face,
                                                      1.0 / math.sqrt(2.0), True, True)),
        },
        "periodization_and_reprojection": {
            "before_reprojection_plain": measured(before_plain),
            "before_reprojection_periodized": measured(before_period),
            "after_reprojection_plain": measured(after_plain),
            "after_reprojection_periodized": measured(after_period),
        },
        "residual_by_pairing_basis_and_mode": decompositions,
        "independent_checks": {"S_action_agreement": action_errors,
                                 "bessel": bessel_crosscheck(r)},
    }
    return report


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline", type=Path, default=DEFAULT_PIPELINE)
    ap.add_argument("--M", type=int, default=28)
    ap.add_argument("--r", type=float, default=6.7439020359331625)
    ap.add_argument("--Y0", type=float, default=0.8)
    ap.add_argument("--n-face", type=int, default=16)
    ap.add_argument("--json-out", type=Path, default=Path("dual_certification_audit_result.json"))
    ns = ap.parse_args()
    if not ns.pipeline.is_dir():
        raise SystemExit(f"pipeline not found: {ns.pipeline}")
    sys.path.insert(0, str(ns.pipeline))
    out = audit(ns.M, ns.r, ns.Y0, ns.n_face)
    ns.json_out.write_text(json.dumps(out, indent=2, allow_nan=False), encoding="utf-8")
    print(ns.json_out.resolve())
    print(json.dumps({"eta_inf": out["residual_by_cusp"]["infinity"]["delta_aut"],
                      "condition_equilibrated": out["conditioning"]["equilibrated"]["condition_2"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
