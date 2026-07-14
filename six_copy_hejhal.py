#!/usr/bin/env python3
"""Exact-coset six-copy Hejhal collocation for Gamma_0(2+i).

The reference-cell generators act by exact permutations of P^1(F_5).  If
gamma_c are right-coset representatives and pi_delta is the gluing
permutation, a subgroup-automorphic function must satisfy

    F_c(p) = F_{pi_delta(c)}(delta p),   F_c = f o gamma_c.

For the ordering used by independent_exclusion.congruence_prototype,
gamma_0 = I and gamma_k = S T_{k-1} for k=1,...,5.  Thus F_0 uses the
infinity-cusp expansion, while the other five components use translated
copies of g=f o S at cusp zero.  The cusp-zero dual lattice is
(1/conj(2+i)) Z[i] = ((2+i)/5) Z[i].

All K_{ir} midpoint/radius pairs come from Arb.  Phases, sample locations,
SVD, and normalization remain floating-point, so this is a diagnostic
collocation experiment, not yet a continuum residual certificate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any, Sequence

import numpy as np

from verified_kir import kir_enclosure
from verified_hejhal_phase1 import equilibrate_with_ledger, svd_record


DEFAULT_PIPELINE = Path(r"C:\Users\Admin\.grok\worktrees\bearings-bianchi-selberg\np9\dual_certification_")
ROOT = Path(__file__).resolve().parent
P = 2.0 + 1.0j
DUAL_SCALE = P / 5.0  # 1/conj(P)


def gaussian_modes(M: int) -> list[tuple[int, int, int]]:
    out = []
    s = math.isqrt(M)
    for a in range(-s, s + 1):
        for b in range(-s, s + 1):
            nn = a * a + b * b
            if 0 < nn <= M:
                out.append((a, b, nn))
    return sorted(out, key=lambda q: (q[2], q[0], q[1]))


def mode_frequencies(modes: Sequence[tuple[int, int, int]], cusp: int):
    beta = np.array([complex(a, b) for a, b, _ in modes], dtype=np.complex128)
    return beta if cusp == 0 else DUAL_SCALE * beta


def basis_mid_rad(
    pts: np.ndarray,
    modes: Sequence[tuple[int, int, int]],
    r: float,
    cusp: int,
    translate: float,
    bits: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Whittaker basis y K_ir(2 pi |mu| y) exp(2 pi i <mu,z>)."""
    freq = mode_frequencies(modes, cusp)
    npt, nm = pts.shape[0], len(modes)
    mid = np.empty((npt, nm), dtype=np.complex128)
    rad = np.empty((npt, nm), dtype=float)
    x1 = pts[:, 0] + float(translate)
    x2 = pts[:, 1]
    y = pts[:, 2]
    for k, mu in enumerate(freq):
        phase = np.exp(2j * math.pi * (mu.real * x1 + mu.imag * x2))
        for j in range(npt):
            q = kir_enclosure(r, 2.0 * math.pi * abs(mu) * float(y[j]), bits)
            mid[j, k] = float(y[j]) * q.midpoint * phase[j]
            rad[j, k] = abs(float(y[j])) * q.radius
    return mid, rad


def component_mid_rad(
    copy: int,
    pts: np.ndarray,
    modes_inf: Sequence[tuple[int, int, int]],
    modes_0: Sequence[tuple[int, int, int]],
    r: float,
    bits: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return rows in the common coefficient vector [a_inf,a_0]."""
    ni, n0 = len(modes_inf), len(modes_0)
    out = np.zeros((pts.shape[0], ni + n0), dtype=np.complex128)
    rad = np.zeros((pts.shape[0], ni + n0), dtype=float)
    if copy == 0:
        q, qr = basis_mid_rad(
            pts, modes_inf, r, cusp=0, translate=0.0, bits=bits
        )
        out[:, :ni], rad[:, :ni] = q, qr
    else:
        q, qr = basis_mid_rad(
            pts, modes_0, r, cusp=1, translate=float(copy - 1), bits=bits
        )
        out[:, ni:], rad[:, ni:] = q, qr
    return out, rad


def transformed_points(pairing, name: str, pts: np.ndarray) -> np.ndarray:
    mat = pairing.GEN[name]
    out = np.zeros_like(pts)
    for j, (x1, x2, y) in enumerate(pts):
        z, yp = pairing.h3_act(mat, complex(x1, x2), float(y))
        out[j] = (z.real, z.imag, yp)
    return out


def interior_samples(n: int, Y0: float, M: int) -> np.ndarray:
    """Deterministic 3-D cloud where the coset automorphy identity holds."""
    g0 = max(2, int(math.ceil(max(int(n), 8) ** (1.0 / 3.0))))
    # Respect the Fourier Nyquist rate.  The earlier g0^3 grid aliased modes
    # when sqrt(M) >= g0, producing a false numerical nullspace for every r.
    gxy = max(2 * math.isqrt(M) + 2, g0)
    gy = max(3, g0)
    ux = (np.arange(gxy) + 0.37) / gxy
    uy = (np.arange(gy) + 0.37) / gy
    xs = 0.9 * (ux - 0.5)
    y_lo = 1.0 / math.sqrt(2.0) + 0.015
    y_hi = max(float(Y0), 1.15)
    ys = y_lo + (y_hi - y_lo) * uy
    return np.array([(x1, x2, y) for x1 in xs for x2 in xs for y in ys], dtype=float)


def build_operator(
    pipeline: Path, M: int, r: float, Y0: float, n_face: int, bits: int,
    sample_kind: str = "interior",
) -> dict[str, Any]:
    if str(pipeline) not in sys.path:
        sys.path.insert(0, str(pipeline))
    ie = pipeline.parent / "independent_exclusion"
    if str(ie) not in sys.path:
        sys.path.insert(0, str(ie))
    import delta_aut_pairing as pairing
    import congruence_prototype as cp

    points, glue, cusp_class = cp.build_gluing("(2+i)")
    expected = [(0, 1), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4)]
    if points != expected or cusp_class != [0, 1, 1, 1, 1, 1]:
        raise RuntimeError(f"unexpected F5 coset ordering: {points}, {cusp_class}")
    modes_inf = gaussian_modes(M)
    # |beta/conj(p)|^2 = N(beta)/N(p).  To give both cusps the same
    # physical frequency cutoff, the zero-cusp Gaussian index must extend to
    # N(p) M = 5M.
    modes_0 = gaussian_modes(5 * M)
    blocks, radii, tags = [], [], []
    for name in ("T1", "R", "TiR", "S"):
        if sample_kind == "face":
            pts = pairing.face_samples(
                name, n_face, Y0, y_min=1.0 / math.sqrt(2.0)
            )
        elif sample_kind == "interior":
            pts = interior_samples(n_face, Y0, M)
        else:
            raise ValueError(f"unknown sample_kind={sample_kind!r}")
        tgt = transformed_points(pairing, name, pts)
        for copy in range(6):
            other = int(glue[name][copy])
            left, lr = component_mid_rad(
                copy, pts, modes_inf, modes_0, r, bits
            )
            right, rr = component_mid_rad(
                other, tgt, modes_inf, modes_0, r, bits
            )
            blocks.append(left - right)
            radii.append(lr + rr)
            tags.extend((name, copy, other, j) for j in range(pts.shape[0]))
    return {
        "mid": np.vstack(blocks), "rad": np.vstack(radii), "tags": tags,
        "modes_inf": modes_inf, "modes_0": modes_0,
        "glue": glue, "points": points,
    }


def reference_matrix(
    pipeline: Path, modes_inf, modes_0, r: float, Y0: float, bits: int
):
    xs = (np.arange(6) + 0.37) / 6.0 - 0.5
    X1, X2 = np.meshgrid(xs, xs, indexing="xy")
    pts = np.column_stack([X1.ravel(), X2.ravel(), np.full(X1.size, Y0)])
    rows = []
    for copy in range(6):
        q, _ = component_mid_rad(copy, pts, modes_inf, modes_0, r, bits)
        rows.append(q)
    return np.vstack(rows)


def solve_operator(op: dict[str, Any], pipeline: Path, r: float, Y0: float, bits: int):
    A, Ar = op["mid"], op["rad"]
    modes_inf, modes_0 = op["modes_inf"], op["modes_0"]
    ni, n0 = len(modes_inf), len(modes_0)
    # Right amplitude scaling at Y0 for each actual cusp lattice.
    amp = np.empty(ni + n0, dtype=float)
    for cusp, modes, offset in ((0, modes_inf, 0), (1, modes_0, ni)):
        freq = mode_frequencies(modes, cusp)
        for k, mu in enumerate(freq):
            q = kir_enclosure(r, 2 * math.pi * abs(mu) * Y0, bits)
            amp[offset + k] = max(abs(Y0 * q.midpoint), q.radius, 1e-300)
    pre = A / amp[None, :]
    # Exact self-identifications can give identically zero equations (for
    # example a translation already built into a cusp Fourier expansion).
    # They carry no information and make multiplicative row equilibration
    # overflow, so remove only rows that are zero to roundoff after amplitude
    # scaling.  Keep them in the final residual report as zero constraints.
    row_peak = np.max(np.abs(pre), axis=1)
    active = row_peak > max(float(np.max(row_peak)) * 1e-14, 1e-280)
    pre_active = pre[active]
    if pre_active.shape[0] < pre_active.shape[1]:
        raise ArithmeticError(
            f"active collocation operator is underdetermined: {pre_active.shape}"
        )
    col = np.maximum(np.linalg.norm(pre_active, axis=0), 1e-300)
    G = pre_active / col[None, :]
    Aeq, left, right, ledger_error = equilibrate_with_ledger(G, 6)
    _, s, vh = np.linalg.svd(Aeq, full_matrices=False)
    veq = vh[-1].conj()
    coeff = (right * veq) / (amp * col)
    # Normalize by the actual six-component field on a reference horosphere.
    Ref = reference_matrix(pipeline, modes_inf, modes_0, r, Y0, bits)
    vals = Ref @ coeff
    scale = float(np.sqrt(np.mean(np.abs(vals) ** 2)))
    if not math.isfinite(scale) or scale <= 1e-300:
        raise ArithmeticError(f"degenerate reference normalization: {scale}")
    coeff /= scale
    vals = Ref @ coeff
    scale_check = float(np.sqrt(np.mean(np.abs(vals) ** 2)))

    residual = A @ coeff
    interval_row_upper = np.abs(residual) + Ar @ np.abs(coeff)
    per_relation = {}
    for name in ("T1", "R", "TiR", "S"):
        idx = np.array([tag[0] == name for tag in op["tags"]], dtype=bool)
        per_relation[name] = {
            "n_rows": int(idx.sum()),
            "max_mid_abs": float(np.max(np.abs(residual[idx]))),
            "rms_mid_abs": float(np.sqrt(np.mean(np.abs(residual[idx]) ** 2))),
            "max_entry_interval_upper": float(np.max(interval_row_upper[idx])),
        }
    physical_rel = np.linalg.norm(residual)
    physical_rel /= max(np.linalg.norm(A, 2) * np.linalg.norm(coeff), 1e-300)
    a_inf, a_0 = coeff[:ni], coeff[ni:]
    index_inf = {(a, b): k for k, (a, b, _nn) in enumerate(modes_inf)}
    lift_pairs = []
    off_zero = []
    for k, (a, b, _nn) in enumerate(modes_0):
        # mu=(a+bi)*(2+i)/5 is integral exactly when both numerators divide 5.
        ur, ui = 2 * a - b, a + 2 * b
        if ur % 5 == 0 and ui % 5 == 0 and (ur // 5, ui // 5) in index_inf:
            lift_pairs.append((index_inf[(ur // 5, ui // 5)], k))
        else:
            off_zero.append(k)
    if lift_pairs:
        vinf = np.array([a_inf[i] for i, _ in lift_pairs])
        vzero = np.array([a_0[k] for _, k in lift_pairs])
        scalar = np.vdot(vinf, vzero) / max(np.vdot(vinf, vinf).real, 1e-300)
        lift_mismatch = np.linalg.norm(vzero - scalar * vinf)
        lift_mismatch /= max(np.linalg.norm(vzero), 1e-300)
    else:
        scalar, lift_mismatch = 0j, float("inf")
    zero_norm = max(np.linalg.norm(a_0), 1e-300)
    off_lattice_fraction = np.linalg.norm(a_0[off_zero]) / zero_norm
    return {
        "coeff": coeff,
        "record": {
            "physical": svd_record(A), "preconditioned": svd_record(pre),
            "column_normalized": svd_record(G), "equilibrated": svd_record(Aeq),
            "equilibration_ledger_relative_error": ledger_error,
            "reference_rms_after_normalization": scale_check,
            "physical_relative_residual": float(physical_rel),
            "sampled_residual_l2": float(np.linalg.norm(residual)),
            "sampled_residual_rms": float(np.sqrt(np.mean(np.abs(residual) ** 2))),
            "sampled_residual_max": float(np.max(np.abs(residual))),
            "sampled_interval_row_upper_max": float(np.max(interval_row_upper)),
            "per_relation": per_relation,
            "sigma_values_tail": [float(x) for x in s[-8:]],
            "amplitude_dynamic_range": float(amp.max() / amp.min()),
            "coefficient_norm_2": float(np.linalg.norm(coeff)),
            "n_rows_total": int(A.shape[0]),
            "n_rows_active": int(active.sum()),
            "n_rows_structurally_zero": int((~active).sum()),
            "left_scale_range": [float(left.min()), float(left.max())],
            "right_scale_range": [float(right.min()), float(right.max())],
            "lifted_level1_structure": {
                "n_matched_integer_zero_cusp_modes": len(lift_pairs),
                "zero_cusp_off_integer_lattice_fraction": float(off_lattice_fraction),
                "matched_coefficient_relative_mismatch_after_scalar": float(lift_mismatch),
                "best_scalar_real": float(np.real(scalar)),
                "best_scalar_imag": float(np.imag(scalar)),
            },
        },
    }


def one_point(
    pipeline: Path, M: int, r: float, Y0: float, n_face: int, bits: int,
    sample_kind: str = "interior",
):
    op = build_operator(pipeline, M, r, Y0, n_face, bits, sample_kind)
    sol = solve_operator(op, pipeline, r, Y0, bits)
    return {
        "r": r, "n_rows": int(op["mid"].shape[0]),
        "n_coefficients": int(op["mid"].shape[1]),
        "n_modes_infinity": len(op["modes_inf"]),
        "n_modes_zero": len(op["modes_0"]),
        "coefficient_order": "[infinity Gaussian N<=M, zero Gaussian N<=5M]",
        "coefficients": [
            {"real": repr(float(z.real)), "imag": repr(float(z.imag))}
            for z in sol["coeff"]
        ],
        **sol["record"],
    }


def parse_scan(text: str):
    a, b, n = text.split(":")
    return np.linspace(float(a), float(b), int(n))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline", type=Path, default=DEFAULT_PIPELINE)
    ap.add_argument("--M", type=int, default=28)
    ap.add_argument("--Y0", type=float, default=0.8)
    ap.add_argument("--n-face", type=int, default=16)
    ap.add_argument("--bits", type=int, default=160)
    ap.add_argument("--r-scan", default="5.5:8.0:11")
    ap.add_argument("--samples", choices=("interior", "face"), default="interior")
    ap.add_argument("--json-out", type=Path, default=Path("six_copy_hejhal_result.json"))
    ns = ap.parse_args()
    rows = []
    for r in parse_scan(ns.r_scan):
        print(f"r={r:.8f}", flush=True)
        rows.append(one_point(
            ns.pipeline, ns.M, float(r), ns.Y0, ns.n_face, ns.bits, ns.samples
        ))
    best = min(rows, key=lambda q: q["sampled_residual_rms"])
    out = {
        "status": "exact F5 coset gluing; Arb K values; floating sampled collocation; not certified eta",
        "parameters": {"M": ns.M, "Y0": ns.Y0, "n_face": ns.n_face,
                       "bits": ns.bits, "r_scan": ns.r_scan},
        "sample_kind": ns.samples,
        "cusp_lattices": {"infinity": "Z[i]", "zero_dual": "(2+i)/5 Z[i]"},
        "rows": rows, "best": best,
        "hard_map_changed": False,
        "remaining_for_certificate": [
            "interval phases and sample coordinates",
            "continuum face-norm enclosure between samples",
            "Fourier truncation and reprojection bounds in the D(K) norm",
            "interval singular/eigenparameter isolation",
            "counting and multi-copy CR bridge",
        ],
    }
    ns.json_out.write_text(json.dumps(out, indent=2, allow_nan=False), encoding="utf-8")
    print(ns.json_out.resolve())
    print(json.dumps({
        "best_r": best["r"], "best_rms": best["sampled_residual_rms"],
        "best_max": best["sampled_residual_max"],
        "best_interval_sample_upper": best["sampled_interval_row_upper_max"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
