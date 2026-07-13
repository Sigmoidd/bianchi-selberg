#!/usr/bin/env python3
"""
Production-path automorphy defect on PAIRINGS faces (Gaussian / Γ₀(2+i)).

Measures a *true jump residual* for a truncated Fourier trial form:

  δ_aut(γ) ≈ max_{p ∈ face samples} |f(p) − f(γ·p)| / scale

using the four generators of independent_exclusion.PAIRINGS / congruence:

  T1 = [[1, 1], [0, 1]]
  R  = [[i, 0], [0, −i]]
  TiR = Ti · R with Ti = [[1, i], [0, 1]]
  S  = [[0, −1], [1, 0]]

and the standard SL(2,ℂ) action on H³ = ℂ × ℝ₊.

Language (honest):
  - This is the production *shape* of δ_aut (pairing faces, true pullback).
  - Coefficients from the collocation near-kernel are still a *model* trial,
    not a certified Maass form.
  - No claim of η ≤ η₀. No hard-map flip. No upper bound from Neumann-free P1.

Usage:
  python delta_aut_pairing.py
  python delta_aut_pairing.py --M 64 --r 6.0 --Y0 0.8 --n-face 24
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from hejhal_conditioning import (  # noqa: E402
    build_single_cusp_V,
    gaussian_modes,
    k_bessel_amp,
    equilibrate,
)


# ---------------------------------------------------------------------------
# Generators (complex entries) — match congruence_prototype.generator_matrices
# ---------------------------------------------------------------------------

def _mat(a: complex, b: complex, c: complex, d: complex) -> np.ndarray:
    return np.array([[a, b], [c, d]], dtype=np.complex128)


GEN: Dict[str, np.ndarray] = {
    "T1": _mat(1.0, 1.0, 0.0, 1.0),
    "R": _mat(1j, 0.0, 0.0, -1j),
    "TiR": _mat(1.0, 1j, 0.0, 1.0) @ _mat(1j, 0.0, 0.0, -1j),  # Ti · R
    "S": _mat(0.0, -1.0, 1.0, 0.0),
}


def h3_act(mat: np.ndarray, z: complex, y: float) -> Tuple[complex, float]:
    """
    SL(2,ℂ) action on upper half-space H³ ≅ ℂ × ℝ₊.

      denom = |c z + d|² + |c|² y²
      z'    = ((a z + b) conj(c z + d) + a conj(c) y²) / denom
      y'    = y / denom
    """
    a, b = complex(mat[0, 0]), complex(mat[0, 1])
    c, d = complex(mat[1, 0]), complex(mat[1, 1])
    cz_d = c * z + d
    denom = abs(cz_d) ** 2 + (abs(c) ** 2) * (y * y)
    if denom < 1e-30:
        return 0j, 1e30
    z_new = ((a * z + b) * np.conj(cz_d) + a * np.conj(c) * (y * y)) / denom
    y_new = y / denom
    return complex(z_new), float(y_new)


def face_samples(name: str, n: int, Y0: float, y_min: float = 0.7071) -> np.ndarray:
    """
    Sample points on computational-cell faces paired by each generator.

    Coordinates: (x1, x2, y) with z = x1 + i x2.
    Faces follow independent_exclusion PAIRINGS tags (x1m/x1p, x2m/x2p, floor).
    """
    n = max(4, int(n))
    g = max(2, int(math.ceil(math.sqrt(n))))
    u = (np.arange(g) + 0.37) / g
    pts: List[Tuple[float, float, float]] = []
    if name == "T1":
        # RIGHT face x1 = 1/2 (paired to LEFT x1=0 via T1)
        for s in u:
            for t in u:
                pts.append((0.5, float(s - 0.5), float(y_min + (Y0 - y_min) * t)))
    elif name == "R":
        # x2 = −1/2 face
        for s in u:
            for t in u:
                pts.append((float(s - 0.5), -0.5, float(y_min + (Y0 - y_min) * t)))
    elif name == "TiR":
        # x2 = +1/2 face
        for s in u:
            for t in u:
                pts.append((float(s - 0.5), 0.5, float(y_min + (Y0 - y_min) * t)))
    elif name == "S":
        # floor hemisphere |z|² + y² = 1, y > 0, |x1|,|x2| ≤ 1/2
        for s in u:
            for t in u:
                x1 = float(s - 0.5)
                x2 = float(t - 0.5)
                r2 = x1 * x1 + x2 * x2
                if r2 >= 1.0 - 1e-12:
                    continue
                y = math.sqrt(max(1.0 - r2, 0.0))
                if y < y_min * 0.5:
                    continue
                pts.append((x1, x2, y))
    else:
        raise ValueError(f"unknown face {name}")
    if not pts:
        # fallback: horosphere samples
        for s in u:
            for t in u:
                pts.append((float(s - 0.5), float(t - 0.5), Y0))
    return np.array(pts, dtype=np.float64)


def eval_fourier(
    pts: np.ndarray,
    modes: Sequence[Tuple[int, int, int]],
    coeffs: np.ndarray,
    r: float,
    theta: float = 0.5,
) -> np.ndarray:
    """f(z,y) = Σ a_β N(β)^θ K_ir(2π|β|y) exp(2π i (a x1 + b x2))."""
    n_pts = pts.shape[0]
    out = np.zeros(n_pts, dtype=np.complex128)
    for k, (a, b, nn) in enumerate(modes):
        ck = coeffs[k]
        if abs(ck) < 1e-30:
            continue
        amp0 = (nn ** theta) if theta != 0 else 1.0
        for j in range(n_pts):
            x1, x2, y = pts[j]
            y_arg = 2.0 * math.pi * math.sqrt(nn) * max(y, 1e-12)
            w = amp0 * k_bessel_amp(y_arg, r)
            phase = np.exp(2j * math.pi * (a * x1 + b * x2))
            out[j] += ck * w * phase
    return out


def near_kernel_coeffs(
    M: int,
    r: float,
    Y0: float,
    theta: float = 0.5,
) -> Tuple[np.ndarray, List[Tuple[int, int, int]], Dict[str, float]]:
    """
    Right singular vector of equilibrated single-cusp collocation V for σ_min.
    Returns coeffs (unit ℓ²), modes, and residual diagnostics.
    """
    V, amps, modes = build_single_cusp_V(M, Y0=Y0, r=r, theta=theta)
    # right amp precond then equilibration (same hygiene as conditioning path)
    F = V / np.maximum(amps[None, :], 1e-300)
    Feq = equilibrate(F, n_iter=6)
    # SVD: V a ≈ 0  ⇔  Feq (D a) ≈ 0 with a = D^{-1} u
    _u, s, vh = np.linalg.svd(Feq, full_matrices=False)
    sig_min = float(s[-1])
    sig_max = float(s[0]) if len(s) else 1.0
    u = vh[-1, :].conj()  # right singular vector in equilibrated coords
    # map back: Feq = Dr F Dc roughly unit; use amp-scaled coefficients
    a = u / np.maximum(amps, 1e-300)
    nrm = np.linalg.norm(a)
    if nrm > 0:
        a = a / nrm
    # discrete residual τ_disc = ||V a||₂ / ||a||_*  (here ||a||_2=1 after scale)
    Va = V @ a
    tau_disc = float(np.linalg.norm(Va))
    rel = sig_min / max(sig_max, 1e-300)
    meta = dict(
        sigma_min=sig_min,
        sigma_max=sig_max,
        rel=rel,
        tau_disc=tau_disc,
        n_modes=float(len(modes)),
        n_pts=float(V.shape[0]),
        kappa_eq=float(sig_max / max(sig_min, 1e-300)),
    )
    return a.astype(np.complex128), modes, meta


def delta_aut_on_pairings(
    M: int = 48,
    r: float = 6.0,
    Y0: float = 0.8,
    theta: float = 0.5,
    n_face: int = 16,
    y_min: float = 1.0 / math.sqrt(2.0),
    use_true_height: bool = True,
) -> Dict[str, Any]:
    """
    Compute per-generator max jump |f − f∘γ| / scale on face samples.

    use_true_height: evaluate f at true pulled-back height (production path).
    If False, clamp pullback height to Y0 (height-matched diagnostic).
    """
    t0 = time.time()
    coeffs, modes, meta = near_kernel_coeffs(M, r, Y0, theta=theta)
    # scale: RMS of |f| on a reference horosphere
    ref = face_samples("T1", n_face, Y0, y_min=y_min)
    # use horosphere y=Y0 strip for scale
    ref_h = np.column_stack(
        [
            (np.linspace(-0.5, 0.5, max(4, n_face // 2), endpoint=False)),
            (np.linspace(-0.5, 0.5, max(4, n_face // 2), endpoint=False)),
            np.full(max(4, n_face // 2), Y0),
        ]
    )
    # proper grid
    g = max(4, int(math.ceil(math.sqrt(n_face))))
    xs = (np.arange(g) + 0.37) / g - 0.5
    X1, X2 = np.meshgrid(xs, xs, indexing="xy")
    ref_h = np.column_stack([X1.ravel(), X2.ravel(), np.full(X1.size, Y0)])
    f_ref = eval_fourier(ref_h, modes, coeffs, r, theta=theta)
    scale = float(np.sqrt(np.mean(np.abs(f_ref) ** 2)))
    scale = max(scale, 1e-30)

    per_gen: Dict[str, Any] = {}
    deltas = []
    for name, mat in GEN.items():
        pts = face_samples(name, n_face, Y0, y_min=y_min)
        f0 = eval_fourier(pts, modes, coeffs, r, theta=theta)
        pts_g = np.zeros_like(pts)
        for j in range(pts.shape[0]):
            z = complex(pts[j, 0], pts[j, 1])
            y = float(pts[j, 2])
            zg, yg = h3_act(mat, z, y)
            if not use_true_height:
                yg = Y0
            pts_g[j, 0] = zg.real
            pts_g[j, 1] = zg.imag
            pts_g[j, 2] = yg
        f1 = eval_fourier(pts_g, modes, coeffs, r, theta=theta)
        jump = np.abs(f0 - f1)
        d_max = float(np.max(jump) / scale)
        d_rms = float(np.sqrt(np.mean(jump ** 2)) / scale)
        per_gen[name] = dict(
            n_samples=int(pts.shape[0]),
            delta_max=d_max,
            delta_rms=d_rms,
            jump_max_abs=float(np.max(jump)),
        )
        deltas.append(d_max)

    delta_aut = float(max(deltas)) if deltas else float("inf")
    # honest production residual package
    tau_disc = float(meta["tau_disc"])
    eta_prod = delta_aut + tau_disc
    return dict(
        M=M,
        r=r,
        Y0=Y0,
        theta=theta,
        y_min=y_min,
        use_true_height=use_true_height,
        scale=scale,
        collocation=meta,
        per_generator=per_gen,
        delta_aut=delta_aut,
        tau_disc=tau_disc,
        eta_production_proxy=eta_prod,
        language=(
            "δ_aut = max_γ max_face |f−f∘γ|/||f||_horosphere from truncated "
            "Fourier near-kernel; τ_disc = ||V a||₂ with ||a||₂=1. "
            "Model trial, not certified Maass form. Not for hard-map flip."
        ),
        seconds=time.time() - t0,
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--M", type=int, default=48)
    p.add_argument("--r", type=float, default=6.0)
    p.add_argument("--Y0", type=float, default=0.8)
    p.add_argument("--theta", type=float, default=0.5)
    p.add_argument("--n-face", type=int, default=16)
    p.add_argument("--height-matched", action="store_true",
                   help="clamp pullback height to Y0 (diagnostic)")
    p.add_argument("--json-out", type=str, default="")
    args = p.parse_args(argv)

    print("=== δ_aut on PAIRINGS faces (production-path residual) ===")
    out = delta_aut_on_pairings(
        M=args.M,
        r=args.r,
        Y0=args.Y0,
        theta=args.theta,
        n_face=args.n_face,
        use_true_height=not args.height_matched,
    )
    print(f"M={out['M']}  r={out['r']}  Y0={out['Y0']}  "
          f"true_height={out['use_true_height']}  ({out['seconds']:.2f}s)")
    col = out["collocation"]
    print(f"  collocation: σ_min={col['sigma_min']:.3e}  rel={col['rel']:.3e}  "
          f"κ_eq={col['kappa_eq']:.1f}  τ_disc={col['tau_disc']:.3e}")
    print(f"  δ_aut (max over generators) = {out['delta_aut']:.6e}")
    print(f"  τ_disc                      = {out['tau_disc']:.6e}")
    print(f"  η_proxy = δ_aut + τ_disc    = {out['eta_production_proxy']:.6e}")
    print("  per generator:")
    for name, rec in out["per_generator"].items():
        print(f"    {name:4s}  δ_max={rec['delta_max']:.4e}  "
              f"δ_rms={rec['delta_rms']:.4e}  n={rec['n_samples']}")
    print(f"  language: {out['language']}")

    if args.json_out:
        path = Path(args.json_out)
        path.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"  wrote {path}")
    else:
        default = _HERE / "delta_aut_result.json"
        default.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"  wrote {default}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
