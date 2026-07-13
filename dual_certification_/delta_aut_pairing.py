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

Trial coefficient modes (morning path 1):
  --mode collocation   single-cusp σ0 collocation near-kernel (baseline)
  --mode multi         stacked multi-pairing jump operator near-kernel
  --mode periodize     multi near-kernel + finite Poincaré sum over gens

Language (honest):
  - Production *shape* of δ_aut (pairing faces, true pullback).
  - Not a certified Maass form. No hard-map flip. No Neumann-free upper bounds.

Usage:
  python delta_aut_pairing.py
  python delta_aut_pairing.py --M 48 --mode multi
  python delta_aut_pairing.py --M 48 --mode periodize --compare
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

# Finite set for Poincaré periodization.
# Horizontal generators first (preserve y); S is optional (height-mixing).
PERIOD_GROUP_H: List[Tuple[str, np.ndarray]] = [
    ("id", _mat(1.0, 0.0, 0.0, 1.0)),
    ("T1", GEN["T1"]),
    ("T1inv", _mat(1.0, -1.0, 0.0, 1.0)),
    ("R", GEN["R"]),
    ("TiR", GEN["TiR"]),
]
PERIOD_GROUP: List[Tuple[str, np.ndarray]] = PERIOD_GROUP_H + [
    ("S", GEN["S"]),
]


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


def mat_inv(mat: np.ndarray) -> np.ndarray:
    """Inverse of SL(2,ℂ) matrix (adjugate since det=1)."""
    a, b = complex(mat[0, 0]), complex(mat[0, 1])
    c, d = complex(mat[1, 0]), complex(mat[1, 1])
    return _mat(d, -b, -c, a)


def face_samples(name: str, n: int, Y0: float, y_min: float = 0.7071) -> np.ndarray:
    """Sample points on computational-cell faces paired by each generator."""
    n = max(4, int(n))
    g = max(2, int(math.ceil(math.sqrt(n))))
    u = (np.arange(g) + 0.37) / g
    pts: List[Tuple[float, float, float]] = []
    if name == "T1":
        for s in u:
            for t in u:
                pts.append((0.5, float(s - 0.5), float(y_min + (Y0 - y_min) * t)))
    elif name == "R":
        for s in u:
            for t in u:
                pts.append((float(s - 0.5), -0.5, float(y_min + (Y0 - y_min) * t)))
    elif name == "TiR":
        for s in u:
            for t in u:
                pts.append((float(s - 0.5), 0.5, float(y_min + (Y0 - y_min) * t)))
    elif name == "S":
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
        for s in u:
            for t in u:
                pts.append((float(s - 0.5), float(t - 0.5), Y0))
    return np.array(pts, dtype=np.float64)


def _mode_row(
    pts: np.ndarray,
    modes: Sequence[Tuple[int, int, int]],
    r: float,
    theta: float,
) -> np.ndarray:
    """V[j,k] = N^θ K_ir(2π√n y_j) exp(2πi (a x1 + b x2))."""
    n_pts = pts.shape[0]
    n_modes = len(modes)
    V = np.zeros((n_pts, n_modes), dtype=np.complex128)
    for k, (a, b, nn) in enumerate(modes):
        amp0 = (nn ** theta) if theta != 0 else 1.0
        for j in range(n_pts):
            x1, x2, y = pts[j]
            y_arg = 2.0 * math.pi * math.sqrt(nn) * max(float(y), 1e-12)
            w = amp0 * k_bessel_amp(y_arg, r)
            phase = np.exp(2j * math.pi * (a * x1 + b * x2))
            V[j, k] = w * phase
    return V


def eval_fourier(
    pts: np.ndarray,
    modes: Sequence[Tuple[int, int, int]],
    coeffs: np.ndarray,
    r: float,
    theta: float = 0.5,
) -> np.ndarray:
    """f(z,y) = Σ a_β N(β)^θ K_ir(2π|β|y) exp(2π i (a x1 + b x2))."""
    return _mode_row(pts, modes, r, theta) @ coeffs


def eval_fourier_periodized(
    pts: np.ndarray,
    modes: Sequence[Tuple[int, int, int]],
    coeffs: np.ndarray,
    r: float,
    theta: float = 0.5,
    group: Optional[Sequence[Tuple[str, np.ndarray]]] = None,
    y_floor: float = 0.5,
) -> np.ndarray:
    """
    Finite Poincaré sum: f_per(p) = Σ_{g ∈ F} f(g^{-1} · p).

    Default F = horizontal generators (id, T1±, R, TiR). Terms whose pullback
    height falls below y_floor are dropped (avoids cusp blow-up under S).
    """
    if group is None:
        group = PERIOD_GROUP_H
    n_pts = pts.shape[0]
    out = np.zeros(n_pts, dtype=np.complex128)
    counts = np.zeros(n_pts, dtype=np.float64)
    for _name, mat in group:
        inv = mat_inv(mat)
        pts_g = np.zeros_like(pts)
        ok = np.ones(n_pts, dtype=bool)
        for j in range(n_pts):
            z = complex(pts[j, 0], pts[j, 1])
            y = float(pts[j, 2])
            zg, yg = h3_act(inv, z, y)
            if yg < y_floor or not math.isfinite(yg):
                ok[j] = False
                continue
            pts_g[j, 0] = zg.real
            pts_g[j, 1] = zg.imag
            pts_g[j, 2] = yg
        if not np.any(ok):
            continue
        vals = eval_fourier(pts_g[ok], modes, coeffs, r, theta=theta)
        out[ok] = out[ok] + vals
        counts[ok] += 1.0
    counts = np.maximum(counts, 1.0)
    return out / counts


def mode_amplitudes(
    modes: Sequence[Tuple[int, int, int]],
    Y0: float,
    r: float,
    theta: float,
) -> np.ndarray:
    amps = np.zeros(len(modes), dtype=np.float64)
    for k, (_a, _b, nn) in enumerate(modes):
        y_arg = 2.0 * math.pi * math.sqrt(nn) * Y0
        amps[k] = (nn ** theta) * k_bessel_amp(y_arg, r)
    return np.maximum(amps, 1e-300)


def near_kernel_collocation(
    M: int,
    r: float,
    Y0: float,
    theta: float = 0.5,
    height_mismatch: float = 0.0,
) -> Tuple[np.ndarray, List[Tuple[int, int, int]], Dict[str, float]]:
    """Right singular vector of single-cusp collocation (baseline).

    height_mismatch=0: height-matched σ0 model (conditioning path).
    height_mismatch=1: true-height pullback weights (production-shaped σ0).
    """
    V, amps, modes = build_single_cusp_V(
        M, Y0=Y0, r=r, theta=theta, height_mismatch=height_mismatch
    )
    F = V / np.maximum(amps[None, :], 1e-300)
    Feq = equilibrate(F, n_iter=6)
    _u, s, vh = np.linalg.svd(Feq, full_matrices=False)
    sig_min = float(s[-1])
    sig_max = float(s[0]) if len(s) else 1.0
    u = vh[-1, :].conj()
    a = u / np.maximum(amps, 1e-300)
    a = a / max(np.linalg.norm(a), 1e-300)
    kind = "collocation_true_height" if height_mismatch > 0 else "collocation"
    return a.astype(np.complex128), modes, dict(
        kind=kind,
        sigma_min=sig_min,
        sigma_max=sig_max,
        rel=sig_min / max(sig_max, 1e-300),
        tau_disc=float(np.linalg.norm(V @ a)),
        n_modes=float(len(modes)),
        n_pts=float(V.shape[0]),
        kappa_eq=float(sig_max / max(sig_min, 1e-300)),
        height_mismatch=float(height_mismatch),
    )


def build_multi_pairing_V(
    M: int,
    r: float,
    Y0: float,
    theta: float = 0.5,
    n_face: int = 16,
    y_min: float = 1.0 / math.sqrt(2.0),
    gens: Sequence[str] = ("R", "TiR", "S", "T1"),
) -> Tuple[np.ndarray, np.ndarray, List[Tuple[int, int, int]]]:
    """
    Stacked jump collocation:

      V_γ[j, β] = φ_β(p_j) − φ_β(γ · p_j)

    over face samples of each generator γ. Near-kernel ⇒ small multi-pairing δ.
    """
    modes = gaussian_modes(M)
    blocks: List[np.ndarray] = []
    for name in gens:
        if name not in GEN:
            continue
        mat = GEN[name]
        pts = face_samples(name, n_face, Y0, y_min=y_min)
        V0 = _mode_row(pts, modes, r, theta)
        pts_g = np.zeros_like(pts)
        for j in range(pts.shape[0]):
            z = complex(pts[j, 0], pts[j, 1])
            y = float(pts[j, 2])
            zg, yg = h3_act(mat, z, y)
            pts_g[j, 0] = zg.real
            pts_g[j, 1] = zg.imag
            pts_g[j, 2] = max(yg, 1e-12)
        V1 = _mode_row(pts_g, modes, r, theta)
        blocks.append(V0 - V1)
    if not blocks:
        raise ValueError("no generator blocks")
    V = np.vstack(blocks)
    # normalize rows for scale-free residual
    V = V / math.sqrt(float(V.shape[0]))
    amps = mode_amplitudes(modes, Y0, r, theta)
    return V, amps, modes


def near_kernel_multi(
    M: int,
    r: float,
    Y0: float,
    theta: float = 0.5,
    n_face: int = 16,
    y_min: float = 1.0 / math.sqrt(2.0),
    jump_weight: float = 1.0,
    height_mismatch: float = 0.0,
) -> Tuple[np.ndarray, List[Tuple[int, int, int]], Dict[str, float]]:
    """
    Hybrid near-kernel: stack single-cusp collocation with multi-pairing jumps.

    Pure jump-only operators are wide (n_face·#gens ≪ n_modes) and have a huge
    null space → garbage δ_aut. Hybrid pins the trial with the Hejhal-like
    collocation residual while penalizing R/TiR/S/T1 face jumps.

    height_mismatch>0: use true-height σ0 collocation pin (production-shaped).
    """
    V_col, amps, modes = build_single_cusp_V(
        M, Y0=Y0, r=r, theta=theta, height_mismatch=height_mismatch
    )
    V_jump, amps_j, modes_j = build_multi_pairing_V(
        M, r, Y0, theta=theta, n_face=n_face, y_min=y_min
    )
    assert modes == modes_j
    # Frobenius-normalize blocks so weights are meaningful
    nc = max(np.linalg.norm(V_col, ord="fro"), 1e-300)
    nj = max(np.linalg.norm(V_jump, ord="fro"), 1e-300)
    V = np.vstack([V_col / nc, (jump_weight * V_jump) / nj])
    F = V / amps[None, :]
    Feq = equilibrate(F, n_iter=6)
    _u, s, vh = np.linalg.svd(Feq, full_matrices=False)
    sig_min = float(s[-1])
    sig_max = float(s[0]) if len(s) else 1.0
    u = vh[-1, :].conj()
    a = u / amps
    a = a / max(np.linalg.norm(a), 1e-300)
    kind = (
        "hybrid_true_height_multi"
        if height_mismatch > 0
        else "hybrid_collocation_multi"
    )
    return a.astype(np.complex128), modes, dict(
        kind=kind,
        sigma_min=sig_min,
        sigma_max=sig_max,
        rel=sig_min / max(sig_max, 1e-300),
        tau_disc=float(np.linalg.norm(V @ a)),
        n_modes=float(len(modes)),
        n_pts=float(V.shape[0]),
        n_col=float(V_col.shape[0]),
        n_jump=float(V_jump.shape[0]),
        kappa_eq=float(sig_max / max(sig_min, 1e-300)),
        jump_weight=float(jump_weight),
        height_mismatch=float(height_mismatch),
    )


def measure_jumps(
    coeffs: np.ndarray,
    modes: Sequence[Tuple[int, int, int]],
    r: float,
    Y0: float,
    theta: float,
    n_face: int,
    y_min: float,
    use_true_height: bool,
    periodize: bool,
) -> Tuple[Dict[str, Any], float, float]:
    """Return (per_gen, delta_aut, scale)."""
    g = max(4, int(math.ceil(math.sqrt(n_face))))
    xs = (np.arange(g) + 0.37) / g - 0.5
    X1, X2 = np.meshgrid(xs, xs, indexing="xy")
    ref_h = np.column_stack([X1.ravel(), X2.ravel(), np.full(X1.size, Y0)])
    eval_fn = eval_fourier_periodized if periodize else eval_fourier
    f_ref = eval_fn(ref_h, modes, coeffs, r, theta=theta)
    scale = float(np.sqrt(np.mean(np.abs(f_ref) ** 2)))
    scale = max(scale, 1e-30)

    per_gen: Dict[str, Any] = {}
    deltas: List[float] = []
    for name, mat in GEN.items():
        pts = face_samples(name, n_face, Y0, y_min=y_min)
        f0 = eval_fn(pts, modes, coeffs, r, theta=theta)
        pts_g = np.zeros_like(pts)
        for j in range(pts.shape[0]):
            z = complex(pts[j, 0], pts[j, 1])
            y = float(pts[j, 2])
            zg, yg = h3_act(mat, z, y)
            if not use_true_height:
                yg = Y0
            pts_g[j, 0] = zg.real
            pts_g[j, 1] = zg.imag
            pts_g[j, 2] = max(yg, 1e-12)
        f1 = eval_fn(pts_g, modes, coeffs, r, theta=theta)
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
    return per_gen, delta_aut, scale


def delta_aut_on_pairings(
    M: int = 48,
    r: float = 6.0,
    Y0: float = 0.8,
    theta: float = 0.5,
    n_face: int = 16,
    y_min: float = 1.0 / math.sqrt(2.0),
    use_true_height: bool = True,
    mode: str = "multi",
    jump_weight: float = 1.0,
) -> Dict[str, Any]:
    """
    Compute per-generator max jump |f − f∘γ| / scale on face samples.

    mode:
      collocation — single-cusp near-kernel (baseline, R/TiR often O(1))
      multi       — hybrid collocation + multi-pairing near-kernel
      periodize   — multi near-kernel + finite Poincaré evaluation
      sigma0      — hybrid with true-height σ0 collocation pin + jumps
      sigma0_per  — sigma0 hybrid + Poincaré evaluation
    """
    t0 = time.time()
    mode = mode.lower().strip()
    periodize = mode in ("periodize", "sigma0_per")
    true_h = 1.0 if mode in ("sigma0", "sigma0_per") else 0.0
    if mode in ("multi", "periodize", "sigma0", "sigma0_per"):
        coeffs, modes, meta = near_kernel_multi(
            M,
            r,
            Y0,
            theta=theta,
            n_face=n_face,
            y_min=y_min,
            jump_weight=jump_weight,
            height_mismatch=true_h,
        )
    else:
        coeffs, modes, meta = near_kernel_collocation(M, r, Y0, theta=theta)

    per_gen, delta_aut, scale = measure_jumps(
        coeffs, modes, r, Y0, theta, n_face, y_min, use_true_height, periodize
    )
    tau_proxy = float(meta["rel"]) * math.sqrt(1.0 + r * r)
    return dict(
        M=M,
        r=r,
        Y0=Y0,
        theta=theta,
        y_min=y_min,
        mode=mode,
        use_true_height=use_true_height,
        scale=scale,
        collocation=meta,
        per_generator=per_gen,
        delta_aut=delta_aut,
        tau_proxy=tau_proxy,
        eta_production_proxy=delta_aut + tau_proxy,
        language=(
            f"mode={mode}: δ_aut = max_γ max_face |f−f∘γ|/‖f‖ on PAIRINGS faces; "
            "true H³ pullback. Not certified Maass form. Not for hard-map flip."
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
    p.add_argument(
        "--mode",
        type=str,
        default="multi",
        choices=("collocation", "multi", "periodize", "sigma0", "sigma0_per"),
        help="trial mode (default multi; sigma0 = true-height hybrid)",
    )
    p.add_argument(
        "--jump-weight",
        type=float,
        default=1.0,
        help="relative weight of multi-pairing block in hybrid mode",
    )
    p.add_argument(
        "--compare",
        action="store_true",
        help="run collocation vs multi vs periodize side-by-side",
    )
    p.add_argument("--height-matched", action="store_true")
    p.add_argument("--json-out", type=str, default="")
    args = p.parse_args(argv)

    modes_run = (
        ["collocation", "multi", "periodize", "sigma0", "sigma0_per"]
        if args.compare
        else [args.mode]
    )
    results: Dict[str, Any] = {}
    print("=== δ_aut on PAIRINGS faces (production-path residual) ===")
    for mode in modes_run:
        out = delta_aut_on_pairings(
            M=args.M,
            r=args.r,
            Y0=args.Y0,
            theta=args.theta,
            n_face=args.n_face,
            use_true_height=not args.height_matched,
            mode=mode,
            jump_weight=args.jump_weight,
        )
        results[mode] = out
        col = out["collocation"]
        print(
            f"\n[{mode}] M={out['M']} r={out['r']} Y0={out['Y0']} "
            f"({out['seconds']:.2f}s)"
        )
        print(
            f"  operator: kind={col.get('kind')}  σ_min={col['sigma_min']:.3e}  "
            f"rel={col['rel']:.3e}  κ_eq={col['kappa_eq']:.1f}"
        )
        print(f"  δ_aut = {out['delta_aut']:.6e}  τ_proxy={out['tau_proxy']:.6e}  "
              f"η≈{out['eta_production_proxy']:.6e}")
        for name, rec in out["per_generator"].items():
            print(
                f"    {name:4s}  δ_max={rec['delta_max']:.4e}  "
                f"δ_rms={rec['delta_rms']:.4e}  n={rec['n_samples']}"
            )

    if args.compare and len(results) >= 2:
        base = results["collocation"]["delta_aut"]
        print("\n=== compare δ_aut factors vs collocation baseline ===")
        for mode, out in results.items():
            d = out["delta_aut"]
            fac = base / d if d > 0 else float("inf")
            print(f"  {mode:12s}  δ={d:.4e}  factor vs collocation={fac:.2f}×")

    payload = results if args.compare else results[modes_run[0]]
    path = Path(args.json_out) if args.json_out else _HERE / "delta_aut_result.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\n  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
