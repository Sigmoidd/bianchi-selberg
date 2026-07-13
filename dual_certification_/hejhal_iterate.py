#!/usr/bin/env python3
"""
Production-shaped Hejhal *iterate*: alternate r-refine and coefficient update.

Not a full nonlinear Hejhal fixed-point for true automorphy. Engineering loop:

  1. Build two-cusp collocation at r; take near-kernel (a_∞, a_0).
  2. Measure rel = σ_min/σ_max and PAIRINGS δ_aut (optional periodize).
  3. Locally minimize rel over r (golden section on collocation residual only).
  4. Recompute near-kernel at best r; re-measure δ_aut.
  5. Optional: periodize Fourier on ∞ cusp samples and L²-project back to
     mode coefficients (one re-projection step).

Tracks residual history. Hard map unchanged.

Usage:
  python hejhal_iterate.py
  python hejhal_iterate.py --M 32 --iters 4 --r0 6.0 --periodize
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

from two_cusp_hejhal_N5 import (  # noqa: E402
    build_block_system,
    coupled_mid_rad,
    precondition_right,
)
from hejhal_conditioning import equilibrate  # noqa: E402
from production_hejhal_residual import two_cusp_near_kernel, residual_at_r  # noqa: E402
from delta_aut_pairing import (  # noqa: E402
    measure_jumps,
    gaussian_modes,
    eval_fourier,
    eval_fourier_periodized,
    face_samples,
)


def collocation_rel(M: int, r: float, Y0: float, theta: float = 0.5) -> float:
    """Fast relative residual of two-cusp operator (no δ_aut)."""
    sysb = build_block_system(M, float(r), Y0, theta=theta)
    mid, rad = coupled_mid_rad(sysb)
    mid_p, _, _ = precondition_right(mid, rad, sysb.w_inf, sysb.w_0)
    col = np.maximum(np.linalg.norm(mid_p, axis=0), 1e-300)
    G = mid_p / col[None, :]
    Geq = equilibrate(G, n_iter=4)
    s = np.linalg.svd(Geq, compute_uv=False)
    return float(s[-1] / max(s[0], 1e-300))


def golden_min_rel(
    M: int,
    r_lo: float,
    r_hi: float,
    Y0: float,
    n_iter: int = 12,
    theta: float = 0.5,
) -> Tuple[float, float]:
    """Golden-section minimize collocation rel on [r_lo, r_hi]."""
    phi = (math.sqrt(5) - 1) / 2
    a, b = r_lo, r_hi
    for _ in range(n_iter):
        r1 = b - phi * (b - a)
        r2 = a + phi * (b - a)
        e1 = collocation_rel(M, r1, Y0, theta)
        e2 = collocation_rel(M, r2, Y0, theta)
        if e1 < e2:
            b = r2
        else:
            a = r1
    r_star = 0.5 * (a + b)
    return r_star, collocation_rel(M, r_star, Y0, theta)


def reproject_periodized(
    a_inf: np.ndarray,
    modes: List[Tuple[int, int, int]],
    r: float,
    Y0: float,
    theta: float = 0.5,
    g: int = 8,
) -> np.ndarray:
    """
    Evaluate periodized f on a torus grid, least-squares fit free Fourier modes.

    Projects the periodized field back to the truncated mode basis (one step
    toward an automorphic coefficient vector).
    """
    xs = (np.arange(g) + 0.37) / g - 0.5
    X1, X2 = np.meshgrid(xs, xs, indexing="xy")
    pts = np.column_stack([X1.ravel(), X2.ravel(), np.full(X1.size, Y0)])
    f = eval_fourier_periodized(pts, modes, a_inf, r, theta=theta)
    # design matrix: mode evaluations at pts
    from delta_aut_pairing import _mode_row

    Phi = _mode_row(pts, modes, r, theta)
    # least squares Phi a ≈ f
    a_new, *_ = np.linalg.lstsq(Phi, f, rcond=None)
    nrm = np.linalg.norm(a_new)
    if nrm > 0:
        a_new = a_new / nrm
    return a_new.astype(np.complex128)


def iterate(
    M: int = 32,
    r0: float = 6.0,
    Y0: float = 0.8,
    n_outer: int = 4,
    r_halfwidth: float = 0.75,
    periodize: bool = True,
    reproject: bool = True,
    n_face: int = 16,
) -> Dict[str, Any]:
    t0 = time.time()
    history: List[Dict[str, Any]] = []
    r = float(r0)
    print(
        f"=== Hejhal iterate  M={M} r0={r0} Y0={Y0}  "
        f"outer={n_outer} periodize={periodize} reproject={reproject} ==="
    )

    a_inf = a_0 = None
    modes = gaussian_modes(M)

    for it in range(n_outer):
        print(f"\n--- outer iter {it}  r={r:.6f} ---", flush=True)
        # 1) r refine on collocation rel
        r_lo, r_hi = max(1.0, r - r_halfwidth), r + r_halfwidth
        r_star, rel_star = golden_min_rel(M, r_lo, r_hi, Y0, n_iter=10)
        print(f"  r-refine: r*={r_star:.6f}  rel={rel_star:.4e}", flush=True)
        r = r_star

        # 2) near-kernel at r*
        a_inf, a_0, meta = two_cusp_near_kernel(M, r, Y0)
        # 3) optional reproject periodized
        if reproject and periodize:
            a_try = reproject_periodized(a_inf, modes, r, Y0)
            # keep if it doesn't explode
            if np.isfinite(a_try).all():
                a_inf = a_try

        # 4) full residual measure
        per_gen, delta_aut, scale = measure_jumps(
            a_inf,
            modes,
            r,
            Y0,
            0.5,
            n_face,
            y_min=1.0 / math.sqrt(2.0),
            use_true_height=True,
            periodize=periodize,
        )
        tau = float(meta["rel"]) * math.sqrt(1.0 + r * r)
        eta = float(delta_aut) + tau
        row = dict(
            iter=it,
            r=float(r),
            rel=float(meta["rel"]),
            kappa_eq=float(meta["kappa_eq"]),
            delta_aut=float(delta_aut),
            tau_proxy=tau,
            eta_proxy=eta,
            scale=float(scale),
            per_generator={k: v["delta_max"] for k, v in per_gen.items()},
            reproject=bool(reproject and periodize),
        )
        history.append(row)
        print(
            f"  δ={delta_aut:.4e}  rel={meta['rel']:.4e}  η={eta:.4e}",
            flush=True,
        )
        # shrink search window as we go
        r_halfwidth *= 0.7

    best = min(history, key=lambda x: x["eta_proxy"])
    first = history[0]
    fac = first["eta_proxy"] / best["eta_proxy"] if best["eta_proxy"] > 0 else float("inf")
    out = dict(
        M=M,
        Y0=Y0,
        r0=r0,
        n_outer=n_outer,
        periodize=periodize,
        reproject=reproject,
        history=history,
        best=best,
        factor_vs_iter0=fac,
        seconds=time.time() - t0,
        language=(
            "Engineering Hejhal-style outer iterate: r-refine on collocation "
            "rel + two-cusp near-kernel + optional periodize reproject. "
            "Not certified Maass. Hard map unchanged."
        ),
    )
    print(
        f"\nBEST iter={best['iter']} r={best['r']:.6f}  "
        f"δ={best['delta_aut']:.4e} η={best['eta_proxy']:.4e}  "
        f"factor vs iter0={fac:.2f}×  ({out['seconds']:.1f}s)"
    )
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--M", type=int, default=32)
    p.add_argument("--r0", type=float, default=6.0)
    p.add_argument("--Y0", type=float, default=0.8)
    p.add_argument("--iters", type=int, default=4)
    p.add_argument("--r-halfwidth", type=float, default=0.75)
    p.add_argument("--periodize", action="store_true", default=True)
    p.add_argument("--no-periodize", action="store_true")
    p.add_argument("--no-reproject", action="store_true")
    p.add_argument("--n-face", type=int, default=16)
    p.add_argument("--json-out", type=str, default="")
    args = p.parse_args(argv)

    periodize = not args.no_periodize
    out = iterate(
        M=args.M,
        r0=args.r0,
        Y0=args.Y0,
        n_outer=args.iters,
        r_halfwidth=args.r_halfwidth,
        periodize=periodize,
        reproject=not args.no_reproject,
        n_face=args.n_face,
    )
    path = Path(args.json_out) if args.json_out else _HERE / "hejhal_iterate_result.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
