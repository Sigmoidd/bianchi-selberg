#!/usr/bin/env python3
"""
Residual scan table: M ↦ (σ_min, rel, τ_disc, κ_eq, n_pts, n_modes).

Morningbrief priority 2 / path (3): quantitative residual vs truncation M
with honest τ_disc := ||V a||₂ (||a||₂=1 near-kernel) and equilibrated κ.

Uses single-cusp collocation (hejhal_conditioning) — same model residual
family as Rung 4 scan, not production Maass residual.

Default M list: 100,200,400,800 (add 1000,1200 via --M).
n_pts = g² with g = 2⌊√M⌋+2 (Nyquist).

Usage:
  python residual_scan_M.py
  python residual_scan_M.py --M 100,200,400 --r 6.0
  python residual_scan_M.py --M 800,1000,1200 --r 6.0 --Y0 0.8
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from hejhal_conditioning import (  # noqa: E402
    build_single_cusp_V,
    equilibrate,
    fit_loglog,
)


def residual_at_M(
    M: int,
    r: float = 6.0,
    Y0: float = 0.8,
    theta: float = 0.5,
) -> Dict[str, Any]:
    t0 = time.time()
    V, amps, modes = build_single_cusp_V(M, Y0=Y0, r=r, theta=theta)
    n_pts, n_modes = V.shape
    F = V / np.maximum(amps[None, :], 1e-300)
    Feq = equilibrate(F, n_iter=6)
    s = np.linalg.svd(Feq, compute_uv=False)
    sig_min = float(s[-1])
    sig_max = float(s[0]) if len(s) else 1.0
    rel = sig_min / max(sig_max, 1e-300)
    # near-kernel in equilibrated coords → physical coeffs
    _u, _s, vh = np.linalg.svd(Feq, full_matrices=False)
    u = vh[-1, :].conj()
    a = u / np.maximum(amps, 1e-300)
    a = a / max(np.linalg.norm(a), 1e-300)
    # ||V a|| for exact SVD right vector is ~0 after amp undo; report both:
    #   tau_svd  = σ_min of equilibrated system (scale-free residual size)
    #   tau_disc = ||V a||₂ (often machine-tiny; kept for API honesty)
    tau_disc = float(np.linalg.norm(V @ a))
    tau_svd = sig_min
    delta_proxy = rel
    tau_proxy = rel * math.sqrt(1.0 + r * r)
    eta_proxy = delta_proxy + tau_proxy
    return dict(
        M=M,
        r=r,
        Y0=Y0,
        n_pts=n_pts,
        n_modes=n_modes,
        sigma_min=sig_min,
        sigma_max=sig_max,
        rel=rel,
        kappa_eq=sig_max / max(sig_min, 1e-300),
        tau_disc=tau_disc,
        tau_svd=tau_svd,
        delta_proxy=delta_proxy,
        tau_proxy=tau_proxy,
        eta_proxy=eta_proxy,
        seconds=time.time() - t0,
    )


def run_scan(
    Ms: Sequence[int],
    r: float = 6.0,
    Y0: float = 0.8,
    theta: float = 0.5,
) -> Dict[str, Any]:
    rows: List[Dict[str, Any]] = []
    for M in Ms:
        print(f"  M={M} ...", flush=True)
        rec = residual_at_M(M, r=r, Y0=Y0, theta=theta)
        rows.append(rec)
        print(
            f"    n_pts={rec['n_pts']} n_modes={rec['n_modes']}  "
            f"rel={rec['rel']:.3e} τ_svd={rec['tau_svd']:.3e}  "
            f"κ_eq={rec['kappa_eq']:.1f}  ({rec['seconds']:.1f}s)",
            flush=True,
        )
    Ms_f = [row["M"] for row in rows]
    rels = [row["rel"] for row in rows]
    taus = [row["tau_disc"] for row in rows]
    _, b_rel = fit_loglog(Ms_f, rels)
    _, b_tau = fit_loglog(Ms_f, taus)
    return dict(
        r=r,
        Y0=Y0,
        theta=theta,
        rows=rows,
        b_loglog_rel=b_rel,
        b_loglog_tau_disc=b_tau,
        language=(
            "Model collocation residual (single-cusp Hejhal-like). "
            "τ_disc = ||V a||₂ with unit-ℓ² near-kernel; rel = σ_min/σ_max "
            "after amp+Sinkhorn equilibration. Not Maass residual."
        ),
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--M",
        type=str,
        default="100,200,400,800",
        help="comma-separated M list (default 100,200,400,800)",
    )
    p.add_argument("--r", type=float, default=6.0)
    p.add_argument("--Y0", type=float, default=0.8)
    p.add_argument("--theta", type=float, default=0.5)
    p.add_argument("--csv-out", type=str, default="")
    args = p.parse_args(argv)

    Ms = [int(x.strip()) for x in args.M.split(",") if x.strip()]
    print("=== Residual scan vs M (path 3 diagnostic) ===")
    print(f"r={args.r}  Y0={args.Y0}  M={Ms}")
    out = run_scan(Ms, r=args.r, Y0=args.Y0, theta=args.theta)
    print(f"\n  log-log slope b(rel)≈{out['b_loglog_rel']:.3f}  "
          f"b(τ_disc)≈{out['b_loglog_tau_disc']:.3f}")
    print(f"  {out['language']}")

    csv_path = Path(args.csv_out) if args.csv_out else _HERE / "residual_scan_M.csv"
    fields = [
        "M", "r", "Y0", "n_pts", "n_modes", "sigma_min", "rel",
        "kappa_eq", "tau_svd", "tau_disc", "delta_proxy", "tau_proxy",
        "eta_proxy", "seconds",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in out["rows"]:
            w.writerow(row)
    print(f"  wrote {csv_path}")

    # short markdown summary
    md = _HERE / "residual_scan_M.md"
    lines = [
        "# Residual scan vs M",
        "",
        f"r={args.r}, Y0={args.Y0}, theta={args.theta}",
        "",
        "| M | n_pts | n_modes | rel | τ_svd | κ_eq | η_proxy | s |",
        "|--:|------:|--------:|----:|------:|-----:|--------:|--:|",
    ]
    for row in out["rows"]:
        lines.append(
            f"| {row['M']} | {row['n_pts']} | {row['n_modes']} | "
            f"{row['rel']:.3e} | {row['tau_svd']:.3e} | {row['kappa_eq']:.1f} | "
            f"{row['eta_proxy']:.3e} | {row['seconds']:.1f} |"
        )
    lines += [
        "",
        f"log-log slopes: b(rel)≈{out['b_loglog_rel']:.3f}, "
        f"b(τ_disc)≈{out['b_loglog_tau_disc']:.3f}",
        "",
        out["language"],
        "",
        "Non-certifying. Hard map unchanged.",
    ]
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  wrote {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
