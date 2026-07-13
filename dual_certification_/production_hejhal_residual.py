#!/usr/bin/env python3
"""
Production-shaped two-cusp Hejhal residual (Γ₀(2+i) model).

Builds the two-cusp coupled collocation operator (Rung 3), extracts the
near-kernel coefficient vector at a candidate r, and measures:

  1. Discrete residual: rel = σ_min/σ_max after amp + equilibration
  2. PAIRINGS face automorphy δ_aut via delta_aut_pairing (true H³ action)
  3. Optional: residual with diagonal reg stripped from reporting

Language
--------
- More production-shaped than single-cusp collocation alone (two cusps + S).
- Still a *model* Fourier/collocation residual, not a certified Maass form.
- Does not flip hard map.

Usage:
  python production_hejhal_residual.py
  python production_hejhal_residual.py --M 48 --r 6.0 --Y0 0.8
  python production_hejhal_residual.py --M 32 --r-scan 5.5:7.0:7
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
    equilibrate_real,
    complex_to_real_mid_rad,
)
from hejhal_conditioning import equilibrate  # noqa: E402
from delta_aut_pairing import (  # noqa: E402
    measure_jumps,
    gaussian_modes,
    eval_fourier_periodized,
    eval_fourier,
)


def two_cusp_near_kernel(
    M: int,
    r: float,
    Y0: float,
    theta: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Right singular vector of equilibrated two-cusp operator.

    Returns (a_inf, a_0, meta) complex coefficient vectors (unit joint ℓ²).
    """
    sysb = build_block_system(M, float(r), Y0, theta=theta)
    mid, rad = coupled_mid_rad(sysb)
    mid_p, rad_p, w = precondition_right(mid, rad, sysb.w_inf, sysb.w_0)
    # column normalize further
    col = np.maximum(np.linalg.norm(mid_p, axis=0), 1e-300)
    G = mid_p / col[None, :]
    Geq = equilibrate(G, n_iter=6)
    _u, s, vh = np.linalg.svd(Geq, full_matrices=False)
    sig_min = float(s[-1])
    sig_max = float(s[0]) if len(s) else 1.0
    u = vh[-1, :].conj()
    # map back through col + amp
    a_full = u / col
    # amp was in precondition_right: mid_p = mid * D^{-1}, so physical a = D^{-1} u_eq roughly
    # w = concat(w_inf, w_0); physical trial coeffs ~ a_full / || ||
    a_full = a_full / max(np.linalg.norm(a_full), 1e-300)
    n = sysb.n
    a_inf = a_full[:n].astype(np.complex128)
    a_0 = a_full[n:].astype(np.complex128)
    # residual of un-equilibrated preconditioned system
    ra = float(np.linalg.norm(mid_p @ (a_full * col)))  # noisy; use rel
    meta = dict(
        M=M,
        r=r,
        Y0=Y0,
        theta=theta,
        n_modes_per_cusp=n,
        sigma_min=sig_min,
        sigma_max=sig_max,
        rel=sig_min / max(sig_max, 1e-300),
        kappa_eq=sig_max / max(sig_min, 1e-300),
        reg=float(sysb.meta.get("reg", 0.0)),
        language=(
            "Two-cusp collocation near-kernel (model S + height-matched "
            "automorphy blocks). reg on V diagonals is present in assembly."
        ),
    )
    return a_inf, a_0, meta


def residual_at_r(
    M: int,
    r: float,
    Y0: float = 0.8,
    theta: float = 0.5,
    n_face: int = 16,
    periodize: bool = True,
) -> Dict[str, Any]:
    t0 = time.time()
    a_inf, a_0, meta = two_cusp_near_kernel(M, r, Y0, theta=theta)
    modes = gaussian_modes(M)
    # measure PAIRINGS jumps on cusp-∞ Fourier trial (primary)
    per_gen, delta_aut, scale = measure_jumps(
        a_inf,
        modes,
        r,
        Y0,
        theta,
        n_face,
        y_min=1.0 / math.sqrt(2.0),
        use_true_height=True,
        periodize=periodize,
    )
    # also report cusp-0 jumps
    per_gen0, delta0, _ = measure_jumps(
        a_0,
        modes,
        r,
        Y0,
        theta,
        n_face,
        y_min=1.0 / math.sqrt(2.0),
        use_true_height=True,
        periodize=periodize,
    )
    tau_proxy = float(meta["rel"]) * math.sqrt(1.0 + r * r)
    eta = float(delta_aut) + tau_proxy
    return dict(
        M=M,
        r=r,
        Y0=Y0,
        periodize=periodize,
        two_cusp=meta,
        cusp_inf=dict(
            delta_aut=delta_aut,
            scale=scale,
            per_generator=per_gen,
            a_norm=float(np.linalg.norm(a_inf)),
        ),
        cusp_0=dict(
            delta_aut=delta0,
            per_generator=per_gen0,
            a_norm=float(np.linalg.norm(a_0)),
        ),
        tau_proxy=tau_proxy,
        eta_proxy=eta,
        seconds=time.time() - t0,
        language=(
            "Production-shaped residual: two-cusp near-kernel + true H³ "
            "PAIRINGS δ_aut. Model Fourier, not certified Maass. "
            "Hard map unchanged."
        ),
    )


def r_scan(
    M: int,
    r_min: float,
    r_max: float,
    n_r: int,
    Y0: float = 0.8,
    periodize: bool = True,
) -> Dict[str, Any]:
    rs = np.linspace(r_min, r_max, n_r)
    rows = []
    best = None
    for r in rs:
        print(f"  r={r:.4f} ...", flush=True)
        rec = residual_at_r(M, float(r), Y0=Y0, periodize=periodize)
        rows.append(
            dict(
                r=float(r),
                delta_aut=rec["cusp_inf"]["delta_aut"],
                delta_0=rec["cusp_0"]["delta_aut"],
                rel=rec["two_cusp"]["rel"],
                tau_proxy=rec["tau_proxy"],
                eta_proxy=rec["eta_proxy"],
                kappa_eq=rec["two_cusp"]["kappa_eq"],
                seconds=rec["seconds"],
            )
        )
        print(
            f"    δ_inf={rec['cusp_inf']['delta_aut']:.3e}  "
            f"rel={rec['two_cusp']['rel']:.3e}  η={rec['eta_proxy']:.3e}",
            flush=True,
        )
        if best is None or rec["eta_proxy"] < best["eta_proxy"]:
            best = rec
    return dict(rows=rows, best=best, r_grid=[float(x) for x in rs])


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--M", type=int, default=48)
    p.add_argument("--r", type=float, default=6.0)
    p.add_argument("--Y0", type=float, default=0.8)
    p.add_argument("--n-face", type=int, default=16)
    p.add_argument("--no-periodize", action="store_true")
    p.add_argument(
        "--r-scan",
        type=str,
        default="",
        help="optional r_min:r_max:n  e.g. 5.5:7.0:7",
    )
    p.add_argument("--json-out", type=str, default="")
    args = p.parse_args(argv)

    print("=== Production-shaped two-cusp Hejhal residual ===")
    periodize = not args.no_periodize
    if args.r_scan:
        parts = args.r_scan.split(":")
        r_min, r_max, n_r = float(parts[0]), float(parts[1]), int(parts[2])
        print(f"M={args.M} Y0={args.Y0} scan r∈[{r_min},{r_max}] n={n_r}")
        out = r_scan(args.M, r_min, r_max, n_r, Y0=args.Y0, periodize=periodize)
        best = out["best"]
        print(
            f"\nbest r≈{best['r']:.4f}  δ_inf={best['cusp_inf']['delta_aut']:.4e}  "
            f"η={best['eta_proxy']:.4e}"
        )
    else:
        print(f"M={args.M} r={args.r} Y0={args.Y0} periodize={periodize}")
        out = residual_at_r(
            args.M,
            args.r,
            Y0=args.Y0,
            n_face=args.n_face,
            periodize=periodize,
        )
        print(f"  two-cusp rel={out['two_cusp']['rel']:.4e}  κ={out['two_cusp']['kappa_eq']:.1f}")
        print(f"  δ_aut(∞)={out['cusp_inf']['delta_aut']:.6e}")
        print(f"  δ_aut(0)={out['cusp_0']['delta_aut']:.6e}")
        print(f"  τ_proxy={out['tau_proxy']:.6e}  η≈{out['eta_proxy']:.6e}")
        print("  per-gen (∞):")
        for g, rec in out["cusp_inf"]["per_generator"].items():
            print(f"    {g:4s}  δ_max={rec['delta_max']:.4e}")
        print(f"  ({out['seconds']:.2f}s)")

    path = Path(args.json_out) if args.json_out else _HERE / "production_hejhal_residual.json"
    path.write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
