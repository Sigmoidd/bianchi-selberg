#!/usr/bin/env python3
"""
Larger Poincaré periodization + conforming V_h^{P1,per} residual diagnostics.

Path A — expand finite group for Poincaré sum:
  words of length ≤ L in {T1±, R, TiR, S} (default L=2), drop terms with
  pullback height < y_floor, re-measure PAIRINGS δ_aut.

Path B — conforming CR residual language:
  first positive eigenpair of V_h^{P1,per} (assemble_level_p):
    Rayleigh λ_h, J_h^{cross} ≡ 0 by construction, GLB sketch.
  This is the correct trial-space flavor for dual upper bounds
  (vs free Neumann lower-bound flavor).

Language: engineering. Not certified Maass / dual GREEN. Hard map unchanged.

Usage:
  python larger_poincare_residual.py
  python larger_poincare_residual.py --M 40 --word-len 2 --skip-fem
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

from delta_aut_pairing import (  # noqa: E402
    GEN,
    PERIOD_GROUP_H,
    _mat,
    mat_inv,
    near_kernel_multi,
    measure_jumps,
    eval_fourier,
    h3_act,
)
from cr_glb_p1_per import cr_glb, reference_h_max  # noqa: E402
from v_h_p1_per_spectrum import (  # noqa: E402
    assemble_periodic,
    float_neumann_eigs,
)
from scipy.sparse.linalg import eigsh  # noqa: E402


def mat_mul(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return A @ B


def expand_group(word_len: int = 2, include_S: bool = True) -> List[Tuple[str, np.ndarray]]:
    """
    Finite set of products of generators of length ≤ word_len.

    Generators: id, T1, T1inv, R, TiR, (+ S if include_S).
    """
    gens: List[Tuple[str, np.ndarray]] = [
        ("T1", GEN["T1"]),
        ("T1inv", _mat(1.0, -1.0, 0.0, 1.0)),
        ("R", GEN["R"]),
        ("TiR", GEN["TiR"]),
    ]
    if include_S:
        gens.append(("S", GEN["S"]))

    # start with identity
    group: Dict[str, np.ndarray] = {"id": _mat(1.0, 0.0, 0.0, 1.0)}
    frontier = [("id", group["id"])]
    for _ in range(word_len):
        new_front = []
        for name, mat in frontier:
            for gname, gmat in gens:
                prod = mat_mul(mat, gmat)
                key = f"{name}*{gname}" if name != "id" else gname
                # dedupe by matrix entries (rounded)
                sig = tuple(np.round(prod.real.ravel(), 10)) + tuple(
                    np.round(prod.imag.ravel(), 10)
                )
                # check existing
                found = False
                for ename, emat in group.items():
                    esig = tuple(np.round(emat.real.ravel(), 10)) + tuple(
                        np.round(emat.imag.ravel(), 10)
                    )
                    if esig == sig:
                        found = True
                        break
                if not found:
                    group[key] = prod
                    new_front.append((key, prod))
        frontier = new_front
    return list(group.items())


def eval_fourier_periodized_group(
    pts: np.ndarray,
    modes,
    coeffs: np.ndarray,
    r: float,
    theta: float,
    group: Sequence[Tuple[str, np.ndarray]],
    y_floor: float = 0.5,
) -> np.ndarray:
    """Poincaré sum over an arbitrary finite group list."""
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


def measure_jumps_group(
    coeffs,
    modes,
    r: float,
    Y0: float,
    theta: float,
    n_face: int,
    group: Sequence[Tuple[str, np.ndarray]],
    y_min: float = 1.0 / math.sqrt(2.0),
) -> Tuple[Dict[str, Any], float, float]:
    """Like measure_jumps but with custom Poincaré group for evaluation."""
    from delta_aut_pairing import face_samples, GEN as GENS

    g = max(4, int(math.ceil(math.sqrt(n_face))))
    xs = (np.arange(g) + 0.37) / g - 0.5
    X1, X2 = np.meshgrid(xs, xs, indexing="xy")
    ref_h = np.column_stack([X1.ravel(), X2.ravel(), np.full(X1.size, Y0)])
    f_ref = eval_fourier_periodized_group(
        ref_h, modes, coeffs, r, theta, group
    )
    scale = float(np.sqrt(np.mean(np.abs(f_ref) ** 2)))
    scale = max(scale, 1e-30)

    per_gen: Dict[str, Any] = {}
    deltas: List[float] = []
    for name, mat in GENS.items():
        pts = face_samples(name, n_face, Y0, y_min=y_min)
        f0 = eval_fourier_periodized_group(pts, modes, coeffs, r, theta, group)
        pts_g = np.zeros_like(pts)
        for j in range(pts.shape[0]):
            z = complex(pts[j, 0], pts[j, 1])
            y = float(pts[j, 2])
            zg, yg = h3_act(mat, z, y)
            pts_g[j, 0] = zg.real
            pts_g[j, 1] = zg.imag
            pts_g[j, 2] = max(yg, 1e-12)
        f1 = eval_fourier_periodized_group(pts_g, modes, coeffs, r, theta, group)
        jump = np.abs(f0 - f1)
        d_max = float(np.max(jump) / scale)
        d_rms = float(np.sqrt(np.mean(jump ** 2)) / scale)
        per_gen[name] = dict(
            n_samples=int(pts.shape[0]),
            delta_max=d_max,
            delta_rms=d_rms,
        )
        deltas.append(d_max)
    return per_gen, float(max(deltas)) if deltas else float("inf"), scale


def poincare_compare(
    M: int = 40,
    r: float = 6.0,
    Y0: float = 0.8,
    n_face: int = 16,
    jump_weight: float = 2.0,
) -> Dict[str, Any]:
    """Compare periodization group sizes: H-only vs words≤1 vs words≤2."""
    t0 = time.time()
    coeffs, modes, meta = near_kernel_multi(
        M, r, Y0, n_face=n_face, jump_weight=jump_weight
    )
    groups = {
        "H_only": PERIOD_GROUP_H,
        "words_leq_1": expand_group(1, include_S=True),
        "words_leq_2": expand_group(2, include_S=True),
    }
    results = {}
    for label, grp in groups.items():
        print(f"  periodize group={label}  |F|={len(grp)} ...", flush=True)
        per_gen, d_aut, scale = measure_jumps_group(
            coeffs, modes, r, Y0, 0.5, n_face, grp
        )
        results[label] = dict(
            group_size=len(grp),
            delta_aut=d_aut,
            scale=scale,
            per_generator=per_gen,
        )
        print(f"    δ_aut={d_aut:.4e}", flush=True)
    base = results["H_only"]["delta_aut"]
    for label, rec in results.items():
        rec["factor_vs_H"] = base / rec["delta_aut"] if rec["delta_aut"] > 0 else float("inf")
    return dict(
        M=M,
        r=r,
        Y0=Y0,
        collocation=meta,
        groups=results,
        seconds=time.time() - t0,
        language=(
            "Larger finite Poincaré periodization of hybrid multi near-kernel. "
            "Engineering residual, not certified Maass form."
        ),
    )


def conforming_fem_residual(
    N1: int = 4,
    N2: int = 2,
    N3: int = 2,
    neigs: int = 4,
) -> Dict[str, Any]:
    """
    First positive eigenmode on V_h^{P1,per}: Rayleigh + GLB + J_h language.

    Cross-copy jumps vanish by construction (periodic space).
    """
    t0 = time.time()
    sys.path.insert(0, str(_HERE.parent / "independent_exclusion"))
    from cr_prototype import KAPPA1 as K1  # noqa: E402

    data = assemble_periodic(N1, N2, N3, verbose=False)
    Q, M = data["Q"], data["M"]
    h_max = reference_h_max(N1, N2, N3)
    k = min(neigs, Q.shape[0] - 2)
    vals, vecs = eigsh(Q, k=k, M=M, sigma=0.0, which="LM", tol=1e-8, maxiter=400)
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]

    modes = []
    for i, lam in enumerate(vals):
        v = vecs[:, i]
        num = float(np.real(v.conj() @ (Q @ v)))
        den = float(np.real(v.conj() @ (M @ v)))
        ray = num / max(den, 1e-300)
        glb = cr_glb(float(lam), h_max, float(K1))
        modes.append(
            dict(
                index=i,
                lam_h=float(lam),
                rayleigh=ray,
                glb=glb,
                is_positive=float(lam) > 1e-6,
            )
        )
    first_pos = next((m for m in modes if m["is_positive"]), None)
    return dict(
        language=(
            "Conforming V_h^{P1,per} first-positive CR eigenmode: discrete "
            "Rayleigh is upper-flavor for the CR operator; GLB is lower-flavor "
            "under CR theory hypotheses. J_h^{cross}≡0 by construction. "
            "NOT certified continuum λ₁. Free Neumann is lower-bound flavor only."
        ),
        mesh=dict(N1=N1, N2=N2, N3=N3),
        n_dofs=data["n_dofs"],
        h_max=h_max,
        KAPPA1=float(K1),
        modes=modes,
        first_positive=first_pos,
        J_h_cross_by_construction=0.0,
        seconds=time.time() - t0,
        dual_language_pointer="dual_glb_language.md",
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--M", type=int, default=40)
    p.add_argument("--r", type=float, default=6.0)
    p.add_argument("--Y0", type=float, default=0.8)
    p.add_argument("--n-face", type=int, default=16)
    p.add_argument("--jump-weight", type=float, default=2.0)
    p.add_argument("--skip-poincare", action="store_true")
    p.add_argument("--skip-fem", action="store_true")
    p.add_argument("--N1", type=int, default=4)
    p.add_argument("--N2", type=int, default=2)
    p.add_argument("--N3", type=int, default=2)
    p.add_argument("--json-out", type=str, default="")
    args = p.parse_args(argv)

    out: Dict[str, Any] = dict(language="larger Poincaré + conforming residual")
    if not args.skip_poincare:
        print("=== Larger Poincaré periodization ===")
        out["poincare"] = poincare_compare(
            M=args.M,
            r=args.r,
            Y0=args.Y0,
            n_face=args.n_face,
            jump_weight=args.jump_weight,
        )
        print(
            f"({out['poincare']['seconds']:.1f}s)  "
            f"best δ among groups: "
            f"{min(g['delta_aut'] for g in out['poincare']['groups'].values()):.4e}"
        )

    if not args.skip_fem:
        print("\n=== Conforming V_h^{P1,per} residual language ===")
        out["conforming"] = conforming_fem_residual(
            N1=args.N1, N2=args.N2, N3=args.N3
        )
        fp = out["conforming"]["first_positive"]
        if fp:
            print(
                f"  first pos λ_h={fp['lam_h']:.6f}  "
                f"Rayleigh={fp['rayleigh']:.6f}  GLB={fp['glb']:.6f}"
            )
            print(f"  J_h_cross=0 (by construction)")
        print(f"  ({out['conforming']['seconds']:.2f}s)")

    path = Path(args.json_out) if args.json_out else _HERE / "larger_poincare_result.json"
    path.write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
