#!/usr/bin/env python3
"""
Engineering spectrum on V_h^{P1,per} = multi-copy CR with cross-copy pairings.

Reuses independent_exclusion.congruence_prototype.assemble_level_p, which
already builds the periodic (pairing-conforming) CR space for Γ₀(2+i):

  - cross-copy PAIRINGS faces identified (conforming)
  - self-identifications relaxed (Neumann)
  - top faces free for Neumann; zeroed for Dirichlet diagnostic

Language (honest)
-----------------
- Float Rayleigh quotients on V_h^{P1,per} are engineering diagnostics.
- Because the trial space is a *conforming* subspace of H¹(Γ₀\\H³) (up to
  CR nonconformity of face means — same as the exclusion pipeline), the
  discrete eigenvalues are *upper-bound flavor* for the CR discrete operator,
  NOT a certified dual upper bound for the true continuum first eigenvalue
  until CR GLB / Rump / residual path close.
- Free single-cell Neumann without pairings remains lower-bound flavor only.
- Does NOT flip hard map / rung4_certified.

Usage:
  python v_h_p1_per_spectrum.py
  python v_h_p1_per_spectrum.py --N1 4 --N2 2 --N3 2 --neigs 8
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh, splu

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
_IE = _ROOT / "independent_exclusion"
sys.path.insert(0, str(_IE))

from congruence_prototype import (  # noqa: E402
    INDEX,
    LEVEL,
    NP,
    assemble_level_p,
    set_level,
)


def top_dofs_from_traces(t_inf: np.ndarray, t_0: np.ndarray, tol: float = 1e-14) -> np.ndarray:
    """Periodic dofs with nonzero top-face mass (artificial boundary y=Y)."""
    t = np.abs(t_inf) + np.abs(t_0)
    return np.where(t > tol)[0]


def float_neumann_eigs(
    Q: csr_matrix,
    M: csr_matrix,
    neigs: int = 8,
) -> np.ndarray:
    """Smallest neigs generalized eigenvalues Qv=λMv (free top)."""
    k = min(neigs, Q.shape[0] - 2)
    # shift-invert around 0 for near-null + first positive modes
    try:
        vals = eigsh(
            Q,
            k=k,
            M=M,
            sigma=0.0,
            which="LM",
            tol=1e-8,
            maxiter=400,
            return_eigenvectors=False,
        )
    except Exception:
        # fallback: LOBPCG-style via dense if small
        if Q.shape[0] <= 800:
            from scipy.linalg import eigh

            vals = eigh(Q.toarray(), M.toarray(), eigvals_only=True)
            vals = vals[:k]
        else:
            raise
    return np.sort(np.real(vals))


def float_dirichlet_eigs(
    Q: csr_matrix,
    M: csr_matrix,
    top: np.ndarray,
    neigs: int = 8,
) -> np.ndarray:
    """Smallest neigs with Dirichlet zero on top-face dofs."""
    n = Q.shape[0]
    free = np.ones(n, dtype=bool)
    free[top] = False
    idx = np.where(free)[0]
    if idx.size < neigs + 2:
        raise RuntimeError("too few free dofs after Dirichlet")
    Qf = Q[idx][:, idx]
    Mf = M[idx][:, idx]
    k = min(neigs, idx.size - 2)
    try:
        vals = eigsh(
            Qf,
            k=k,
            M=Mf,
            sigma=0.0,
            which="LM",
            tol=1e-8,
            maxiter=400,
            return_eigenvectors=False,
        )
    except Exception:
        if idx.size <= 800:
            from scipy.linalg import eigh

            vals = eigh(Qf.toarray(), Mf.toarray(), eigvals_only=True)[:k]
        else:
            raise
    return np.sort(np.real(vals))


def assemble_periodic(
    N1: int = 4,
    N2: int = 2,
    N3: int = 2,
    verbose: bool = True,
) -> Dict[str, Any]:
    set_level("(2+i)")
    Q, M, avec, t_inf, t_0 = assemble_level_p(N1, N2, N3, verbose=verbose)
    top = top_dofs_from_traces(t_inf, t_0)
    one = np.ones(Q.shape[0])
    return dict(
        Q=Q,
        M=M,
        avec=avec,
        t_inf=t_inf,
        t_0=t_0,
        top_dofs=top,
        n_dofs=Q.shape[0],
        n_top=int(top.size),
        vol_check=float(one @ (M @ one)),
        Q1_norm=float(np.abs(Q @ one).max()),
        t_inf_sum=float(t_inf @ one),
        t_0_sum=float(t_0 @ one),
        N=(N1, N2, N3),
        level=LEVEL,
        NP=NP,
        INDEX=INDEX,
    )


def run_spectrum(
    N1: int = 4,
    N2: int = 2,
    N3: int = 2,
    neigs: int = 8,
    verbose: bool = True,
) -> Dict[str, Any]:
    t0 = time.time()
    data = assemble_periodic(N1, N2, N3, verbose=verbose)
    Q, M = data["Q"], data["M"]
    top = data["top_dofs"]

    if verbose:
        print(
            f"\n=== V_h^{{P1,per}} spectrum (engineering) ===\n"
            f"level={data['level']} NP={data['NP']} copies={data['INDEX']}\n"
            f"mesh {N1}x{N2}x{N3}  dofs={data['n_dofs']}  top={data['n_top']}\n"
            f"1'M1={data['vol_check']:.6f}  |Q1|={data['Q1_norm']:.2e}  "
            f"t∞(1)={data['t_inf_sum']:.4f}  t0(1)={data['t_0_sum']:.4f}"
        )

    lam_N = float_neumann_eigs(Q, M, neigs=neigs)
    lam_D = float_dirichlet_eigs(Q, M, top, neigs=neigs)
    brackets = [
        (k + 1, float(lam_N[k]), float(lam_D[k]))
        for k in range(min(len(lam_N), len(lam_D)))
    ]

    # first positive Neumann mode (skip near-null constant)
    pos_N = [float(x) for x in lam_N if x > 1e-6]
    first_pos = pos_N[0] if pos_N else float("nan")

    out = dict(
        language=(
            "Float CR spectrum on multi-copy pairing-conforming space "
            "V_h^{P1,per} for Γ₀(2+i). Engineering diagnostic only — not "
            "certified dual upper bound; not free-Neumann lower-bound flavor. "
            "Hard map unchanged."
        ),
        level=data["level"],
        NP=data["NP"],
        INDEX=data["INDEX"],
        mesh=dict(N1=N1, N2=N2, N3=N3),
        n_dofs=data["n_dofs"],
        n_top=data["n_top"],
        vol_check=data["vol_check"],
        Q1_norm=data["Q1_norm"],
        t_inf_sum=data["t_inf_sum"],
        t_0_sum=data["t_0_sum"],
        lam_N=[float(x) for x in lam_N],
        lam_D=[float(x) for x in lam_D],
        brackets=brackets,
        first_positive_Neumann=first_pos,
        seconds=time.time() - t0,
        warnings=[
            "CR is nonconforming; continuum upper bound needs CR theory / GLB.",
            "Do not claim certified λ1 or dual GREEN from this spectrum alone.",
        ],
    )
    if verbose:
        print(f"\n{'k':>4}  {'λ^N':>12}  {'λ^D':>12}")
        for k, lo, hi in brackets:
            print(f"{k:4d}  {lo:12.6f}  {hi:12.6f}")
        print(f"\nfirst positive Neumann ≈ {first_pos:.6f}")
        print(f"({out['seconds']:.2f}s)")
        for w in out["warnings"]:
            print(f"  WARN: {w}")
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--N1", type=int, default=4)
    p.add_argument("--N2", type=int, default=2)
    p.add_argument("--N3", type=int, default=2)
    p.add_argument("--neigs", type=int, default=8)
    p.add_argument("--json-out", type=str, default="")
    args = p.parse_args(argv)

    out = run_spectrum(args.N1, args.N2, args.N3, neigs=args.neigs)
    path = Path(args.json_out) if args.json_out else _HERE / "v_h_p1_per_spectrum.json"
    # drop non-serializable
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {path}")
    # infrastructure checks
    ok = (
        out["Q1_norm"] < 1e-6
        and abs(out["lam_N"][0]) < 1e-4
        and out["brackets"][0][1] <= out["brackets"][0][2] + 1e-8
    )
    print("PASS infrastructure" if ok else "FAIL infrastructure checks")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
