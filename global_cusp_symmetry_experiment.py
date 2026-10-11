#!/usr/bin/env python3
r"""EXPERIMENTAL: effect of the Q-parity sector decomposition on the Track-B physical Hejhal matrix.

This script assembles the repository's own physical rows (``assemble_physical_rows``) at the
frozen spectral parameter r = 6.7439020359331625 and compares singular-value structure of

    A            full 144 x n physical matrix (no symmetry),
    A B_{+}, A B_{-}   the same rows restricted to the even / odd sector (one column per C4-orbit).

Only midpoint (floating) SVDs are used.  Nothing here is a certificate: Arb is used to build
the entries, but ranks and singular values come from numpy.  Nothing produced here may be used
to flip a certification flag.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from flint import arb, ctx

import global_quarter_turn as gq
from track_b_two_cusp_hejhal import ValidatedWhittaker, assemble_physical_rows, midpoint_matrix

R_FROZEN = "6.7439020359331625"


def sector_matrix(cutoff: int, eps: int, n: int) -> np.ndarray:
    cols = gq.sector_basis(cutoff, eps)
    b = np.zeros((n, len(cols)))
    for j, col in enumerate(cols):
        for idx, v in col.items():
            b[idx, j] = v
    return b


def sv(m: np.ndarray) -> np.ndarray:
    return np.linalg.svd(m, compute_uv=False)


def numerical_rank(m: np.ndarray) -> int:
    tol = max(m.shape) * np.finfo(float).eps * max(np.linalg.norm(m, 2), 1.0)   # repo's _rank rule
    return int(np.linalg.matrix_rank(m, tol=tol))


def equilibrate_columns(m: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(m, axis=0)
    norms[norms == 0] = 1.0
    return m / norms


def analyse(cutoff: int, bits: int, points: int = 6) -> dict:
    ctx.prec = bits
    backend = ValidatedWhittaker(arb(R_FROZEN))
    t0 = time.time()
    rows, ledger, modes_inf, modes_zero = assemble_physical_rows(cutoff, points, backend, None)
    a = midpoint_matrix(rows)
    n = a.shape[1]
    norm_idx = next(i for i, m in enumerate(modes_inf) if m[:2] == (1, 0))
    nrm = np.zeros((1, n))
    nrm[0, norm_idx] = 1.0
    out = {"cutoff_M": cutoff, "rows": a.shape[0], "unknowns_full": n, "assemble_seconds": round(time.time() - t0, 1)}
    full = np.vstack([a, nrm])
    out["full"] = {
        "rank_repo_rule": numerical_rank(full),
        "rank_column_equilibrated_tol_1e-10": int(np.linalg.matrix_rank(equilibrate_columns(full), tol=1e-10)),
        "smallest_5_sv_equilibrated": [float(x) for x in sorted(sv(equilibrate_columns(full)))[:5]],
    }
    for eps, name in ((1, "even"), (-1, "odd")):
        b = sector_matrix(cutoff, eps, n)
        sec = np.vstack([a @ b, nrm @ b])
        s_eq = sv(equilibrate_columns(sec))
        out[name] = {
            "unknowns": b.shape[1],
            "rank_repo_rule": numerical_rank(sec),
            "rank_column_equilibrated_tol_1e-10": int(np.linalg.matrix_rank(equilibrate_columns(sec), tol=1e-10)),
            "smallest_5_sv_equilibrated": [float(x) for x in sorted(s_eq)[:5]],
            "rank_deficiency_repo_rule": b.shape[1] - numerical_rank(sec),
        }
    # normalized least squares  min ||A c||_2  s.t.  c_{(1,0)} = 1   (unprojected and per sector)
    def normalized_lstsq(mat: np.ndarray, norm_col: int):
        rest = [j for j in range(mat.shape[1]) if j != norm_col]
        rhs = -mat[:, norm_col]
        sol, *_ = np.linalg.lstsq(mat[:, rest], rhs, rcond=1e-13)
        c = np.zeros(mat.shape[1], dtype=complex)
        c[norm_col] = 1.0
        c[rest] = sol
        return c, float(np.linalg.norm(mat @ c))
    c_full, res_full = normalized_lstsq(a, norm_idx)
    keys = [('i', m[0], m[1]) for m in modes_inf] + [('z', m[0], m[1]) for m in modes_zero]
    solutions = {'unprojected': dict(zip(keys, c_full))}
    out["full"]["normalized_lstsq_residual_L2_rows"] = res_full
    sector_res = {}
    for eps, name in ((1, "even"), (-1, "odd")):
        b = sector_matrix(cutoff, eps, n)
        nc = int(np.nonzero(nrm @ b)[1][0])
        _c, res = normalized_lstsq(a @ b, nc)
        solutions[name] = dict(zip(keys, b @ _c))
        sector_res[name] = res
        out[name]["normalized_lstsq_residual_L2_rows"] = res
    # same-trial comparison: project the unprojected least-squares trial onto each parity
    # (U_Q c)_nu = c_{i nu};  P_eps = (I + eps U_Q)/2 on coefficient vectors
    perm = np.zeros((n, n))
    for offset, modes in ((0, modes_inf), (len(modes_inf), modes_zero)):
        index = {(m[0], m[1]): k for k, m in enumerate(modes)}
        for k, m in enumerate(modes):
            perm[offset + k, offset + index[gq.rot_i((m[0], m[1]))]] = 1.0
    same = {}
    for eps, name in ((1, "even"), (-1, "odd")):
        proj = 0.5 * (np.eye(n) + eps * perm)
        cp = proj @ c_full
        same[name] = {
            "coefficient_norm_fraction": float(np.linalg.norm(cp) / np.linalg.norm(c_full)),
            "residual_A_P_c": float(np.linalg.norm(a @ cp)),
            "residual_over_coefficient_norm": float(np.linalg.norm(a @ cp) / max(np.linalg.norm(cp), 1e-300)),
        }
    same["unprojected"] = {
        "residual_A_c": res_full,
        "residual_over_coefficient_norm": res_full / float(np.linalg.norm(c_full)),
    }
    out["same_trial_projection"] = same
    out["_solutions"] = solutions
    return out


def main() -> int:
    bits = 128
    cutoffs = [int(x) for x in sys.argv[1:]] or [1, 2, 3, 4]
    results = [analyse(m, bits) for m in cutoffs]
    drift = {}
    for name in ("unprojected", "even", "odd"):
        drift[name] = []
        for lo, hi in zip(results, results[1:]):
            small, big = lo["_solutions"][name], hi["_solutions"][name]
            common = [k for k in small if k in big]
            drift[name].append({"from_M": lo["cutoff_M"], "to_M": hi["cutoff_M"],
                                "l2_difference_on_common_modes": float(np.linalg.norm(
                                    np.array([small[k] - big[k] for k in common]))),
                                "l2_norm_of_smaller_cutoff_solution": float(np.linalg.norm(
                                    np.array([small[k] for k in common])))})
    for r in results:
        r.pop("_solutions")
    Path("global_cusp_symmetry_experiment_result.json").write_text(json.dumps(
        {"tier": "EXPERIMENTAL (floating SVD of Arb midpoints; not a certificate)",
         "r": R_FROZEN, "bits": bits, "results": results, "cutoff_drift": drift}, indent=2), encoding="utf-8")
    for r in results:
        print(json.dumps(r, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
