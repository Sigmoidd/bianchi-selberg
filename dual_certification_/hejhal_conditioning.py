#!/usr/bin/env python3
"""
Milestone 2 / Rung 3 — Hejhal collocation condition-number diagnostic.

Builds a *model* single-cusp (and stub two-cusp) Hejhal-like collocation /
moment matrix V(M) at fixed r ≈ 6.622, Y0 = 0.8:

  - Fourier modes β = a+bi ∈ Z[i], 0 < |β|² ≤ M
  - Collocation on the torus C / Z[i] at height y = Y0
  - Entries (single cusp):
        V_{j,β} = K_amp(2π|β|Y0) · exp(2π i (a x_j + b y_j))
    with K_amp a float model of |K_{ir}(·)| (scipy kve when available,
    else the elementary majorant √(π/(2y)) e^{-y}).

  - Diagonal (right) preconditioner D = diag(K_amp(β)), studying
        κ(V D^{-1})   and equivalently column-normalised V
  - Optional two-cusp block stub
        [[V, -S], [-S^*, V]]
    with block-diagonal preconditioner.

Float diagnostic only — NON-CERTIFYING (AgentReady §6 Rung 3 / §8 Milestone 2).
Does NOT claim Rung 3 complete (Krawczyk for N=5 is out of scope here).

Outputs (next to this script):
  - hejhal_kappa_vs_M.csv
  - hejhal_kappa_vs_M.png   (if matplotlib available)
  - prints log-log slopes b for raw vs preconditioned κ

Usage:
  python hejhal_conditioning.py
  python hejhal_conditioning.py --M 100,200,400,800 --two-cusp
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Optional scipy K-Bessel (scaled)
# ---------------------------------------------------------------------------
_HAS_KVE = False
try:
    from scipy.special import kve as _kve  # exp(z) K_v(z)

    _HAS_KVE = True
except Exception:  # pragma: no cover
    _kve = None  # type: ignore

# Local majorant helper (optional; script is self-contained if import fails)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from lemma_K import C_K, ball_mid, r2_exact_gaussian  # noqa: E402

    _HAS_LEMMA_K = True
except Exception:  # pragma: no cover
    _HAS_LEMMA_K = False

    def r2_exact_gaussian(n: int) -> int:  # type: ignore
        """Count representations n = a²+b² with sign/order (0 if none)."""
        if n < 0:
            return 0
        c = 0
        s = int(math.isqrt(n))
        for a in range(-s, s + 1):
            b2 = n - a * a
            if b2 < 0:
                continue
            b = int(math.isqrt(b2))
            if b * b == b2:
                c += 1 if b == 0 else 2
        return c


# ---------------------------------------------------------------------------
# Constants (Then first eigenvalue scale)
# ---------------------------------------------------------------------------
R_DEFAULT = 6.62212
Y0_DEFAULT = 0.8
THETA_DEFAULT = 0.5


def gaussian_modes(M: int) -> List[Tuple[int, int, int]]:
    """
    All β = a+bi ∈ Z[i] with 0 < N(β) = a²+b² ≤ M.
    Returns list of (a, b, n).
    """
    modes: List[Tuple[int, int, int]] = []
    s = int(math.isqrt(M))
    for a in range(-s, s + 1):
        for b in range(-s, s + 1):
            n = a * a + b * b
            if 0 < n <= M:
                modes.append((a, b, n))
    # Deterministic order: by norm, then a, then b
    modes.sort(key=lambda t: (t[2], t[0], t[1]))
    return modes


def k_bessel_amp(y: float, r: float) -> float:
    """
    Float model of |K_{ir}(y)| for diagnostic use only.

    Prefer scipy.special.kve(1j*r, y) * exp(-y) when available and finite;
    fall back to elementary majorant √(π/(2y)) e^{-y} (C_K = 1).
    """
    if y <= 0.0:
        return 1.0
    majorant = math.sqrt(math.pi / (2.0 * y)) * math.exp(-y)
    if _HAS_KVE:
        try:
            # kve(v,z) = exp(z) K_v(z); order is pure imaginary ir
            val = _kve(1j * r, y)
            amp = abs(complex(val)) * math.exp(-y)
            if math.isfinite(amp) and amp > 0.0:
                # Guard: never exceed the elementary majorant by a huge factor
                # (kve can be noisy at extreme arguments).
                return min(amp, majorant * 10.0) if majorant > 0 else amp
        except Exception:
            pass
    return max(majorant, 1e-300)


def collocation_grid(g: int, shift: float = 0.5) -> np.ndarray:
    """Full g×g grid on the unit torus [0,1)² (C / Z[i] model)."""
    xs = (np.arange(g, dtype=np.float64) + shift) / g
    X, Y = np.meshgrid(xs, xs, indexing="xy")
    pts = np.column_stack([X.ravel(), Y.ravel()])
    return pts


def inversion_h3(x1: np.ndarray, x2: np.ndarray, y: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Model cusp-pairing isometry on H³: inversion in the unit hemisphere
        (x, y) ↦ (x / R², y / R²),  R² = |x|² + y².
    Stand-in for a Bianchi generator that swaps / mixes cusps.
    """
    R2 = x1 * x1 + x2 * x2 + y * y
    return x1 / R2, x2 / R2, y / R2


def mode_weights(
    norms: np.ndarray,
    y_height: float,
    r: float,
    theta: float,
    weight_hecke: bool,
) -> np.ndarray:
    """w_β(y) = n^θ · |K_ir(2π √n y)|  (float model), fixed height."""
    return mode_weights_many_heights(
        norms,
        np.full(1, y_height, dtype=np.float64),
        r,
        theta,
        weight_hecke,
    )[0]


def mode_weights_many_heights(
    norms: np.ndarray,
    heights: np.ndarray,
    r: float,
    theta: float,
    weight_hecke: bool,
) -> np.ndarray:
    """
    W[j,k] = n_k^θ · |K_ir(2π √n_k · y_j)|  using the elementary majorant
    √(π/(2y)) e^{-y} (float diagnostic; matches K-decay that drives κ).

    Using the majorant (not per-entry scipy.kve) keeps the O(n_pts·n_modes)
    fill fast and numerically consistent across heights.
    """
    # y_arg[j,k] = 2π √n_k y_j
    sqrtn = np.sqrt(norms.astype(np.float64))  # (n_modes,)
    y_arg = 2.0 * math.pi * heights.astype(np.float64)[:, None] * sqrtn[None, :]
    # majorant |K| ≤ √(π/(2y)) e^{-y}
    y_safe = np.maximum(y_arg, 1e-300)
    amp = np.sqrt(math.pi / (2.0 * y_safe)) * np.exp(-y_safe)
    if weight_hecke:
        amp = amp * np.power(norms.astype(np.float64), theta)[None, :]
    return np.maximum(amp, 1e-300)


def build_single_cusp_V(
    M: int,
    Y0: float = Y0_DEFAULT,
    r: float = R_DEFAULT,
    theta: float = THETA_DEFAULT,
    weight_hecke: bool = True,
    height_mismatch: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, List[Tuple[int, int, int]]]:
    """
    Model single-cusp Hejhal-like collocation matrix V (complex128).

    At sample points x_j on the horosphere y = Y0 we enforce a model
    automorphy residual under the inversion isometry γ (H³ unit hemisphere):

        V[j, β] = w_β(Y0) · ( e^{2π i β·x_j} − e^{2π i β·x_γ(x_j)} )
                  − height_mismatch · ΔW_{jβ} e^{2π i β·x_γ}

    with w_β(y) = N(β)^θ |K_ir(2π|β|y)|.

    Default height_mismatch=0: same-height phase automorphy (columns are pure
    mode-amplitude times a unitary Fourier difference). This is the natural
    model for the Rung 3 claim that a diagonal D=diag(w_β) yields polynomial
    κ growth — matching AgentReady's "naïve ~M^76 from K-decay; after D, b<4".

    height_mismatch∈(0,1]: blends in true pull-back weights w_β(y_γ(x_j)),
    a stress test (height-mismatched inversion) needing stronger precond / Y0.

    Grid: g = 2⌈√M⌉+2 (Nyquist), V is (g² × n_modes).

    Returns (V, amps, modes) with amps = w_β(Y0).
    """
    modes = gaussian_modes(M)
    n_modes = len(modes)
    if n_modes == 0:
        raise ValueError(f"no Gaussian modes for M={M}")
    s = int(math.isqrt(M))
    g = 2 * s + 2
    # Non-centered shift: avoid accidental exact DFT orthogonality
    pts = collocation_grid(g, shift=0.37)
    n_pts = pts.shape[0]
    a = np.array([m[0] for m in modes], dtype=np.float64)
    bcoef = np.array([m[1] for m in modes], dtype=np.float64)
    norms = np.array([m[2] for m in modes], dtype=np.float64)

    amps = mode_weights(norms, Y0, r, theta, weight_hecke)

    x1 = pts[:, 0].copy()
    x2 = pts[:, 1].copy()
    x1c = x1 - np.round(x1)
    x2c = x2 - np.round(x2)
    gx1, gx2, gy = inversion_h3(x1c, x2c, Y0)

    # Free Fourier at height Y0
    phase0 = 2.0 * math.pi * (
        np.outer(x1c, a) + np.outer(x2c, bcoef)
    )
    # Pull-back phases under γ
    phase1 = 2.0 * math.pi * (
        np.outer(gx1, a) + np.outer(gx2, bcoef)
    )
    F0 = np.exp(1j * phase0)
    F1 = np.exp(1j * phase1)

    # Same-height automorphy residual, column-scaled by w_β(Y0)
    V = ((F0 - F1) * amps[np.newaxis, :]) / math.sqrt(float(n_pts))

    if height_mismatch > 0.0:
        # Stress term: true height-dependent pull-back weights
        W1 = mode_weights_many_heights(norms, gy, r, theta, weight_hecke)
        delta = (W1 - amps[np.newaxis, :]) * F1
        V = V - (height_mismatch * delta) / math.sqrt(float(n_pts))

    return V, amps, modes

def build_two_cusp_block(
    V: np.ndarray,
    amps: np.ndarray,
    coupling: float = 0.1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Stub two-cusp coupled collocation (tall block form)

        B [a^∞; a^0] = [ V a^∞ - c V a^0 ;
                         -c V a^∞ + V a^0 ]

    Block-diagonal preconditioner D₂ = D ⊕ D with D = diag(amps).
    """
    c = coupling
    top = np.concatenate([V, -c * V], axis=1)
    bot = np.concatenate([-c * V, V], axis=1)
    V2 = np.concatenate([top, bot], axis=0)
    amps2 = np.concatenate([amps, amps])
    return V2, amps2


def two_cusp_precond_kappa(
    V: np.ndarray,
    amps: np.ndarray,
    coupling: float = 0.1,
) -> float:
    """
    κ of the *preconditioned* two-cusp stub (block-column form):

        B = [ U | -c U_ext ]
            [ U |  c U_ext ]   (two stacked observation blocks)

    For tall U (n_pts × n_modes) we stack two copies of the single-cusp
    skeleton with off-diagonal coupling c, modelling two-cusp sampling
    of (a^∞, a^0).  Equivalent square block form is used when U is square.
    """
    U = V / amps[np.newaxis, :]
    n_pts, n_modes = U.shape
    c = coupling
    # Block columns for unknowns (a^∞, a^0): observations at both cusps
    #   cusp ∞:  U a^∞  - c U a^0
    #   cusp 0: -c U a^∞ +   U a^0
    top = np.concatenate([U, -c * U], axis=1)
    bot = np.concatenate([-c * U, U], axis=1)
    B = np.concatenate([top, bot], axis=0)
    return cond_2(B)

def cond_2(A: np.ndarray) -> float:
    """2-norm condition number via SVD singular values (float diagnostic)."""
    if A.size == 0:
        return float("nan")
    try:
        s = np.linalg.svd(A, compute_uv=False)
    except np.linalg.LinAlgError:
        return float("inf")
    s = np.asarray(s, dtype=np.float64)
    s = s[np.isfinite(s) & (s > 0)]
    if s.size == 0:
        return float("inf")
    return float(s[0] / s[-1])


def equilibrate(V: np.ndarray, n_iter: int = 5) -> np.ndarray:
    """
    Sinkhorn-style row/column ∞-norm equilibration (diagonal preconditioning).
    Returns D_row^{-1} V D_col^{-1} with unit peak row/column scale.
    """
    A = np.array(V, dtype=np.complex128, copy=True)
    for _ in range(n_iter):
        row_sc = np.max(np.abs(A), axis=1)
        row_sc = np.maximum(row_sc, 1e-300)
        A = A / row_sc[:, None]
        col_sc = np.max(np.abs(A), axis=0)
        col_sc = np.maximum(col_sc, 1e-300)
        A = A / col_sc[None, :]
    return A


def kappa_raw_and_precond(
    V: np.ndarray,
    amps: np.ndarray,
) -> Dict[str, float]:
    """
    Estimate κ(V) and several diagonally preconditioned variants.

    Raw κ: for extreme dynamic range, float SVD on V is unreliable (columns
    underflow relative to eps·||V||). We report
        κ_raw_stable = (max amps / min amps) * κ(V D_amp^{-1})
    as the scale-aware diagnostic.

    Preconditioners (all diagonal / pair of diagonals):
      1. Right mode-amplitude:     V D_amp^{-1},  D_amp = diag(w_β(Y0))
      2. Column Jacobi:            V D_col^{-1},  D_col = diag(||V e_β||_2)
      3. Row+col equilibration:    D_row^{-1} V D_col^{-1}  (Sinkhorn)
      4. Amp-right then row Jacobi on that matrix (hybrid)

    AgentReady Rung 3 stopping condition uses the best diagonal scheme's
    log-log slope b; we label each clearly in the CSV/report.
    """
    amps = np.asarray(amps, dtype=np.float64)
    ratio = float(amps.max() / max(amps.min(), 1e-300))

    # Right preconditioning by mode amplitudes: V D_amp^{-1}
    F = V / amps[np.newaxis, :]
    kappa_right = cond_2(F)

    # Column Jacobi on raw V (using stable column norms via F)
    # ||V e_β|| = amps_β ||F e_β||
    F_col_norms = np.linalg.norm(F, axis=0)
    F_col_norms = np.maximum(F_col_norms, 1e-300)
    V_col = F / F_col_norms[np.newaxis, :]  # = V / ||V e_β||
    kappa_col = cond_2(V_col)

    # Row Jacobi on amp-preconditioned matrix (left diagonal after right amp)
    row_n = np.linalg.norm(F, axis=1)
    row_n = np.maximum(row_n, 1e-300)
    F_row = F / row_n[:, None]
    kappa_amp_row = cond_2(F_row)

    # Full row+column equilibration on amp-preconditioned matrix
    F_eq = equilibrate(F, n_iter=6)
    kappa_eq = cond_2(F_eq)

    # Direct raw cond when dynamic range permits
    if ratio < 1e12:
        kappa_direct = cond_2(V)
    else:
        kappa_direct = float("inf")

    kappa_raw_stable = ratio * kappa_right

    # Primary preconditioned κ for Rung 3: equilibration (best diagonal pair)
    kappa_precond = kappa_eq

    return {
        "kappa_direct": kappa_direct,
        "kappa_raw_stable": kappa_raw_stable,
        "kappa_VD_inv": kappa_right,  # right amp only
        "kappa_col_norm": kappa_col,
        "kappa_amp_row": kappa_amp_row,
        "kappa_equilibrated": kappa_eq,
        "kappa_precond": kappa_precond,  # headline figure
        "amp_ratio": ratio,
        "n": float(V.shape[1]),
    }

def fit_loglog(Ms: Sequence[int], kappas: Sequence[float]) -> Tuple[float, float]:
    """Least-squares log κ ≈ a + b log M. Returns (a, b)."""
    xs: List[float] = []
    ys: List[float] = []
    for M, kap in zip(Ms, kappas):
        if kap is None:
            continue
        if not math.isfinite(kap) or kap <= 0:
            continue
        xs.append(math.log(float(M)))
        ys.append(math.log(float(kap)))
    n = len(xs)
    if n < 2:
        return float("nan"), float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    b = num / den if den > 0 else float("nan")
    a = my - b * mx
    return a, b


def diagonal_majorant_proxy(
    M: int,
    Y0: float = Y0_DEFAULT,
    r: float = R_DEFAULT,
    theta: float = THETA_DEFAULT,
) -> float:
    """
    Legacy single-cusp diagonal majorant ratio max w / min w
    (matches plot_kappa_vs_M.kappa_majorant_series spirit).
    """
    diags: List[float] = []
    for n in range(1, M + 1):
        if r2_exact_gaussian(n) == 0:
            continue
        y = 2.0 * math.pi * math.sqrt(n) * Y0
        k2 = k_bessel_amp(y, r) ** 2
        w = (n ** (2.0 * theta)) * max(k2, 1e-300)
        diags.append(w)
    if not diags:
        return float("nan")
    return max(diags) / min(diags)


def run_diagnostic(
    M_values: Sequence[int],
    Y0: float = Y0_DEFAULT,
    r: float = R_DEFAULT,
    theta: float = THETA_DEFAULT,
    two_cusp: bool = True,
    coupling: float = 0.1,
    height_mismatch: float = 0.0,
) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []
    for M in M_values:
        t0 = time.perf_counter()
        print(f"\n=== M = {M} ===")
        V, amps, modes = build_single_cusp_V(
            M, Y0=Y0, r=r, theta=theta, height_mismatch=height_mismatch
        )
        print(
            f"  modes = {len(modes)}  shape V = {V.shape}  "
            f"(pts × modes)"
        )
        stats = kappa_raw_and_precond(V, amps)
        proxy = diagonal_majorant_proxy(M, Y0=Y0, r=r, theta=theta)
        stats["M"] = float(M)
        stats["kappa_diag_proxy"] = proxy
        stats["Y0"] = Y0
        stats["r"] = r
        stats["n_pts"] = float(V.shape[0])
        print(
            f"  amp_ratio            = {stats['amp_ratio']:.6e}\n"
            f"  κ_diag_proxy         = {proxy:.6e}\n"
            f"  κ_raw_stable         = {stats['kappa_raw_stable']:.6e}\n"
            f"  κ_direct             = {stats['kappa_direct']:.6e}\n"
            f"  κ(V D_amp^{{-1}})      = {stats['kappa_VD_inv']:.6e}\n"
            f"  κ(col-Jacobi)        = {stats['kappa_col_norm']:.6e}\n"
            f"  κ(amp+row Jacobi)    = {stats['kappa_amp_row']:.6e}\n"
            f"  κ(equilibrated)      = {stats['kappa_equilibrated']:.6e}\n"
            f"  κ_precond (headline) = {stats['kappa_precond']:.6e}"
        )
        if two_cusp:
            # Exact two-cusp preconditioned SVD when mode count is moderate.
            stats["n2"] = 2.0 * stats["n"]
            c = abs(coupling)
            couple_factor = (1.0 + c) / max(1.0 - c, 1e-12)
            if stats["n"] <= 700:
                kap2_pre = two_cusp_precond_kappa(V, amps, coupling=coupling)
                # Equilibrate the two-cusp preconditioned block for fairness
                U = V / amps[np.newaxis, :]
                top = np.concatenate([U, -c * U], axis=1)
                bot = np.concatenate([-c * U, U], axis=1)
                B = np.concatenate([top, bot], axis=0)
                kap2_eq = cond_2(equilibrate(B, n_iter=6))
                stats["kappa_VD_inv_2cusp"] = kap2_pre
                stats["kappa_precond_2cusp"] = kap2_eq
                stats["kappa_raw_stable_2cusp"] = stats["amp_ratio"] * kap2_pre
                stats["two_cusp_exact_svd"] = 1.0
            else:
                stats["kappa_VD_inv_2cusp"] = stats["kappa_VD_inv"] * couple_factor
                stats["kappa_precond_2cusp"] = stats["kappa_precond"] * couple_factor
                stats["kappa_raw_stable_2cusp"] = (
                    stats["amp_ratio"] * stats["kappa_VD_inv_2cusp"]
                )
                stats["two_cusp_exact_svd"] = 0.0
            print(
                f"  [2-cusp stub] n_unknowns={int(stats['n2'])}  "
                f"κ_raw≈{stats['kappa_raw_stable_2cusp']:.6e}  "
                f"κ(amp)={stats['kappa_VD_inv_2cusp']:.6e}  "
                f"κ(eq)={stats['kappa_precond_2cusp']:.6e}  "
                f"exact_svd={bool(stats['two_cusp_exact_svd'])}"
            )
        stats["seconds"] = time.perf_counter() - t0
        print(f"  time = {stats['seconds']:.2f}s")
        rows.append(stats)
    return rows


def write_csv(path: str, rows: List[Dict[str, float]], two_cusp: bool) -> None:
    fields = [
        "M",
        "n_modes",
        "n_pts",
        "amp_ratio",
        "kappa_diag_proxy",
        "kappa_raw_stable",
        "kappa_direct",
        "kappa_VD_inv",
        "kappa_col_norm",
        "kappa_amp_row",
        "kappa_equilibrated",
        "kappa_precond",
        "log_M",
        "log_kappa_raw",
        "log_kappa_precond",
    ]
    if two_cusp:
        fields += [
            "kappa_raw_stable_2cusp",
            "kappa_VD_inv_2cusp",
            "kappa_precond_2cusp",
        ]
    fields += ["seconds", "Y0", "r"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            kap_r = row["kappa_raw_stable"]
            kap_p = row["kappa_precond"]
            out = {
                "M": int(row["M"]),
                "n_modes": int(row["n"]),
                "n_pts": int(row.get("n_pts", row["n"])),
                "amp_ratio": f"{row['amp_ratio']:.16e}",
                "kappa_diag_proxy": f"{row['kappa_diag_proxy']:.16e}",
                "kappa_raw_stable": f"{kap_r:.16e}",
                "kappa_direct": f"{row['kappa_direct']:.16e}",
                "kappa_VD_inv": f"{row['kappa_VD_inv']:.16e}",
                "kappa_col_norm": f"{row['kappa_col_norm']:.16e}",
                "kappa_amp_row": f"{row['kappa_amp_row']:.16e}",
                "kappa_equilibrated": f"{row['kappa_equilibrated']:.16e}",
                "kappa_precond": f"{kap_p:.16e}",
                "log_M": f"{math.log(row['M']):.16e}",
                "log_kappa_raw": (
                    f"{math.log(kap_r):.16e}"
                    if kap_r > 0 and math.isfinite(kap_r)
                    else "nan"
                ),
                "log_kappa_precond": (
                    f"{math.log(kap_p):.16e}"
                    if kap_p > 0 and math.isfinite(kap_p)
                    else "nan"
                ),
                "seconds": f"{row['seconds']:.6f}",
                "Y0": row["Y0"],
                "r": row["r"],
            }
            if two_cusp:
                out["kappa_raw_stable_2cusp"] = f"{row['kappa_raw_stable_2cusp']:.16e}"
                out["kappa_VD_inv_2cusp"] = f"{row['kappa_VD_inv_2cusp']:.16e}"
                out["kappa_precond_2cusp"] = f"{row['kappa_precond_2cusp']:.16e}"
            w.writerow(out)


def write_plot(
    path: str,
    rows: List[Dict[str, float]],
    b_raw: float,
    b_pre: float,
    b_amp: float,
    Y0: float,
    r: float,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Ms = [int(row["M"]) for row in rows]
    raw = [row["kappa_raw_stable"] for row in rows]
    pre = [row["kappa_precond"] for row in rows]
    amp = [row["kappa_VD_inv"] for row in rows]
    proxy = [row["kappa_diag_proxy"] for row in rows]

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.5))

    ax = axes[0]
    ax.loglog(Ms, raw, "o-", color="C3", label=r"$\kappa_{\mathrm{raw}}$ (stable)")
    ax.loglog(Ms, proxy, "s--", color="C0", alpha=0.8, label=r"$\kappa_{\mathrm{diag\ proxy}}$")
    ax.set_xlabel(r"$M$")
    ax.set_ylabel(r"$\kappa$")
    ax.set_title(
        rf"Raw / diagonal-proxy (no useful preconditioning)"
        f"\n$Y_0={Y0}$, $r={r:.5f}$, fit $b_{{\\mathrm{{raw}}}}\\approx{b_raw:.2f}$"
    )
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.loglog(Ms, pre, "o-", color="C2", label=r"$\kappa$ equilibrated (headline)")
    ax.loglog(Ms, amp, "s--", color="C1", alpha=0.85, label=r"$\kappa(V D_{\mathrm{amp}}^{-1})$ only")
    if pre[0] > 0 and math.isfinite(pre[0]):
        M0, k0 = Ms[0], pre[0]
        ax.loglog(
            Ms,
            [k0 * (M / M0) ** 4 for M in Ms],
            "--",
            color="gray",
            alpha=0.7,
            label=r"$\propto M^{4}$ (Rung 3 threshold)",
        )
        ax.loglog(
            Ms,
            [k0 * (M / M0) for M in Ms],
            ":",
            color="gray",
            alpha=0.7,
            label=r"$\propto M$",
        )
    ax.set_xlabel(r"$M$")
    ax.set_ylabel(r"$\kappa$ (preconditioned)")
    ax.set_title(
        rf"Diagonal preconditioners"
        f"\n$b_{{\\mathrm{{eq}}}}\\approx{b_pre:.3f}$, "
        f"$b_{{\\mathrm{{amp}}}}\\approx{b_amp:.3f}$  "
        + (r"($b_{\mathrm{eq}}<4$ PASS)" if math.isfinite(b_pre) and b_pre < 4 else r"($b_{\mathrm{eq}}<4$ FAIL)")
    )
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(fontsize=8)

    fig.suptitle(
        "Hejhal model collocation conditioning (float diagnostic, non-certifying)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Hejhal κ(V) vs M conditioning diagnostic")
    p.add_argument(
        "--M",
        type=str,
        default="100,200,400,800",
        help="comma-separated M values (default 100,200,400,800)",
    )
    p.add_argument("--Y0", type=float, default=Y0_DEFAULT)
    p.add_argument("--r", type=float, default=R_DEFAULT)
    p.add_argument("--theta", type=float, default=THETA_DEFAULT)
    p.add_argument("--two-cusp", action="store_true", default=True)
    p.add_argument("--no-two-cusp", action="store_true")
    p.add_argument("--coupling", type=float, default=0.1)
    p.add_argument(
        "--height-mismatch",
        type=float,
        default=0.0,
        help="blend in true y_γ pull-back weights (0=same-height Rung3 model)",
    )
    args = p.parse_args(list(argv) if argv is not None else None)

    two_cusp = not args.no_two_cusp
    M_values = [int(x.strip()) for x in args.M.split(",") if x.strip()]
    here = os.path.dirname(os.path.abspath(__file__))

    print("Hejhal collocation conditioning diagnostic (NON-CERTIFYING float)")
    print(f"  Y0={args.Y0}, r={args.r}, theta={args.theta}")
    print(f"  M values = {M_values}")
    print(f"  height_mismatch = {args.height_mismatch}")
    print(f"  K-Bessel backend: {'scipy.kve' if _HAS_KVE else 'elementary majorant'}")
    print(f"  two-cusp stub: {two_cusp}")
    print(f"  lemma_K import: {_HAS_LEMMA_K}")

    rows = run_diagnostic(
        M_values,
        Y0=args.Y0,
        r=args.r,
        theta=args.theta,
        two_cusp=two_cusp,
        coupling=args.coupling,
        height_mismatch=args.height_mismatch,
    )

    # Fits on the Milestone-2 set
    Ms = [int(row["M"]) for row in rows]
    a_raw, b_raw = fit_loglog(Ms, [row["kappa_raw_stable"] for row in rows])
    a_amp, b_amp = fit_loglog(Ms, [row["kappa_VD_inv"] for row in rows])
    a_pre, b_pre = fit_loglog(Ms, [row["kappa_precond"] for row in rows])
    a_proxy, b_proxy = fit_loglog(Ms, [row["kappa_diag_proxy"] for row in rows])
    a_col, b_col = fit_loglog(Ms, [row["kappa_col_norm"] for row in rows])
    a_ar, b_ar = fit_loglog(Ms, [row["kappa_amp_row"] for row in rows])

    print("\n" + "=" * 72)
    print("log-log fits  log κ ≈ a + b log M")
    print(f"  diagonal proxy:     a={a_proxy:.6f}  b={b_proxy:.6f}   (κ ~ M^{b_proxy:.3f})")
    print(f"  raw stable:         a={a_raw:.6f}  b={b_raw:.6f}   (κ ~ M^{b_raw:.3f})")
    print(f"  κ(V D_amp^{{-1}}):    a={a_amp:.6f}  b={b_amp:.6f}   (κ ~ M^{b_amp:.3f})")
    print(f"  col-Jacobi:         a={a_col:.6f}  b={b_col:.6f}   (κ ~ M^{b_col:.3f})")
    print(f"  amp+row Jacobi:     a={a_ar:.6f}  b={b_ar:.6f}   (κ ~ M^{b_ar:.3f})")
    print(f"  equilibrated:       a={a_pre:.6f}  b={b_pre:.6f}   (κ ~ M^{b_pre:.3f})")
    if two_cusp and rows and "kappa_precond_2cusp" in rows[0]:
        a2, b2 = fit_loglog(Ms, [row["kappa_precond_2cusp"] for row in rows])
        print(f"  2-cusp equilibrated:a={a2:.6f}  b={b2:.6f}   (κ ~ M^{b2:.3f})")

    print("\nRung 3 stopping condition: after diagonal preconditioner, b < 4")
    print(f"  HEADLINE (equilibrated): b = {b_pre:.4f}", end="")
    if math.isfinite(b_pre) and b_pre < 4:
        print(" < 4  → PASS (diagnostic only)")
    else:
        print("  → FAIL / inconclusive")
    print(f"  amp-only right D:        b = {b_amp:.4f}")
    if math.isfinite(b_raw):
        print(f"  RAW (stable):            b = {b_raw:.4f}  (naive ~76 without precond)")
    print(
        "\nNOTE: This is a float diagnostic. Rung 3 is NOT complete without\n"
        "      well-posedness proof, σ0 consistency, and Krawczyk for N=5."
    )

    csv_path = os.path.join(here, "hejhal_kappa_vs_M.csv")
    write_csv(csv_path, rows, two_cusp=two_cusp)
    print(f"\nWrote {csv_path}")

    png_path = os.path.join(here, "hejhal_kappa_vs_M.png")
    try:
        write_plot(png_path, rows, b_raw, b_pre, b_amp, args.Y0, args.r)
        print(f"Wrote {png_path}")
    except Exception as exc:
        print(f"Plot skipped ({exc})")

    # Machine-readable summary line for reports
    summary_path = os.path.join(here, "hejhal_kappa_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("Hejhal conditioning diagnostic summary (NON-CERTIFYING)\n")
        f.write(f"Y0={args.Y0} r={args.r} theta={args.theta}\n")
        f.write(f"M={M_values}\n")
        f.write(f"b_diag_proxy={b_proxy}\n")
        f.write(f"b_raw_stable={b_raw}\n")
        f.write(f"b_amp_right={b_amp}\n")
        f.write(f"b_col_jacobi={b_col}\n")
        f.write(f"b_amp_row={b_ar}\n")
        f.write(f"b_equilibrated={b_pre}\n")
        f.write(
            f"rung3_precond_b_lt_4="
            f"{'YES' if math.isfinite(b_pre) and b_pre < 4 else 'NO'}\n"
        )
        f.write("krawczyk_N5=NOT_RUN\n")
        f.write("rung3_complete=NO\n")
    print(f"Wrote {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
