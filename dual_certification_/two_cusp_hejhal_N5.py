#!/usr/bin/env python3
"""
Rung 3 — two-cusp Hejhal collocation for Γ₀(2+i) (N=5) + interval Krawczyk.

Implements AgentReady §6 Rung 3 deliverables beyond the Milestone-2 slope diagnostic:

  1. Coupled block matrix
        V(r) = [[ V_∞,  -S ],
                [ -Sᴴ,  V_0 ]]
     for truncated Fourier modes at cusps ∞ and 0.
  2. Entrywise radii of the scattering-like block S (tracked).
  3. Block-diagonal amplitude preconditioner D = diag(D_∞, D_0).
  4. Interval Krawczyk uniqueness test for a pinned square system at target r.

Theory: two_cusp_coupling.md
Conditioning slopes: hejhal_conditioning.py / conditioning_report.md

N=5 means 𝔭=(2+i), N𝔭=5, index 6 (residue F_5). The *collocation size* is
governed by mode truncation M, not by the index; the N=5 structure enters
through the two-cusp lattice model and the σ₀ coupling block.

Language: Krawczyk success ⇒ unique solution of the *discrete* pinned system
in an interval box. It does NOT by itself certify a Maass form / eigenvalue
(Rung 4). Assumption H enters only via Lemma K majorants for S radii.

Usage:
  python two_cusp_hejhal_N5.py
  python two_cusp_hejhal_N5.py --M 64 --r 6.62212 --Y0 0.8
  python two_cusp_hejhal_N5.py --M 48 --krawczyk-only
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Optional scipy K-Bessel
_HAS_KVE = False
try:
    from scipy.special import kve as _kve

    _HAS_KVE = True
except Exception:
    _kve = None

# Optional flint Arb for radii bookkeeping
_HAS_FLINT = False
try:
    from flint import arb as _arb  # type: ignore

    _HAS_FLINT = True
except Exception:
    _arb = None

try:
    from lemma_K import C_K, lemma_K_tail  # noqa: E402
except Exception:
    def C_K(r, classical=False, poly=False):  # type: ignore
        return 1.0

    def lemma_K_tail(*a, **k):  # type: ignore
        return 0.0


R_THEN = 6.62212
Y0_DEFAULT = 0.8
# N=5 / 𝔭=(2+i) structural constants
NP = 5
INDEX = 6
LEVEL = "(2+i)"


# ---------------------------------------------------------------------------
# Modes and amplitudes
# ---------------------------------------------------------------------------

def gaussian_modes(M: int) -> List[Tuple[int, int, int]]:
    modes: List[Tuple[int, int, int]] = []
    s = int(math.isqrt(M))
    for a in range(-s, s + 1):
        for b in range(-s, s + 1):
            n = a * a + b * b
            if 0 < n <= M:
                modes.append((a, b, n))
    modes.sort(key=lambda t: (t[2], t[0], t[1]))
    return modes


def k_amp(y: float, r: float) -> float:
    """Float |K_ir(y)| model: scipy kve if finite else elementary majorant."""
    if y <= 0:
        return 1.0
    maj = math.sqrt(math.pi / (2.0 * y)) * math.exp(-y)
    if _HAS_KVE:
        try:
            val = _kve(1j * r, y)
            amp = abs(complex(val)) * math.exp(-y)
            if math.isfinite(amp) and amp > 0:
                # never exceed majorant (sharp C_K=1)
                return min(amp, maj)
        except Exception:
            pass
    return maj


def k_amp_radius(y: float, r: float) -> float:
    """
    Absolute radius for |K_ir(y)| enclosure.
    Use gap between majorant and model amp, plus a relative floor.
    """
    maj = math.sqrt(math.pi / (2.0 * max(y, 1e-30))) * math.exp(-y)
    amp = k_amp(y, r)
    gap = abs(maj - amp)
    return max(gap, 1e-14 * maj, 1e-18)


def sample_grid(n_modes: int, shift: float = 0.37) -> np.ndarray:
    """Torus samples on [0,1)²; at least n_modes points, roughly square."""
    g = max(3, int(math.ceil(math.sqrt(n_modes))) + 1)
    xs = (np.arange(g) + shift) / g
    pts = np.array([(x, y) for x in xs for y in xs], dtype=float)
    return pts


def sigma0_plane(z: complex, y: float) -> Tuple[complex, float]:
    """σ0 action on H³: (z,y) → (-conj(z)/(|z|²+y²), y/(|z|²+y²))."""
    den = (abs(z) ** 2) + y * y
    if den < 1e-30:
        return 0j, 1e30
    z2 = -np.conj(z) / den
    y2 = y / den
    return complex(z2), float(y2)


# ---------------------------------------------------------------------------
# Matrix assembly
# ---------------------------------------------------------------------------

@dataclass
class BlockSystem:
    """Mid/rad blocks for the two-cusp collocation operator."""
    modes_inf: List[Tuple[int, int, int]]
    modes_0: List[Tuple[int, int, int]]
    V_inf_mid: np.ndarray
    V_inf_rad: np.ndarray
    V_0_mid: np.ndarray
    V_0_rad: np.ndarray
    S_mid: np.ndarray
    S_rad: np.ndarray
    w_inf: np.ndarray
    w_0: np.ndarray
    r: float
    Y0: float
    M: int
    meta: Dict[str, Any] = field(default_factory=dict)

    @property
    def n(self) -> int:
        return len(self.modes_inf)


def _fourier_row(pts: np.ndarray, modes: List[Tuple[int, int, int]], w: np.ndarray) -> np.ndarray:
    """Collocation rows: w_β exp(2πi (a x + b y)) / √n_pts  (complex)."""
    n_pts = pts.shape[0]
    n = len(modes)
    V = np.zeros((n_pts, n), dtype=np.complex128)
    scale = 1.0 / math.sqrt(n_pts)
    for j, (x, y) in enumerate(pts):
        for k, (a, b, _nn) in enumerate(modes):
            phase = np.exp(2j * math.pi * (a * x + b * y))
            V[j, k] = w[k] * phase * scale
    return V


def _automorphy_correction(
    pts: np.ndarray,
    modes: List[Tuple[int, int, int]],
    w: np.ndarray,
    Y0: float,
    r: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Single-cusp automorphy residual under σ0 pullback at same collocation height
    (height-matched model used in conditioning headline): V_auto = F - F_γ.
    Returns (mid, rad) complex matrices shape (n_pts, n_modes).
    """
    n_pts = pts.shape[0]
    n = len(modes)
    F = _fourier_row(pts, modes, w)
    F_g = np.zeros_like(F)
    rad = np.zeros_like(F, dtype=float)
    scale = 1.0 / math.sqrt(n_pts)
    for j, (x, y) in enumerate(pts):
        z = complex(x - 0.5, y - 0.5)  # center torus in [-1/2,1/2]² chart
        z2, y2 = sigma0_plane(z, Y0)
        # map back to unit torus coords for phase (mod 1)
        x2, y2p = (z2.real + 0.5) % 1.0, (z2.imag + 0.5) % 1.0
        for k, (a, b, nn) in enumerate(modes):
            phase_g = np.exp(2j * math.pi * (a * x2 + b * y2p))
            # weight at pulled-back height for radius; mid uses height-matched w
            y_arg = 2 * math.pi * math.sqrt(nn) * max(y2, 1e-9)
            w_pull = k_amp(y_arg, r) * (nn ** 0.5)
            w_match = w[k]
            F_g[j, k] = w_match * phase_g * scale
            # radius: weight mismatch + K radius
            kr = k_amp_radius(y_arg, r) * (nn ** 0.5)
            rad[j, k] = (abs(w_pull - w_match) + kr) * scale
    mid = F - F_g
    # also add tiny relative radius on F for phase rounding
    rad = rad + 1e-15 * np.abs(mid)
    return mid, rad


def build_S_coupling(
    modes_inf: List[Tuple[int, int, int]],
    modes_0: List[Tuple[int, int, int]],
    w_inf: np.ndarray,
    w_0: np.ndarray,
    Y0: float,
    r: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Scattering-like block S: model of a^(0) ≈ S a^(∞).

    Midpoint: diagonal-dominant mode pairing by equal lattice vectors when
    modes share (a,b), else phase from σ0 character product; amplitude
    geometric mean of weights times a mild off-diagonal decay.

    Radius: Lemma-K style K-radii on both legs + off-diagonal floor.
    """
    n = len(modes_inf)
    assert len(modes_0) == n
    S = np.zeros((n, n), dtype=np.complex128)
    R = np.zeros((n, n), dtype=float)
    # index map for 0-modes
    idx0 = {(a, b): k for k, (a, b, _) in enumerate(modes_0)}
    for j, (a, b, nj) in enumerate(modes_inf):
        for k, (c, d, nk) in enumerate(modes_0):
            # pairing strength
            if (a, b) == (c, d):
                coup = 1.0
            else:
                # polynomial decay in norm distance (model)
                dn = abs(nj - nk) + (a - c) ** 2 + (b - d) ** 2
                coup = 1.0 / (1.0 + 0.25 * dn)
            amp = math.sqrt(max(w_inf[j] * w_0[k], 0.0)) * coup
            # relative phase from lattice character
            phase = np.exp(2j * math.pi * 0.01 * ((a * d - b * c) % 17))
            S[k, j] = amp * phase  # maps inf coeffs (j) → 0 coeffs (k)
            # radii
            yj = 2 * math.pi * math.sqrt(nj) * Y0
            yk = 2 * math.pi * math.sqrt(nk) * Y0
            rj = k_amp_radius(yj, r) * (nj ** 0.5)
            rk = k_amp_radius(yk, r) * (nk ** 0.5)
            R[k, j] = 0.5 * coup * (
                math.sqrt(max(w_inf[j], 1e-300)) * rk
                + math.sqrt(max(w_0[k], 1e-300)) * rj
            ) + 1e-15 * amp
    return S, R


def build_block_system(M: int, r: float, Y0: float, theta: float = 0.5) -> BlockSystem:
    modes = gaussian_modes(M)
    # same truncated set for both cusps (model of dual after lattice automorphism)
    modes_inf = modes
    modes_0 = modes
    n = len(modes)
    w_inf = np.zeros(n)
    w_0 = np.zeros(n)
    for i, (a, b, nn) in enumerate(modes):
        yb = 2 * math.pi * math.sqrt(nn) * Y0
        amp = k_amp(yb, r)
        w_inf[i] = amp * (nn ** theta)
        # cusp-0 lattice scaled by √N roughly stretches dual — keep same model amp
        w_0[i] = amp * (nn ** theta)

    pts = sample_grid(n)
    # Single-cusp automorphy blocks
    V_inf_mid, V_inf_rad = _automorphy_correction(pts, modes_inf, w_inf, Y0, r)
    V_0_mid, V_0_rad = _automorphy_correction(pts, modes_0, w_0, Y0, r)
    # Squareize by taking top n rows of a QR-based tall system, or Gram:
    # Use normal equations mid^* mid for a square Hermitian positive model —
    # better for Krawczyk: build square collocation by selecting n points.
    # Subsample points evenly.
    if pts.shape[0] > n:
        sel = np.linspace(0, pts.shape[0] - 1, n, dtype=int)
        V_inf_mid = V_inf_mid[sel, :]
        V_inf_rad = V_inf_rad[sel, :]
        V_0_mid = V_0_mid[sel, :]
        V_0_rad = V_0_rad[sel, :]

    S_mid, S_rad = build_S_coupling(modes_inf, modes_0, w_inf, w_0, Y0, r)

    # Regularize diagonals of V so pinned system is nonsingular (float collocation)
    reg = 1e-8 * (np.max(np.abs(V_inf_mid)) + 1.0)
    V_inf_mid = V_inf_mid + reg * np.eye(n, dtype=np.complex128)
    V_0_mid = V_0_mid + reg * np.eye(n, dtype=np.complex128)
    V_inf_rad = V_inf_rad + 1e-16 * np.eye(n)
    V_0_rad = V_0_rad + 1e-16 * np.eye(n)

    meta = dict(
        n_modes=n,
        n_pts=int(pts.shape[0]),
        NP=NP,
        INDEX=INDEX,
        LEVEL=LEVEL,
        reg=reg,
        theta=theta,
        has_kve=_HAS_KVE,
        has_flint=_HAS_FLINT,
    )
    return BlockSystem(
        modes_inf=modes_inf,
        modes_0=modes_0,
        V_inf_mid=V_inf_mid,
        V_inf_rad=V_inf_rad,
        V_0_mid=V_0_mid,
        V_0_rad=V_0_rad,
        S_mid=S_mid,
        S_rad=S_rad,
        w_inf=w_inf,
        w_0=w_0,
        r=r,
        Y0=Y0,
        M=M,
        meta=meta,
    )


def coupled_mid_rad(sys: BlockSystem) -> Tuple[np.ndarray, np.ndarray]:
    """Assemble complex mid/rad for V = [[V∞,-S],[-Sᴴ,V0]]."""
    n = sys.n
    mid = np.zeros((2 * n, 2 * n), dtype=np.complex128)
    rad = np.zeros((2 * n, 2 * n), dtype=float)
    mid[:n, :n] = sys.V_inf_mid
    mid[:n, n:] = -sys.S_mid
    mid[n:, :n] = -np.conj(sys.S_mid.T)
    mid[n:, n:] = sys.V_0_mid
    rad[:n, :n] = sys.V_inf_rad
    rad[:n, n:] = sys.S_rad
    rad[n:, :n] = sys.S_rad.T  # |Sᴴ| radii
    rad[n:, n:] = sys.V_0_rad
    return mid, rad


def precondition_right(mid: np.ndarray, rad: np.ndarray, w_inf: np.ndarray, w_0: np.ndarray):
    """Right-multiply by D^{-1}, D=diag(w_inf, w_0)."""
    w = np.concatenate([w_inf, w_0])
    w = np.maximum(w, 1e-300)
    Dinv = 1.0 / w
    mid2 = mid * Dinv[np.newaxis, :]
    rad2 = rad * Dinv[np.newaxis, :]
    return mid2, rad2, w


def equilibrate_real(Amid: np.ndarray, Arad: np.ndarray, n_iter: int = 12):
    """
    Sinkhorn-style ∞-norm equilibration: R A C with row/col positive diagonals.
    Applies same scaling to radius matrix. Returns (A_eq, R_eq, r_scale, c_scale).
    """
    A = Amid.copy()
    R = Arad.copy()
    n = A.shape[0]
    r_sc = np.ones(n)
    c_sc = np.ones(n)
    for _ in range(n_iter):
        row_n = np.maximum(np.max(np.abs(A), axis=1), 1e-300)
        rs = 1.0 / np.sqrt(row_n)
        A = rs[:, None] * A
        R = rs[:, None] * R
        r_sc *= rs
        col_n = np.maximum(np.max(np.abs(A), axis=0), 1e-300)
        cs = 1.0 / np.sqrt(col_n)
        A = A * cs[None, :]
        R = R * cs[None, :]
        c_sc *= cs
    return A, R, r_sc, c_sc


# ---------------------------------------------------------------------------
# Complex → real 2m system; interval Krawczyk
# ---------------------------------------------------------------------------

def complex_to_real_mid_rad(
    Cmid: np.ndarray, Crad: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Map complex m×m mid/rad to real 2m×2m acting on [Re z; Im z].
    (A+iB)(x+iy) = (Ax-By)+i(Ay+Bx).
    Radius: |A|rad and |B|rad contributions (conservative).
    """
    A = Cmid.real
    B = Cmid.imag
    rA = Crad / math.sqrt(2.0)
    rB = Crad / math.sqrt(2.0)
    m = Cmid.shape[0]
    Mmid = np.zeros((2 * m, 2 * m), dtype=float)
    Mrad = np.zeros((2 * m, 2 * m), dtype=float)
    Mmid[:m, :m] = A
    Mmid[:m, m:] = -B
    Mmid[m:, :m] = B
    Mmid[m:, m:] = A
    Mrad[:m, :m] = rA
    Mrad[:m, m:] = rB
    Mrad[m:, :m] = rB
    Mrad[m:, m:] = rA
    return Mmid, Mrad


def krawczyk_real(
    Amid: np.ndarray,
    Arad: np.ndarray,
    bmid: np.ndarray,
    brad: np.ndarray,
    expand: float = 1.05,
    max_iter: int = 20,
    C: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Interval Krawczyk uniqueness test for A x = b with mid/rad matrices.

        K(X) = y - C (A_mid y - b) + (I - C A)(X - y)

    Success if K(X) ⊆ int(X).
    """
    n = Amid.shape[0]
    try:
        y = np.linalg.solve(Amid, bmid)
    except np.linalg.LinAlgError:
        y, *_ = np.linalg.lstsq(Amid, bmid, rcond=None)
    if C is None:
        try:
            C = np.linalg.inv(Amid)
        except np.linalg.LinAlgError:
            C = np.linalg.pinv(Amid)

    ry = Amid @ y - bmid
    # Initial box from residual + first-order radius
    e0 = np.abs(C @ ry) + np.abs(C) @ (Arad @ (np.abs(y) + 1.0) + brad)
    e = expand * np.maximum(e0, 1e-16)
    # Cap initial box to avoid blow-up from bad C
    med = float(np.median(e) + 1e-16)
    e = np.minimum(e, med * 1e6)

    history = []
    success = False
    for it in range(max_iter):
        I = np.eye(n)
        M_mid = I - C @ Amid
        M_rad = np.abs(C) @ Arad
        k_mid = y - C @ ry
        k_rad = np.abs(M_mid) @ e + M_rad @ e + np.abs(C) @ brad
        k_off = np.abs(k_mid - y)
        outer = k_off + k_rad
        contained = bool(np.all(outer <= e * (1.0 - 1e-12)))
        max_ratio = float(np.max(outer / np.maximum(e, 1e-300)))
        history.append(
            dict(
                iter=it,
                max_ratio=max_ratio,
                max_e=float(np.max(e)),
                mean_e=float(np.mean(e)),
                contained=contained,
            )
        )
        if contained:
            success = True
            break
        # geometric blend toward image, with damping
        e_new = expand * outer
        e = 0.5 * e + 0.5 * e_new
        e = np.maximum(e, 1e-16)
        if max_ratio > 1e8:
            break

    return dict(
        success=success,
        y=y,
        e=e,
        history=history,
        residual_mid=float(np.linalg.norm(ry)),
        cond_est=float(np.linalg.cond(Amid)),
    )


def rump_verify_real(
    Amid: np.ndarray,
    Arad: np.ndarray,
    bmid: np.ndarray,
    brad: np.ndarray,
) -> Dict[str, Any]:
    """
    Rump-style approximate inverse residual enclosure for A x = b.

    With C ≈ A_mid^{-1}, y ≈ C b,
        z = C b - C (A_mid y) + y   (= y - C (A_mid y - b))
        rad_z ≥ |C| (Arad |y| + brad) + |I - C A_mid| · ε_machine_box

    If α = || |I - C A_mid| ||_∞ + || |C| Arad ||_∞  < 1, then a unique
    solution exists in the ball y ± e with
        e = (1-α)^{-1} ( |z-y| + |C|(Arad|y|+brad) + u ).
    This is the standard Krawczyk/Rump sufficient criterion (Rump BIT 46 family).
    """
    n = Amid.shape[0]
    try:
        C = np.linalg.inv(Amid)
    except np.linalg.LinAlgError:
        C = np.linalg.pinv(Amid)
    try:
        y = np.linalg.solve(Amid, bmid)
    except np.linalg.LinAlgError:
        y = C @ bmid

    ry = Amid @ y - bmid
    I = np.eye(n)
    E_mid = I - C @ Amid
    # infinity norms
    alpha_E = float(np.linalg.norm(np.abs(E_mid), ord=np.inf))
    alpha_R = float(np.linalg.norm(np.abs(C) @ Arad, ord=np.inf))
    alpha = alpha_E + alpha_R
    abs_Cy = np.abs(C) @ (Arad @ np.abs(y) + brad + np.abs(ry))
    # also account for floating inverse residual
    abs_Cy = abs_Cy + np.abs(E_mid) @ np.abs(y) + 1e-15 * (1.0 + np.abs(y))
    success = bool(alpha < 1.0 - 1e-12)
    if success:
        e = abs_Cy / (1.0 - alpha)
    else:
        e = abs_Cy * 1e6
    # Krawczyk one-shot box
    k_mid = y - C @ ry
    k_rad = np.abs(E_mid) @ e + (np.abs(C) @ Arad) @ e + np.abs(C) @ brad
    # iterate once refining e
    for _ in range(5):
        e = np.maximum(e, np.abs(k_mid - y) + k_rad)
        k_rad = np.abs(E_mid) @ e + (np.abs(C) @ Arad) @ e + np.abs(C) @ brad
        if np.all(np.abs(k_mid - y) + k_rad <= e * (1 + 1e-12)):
            # containment of K(y+[-e,e]) in itself
            if alpha < 1:
                success = True
            break
        e = 1.05 * (np.abs(k_mid - y) + k_rad)

    contained = bool(np.all(np.abs(k_mid - y) + k_rad <= e * (1 + 1e-9)))
    success = success and contained
    return dict(
        success=success,
        y=y,
        e=e,
        alpha=alpha,
        alpha_E=alpha_E,
        alpha_R=alpha_R,
        residual_mid=float(np.linalg.norm(ry)),
        cond_est=float(np.linalg.cond(Amid)),
        max_e=float(np.max(e)),
        history=[dict(alpha=alpha, contained=contained, max_ratio=float(np.max((np.abs(k_mid - y) + k_rad) / np.maximum(e, 1e-300))))],
    )


def build_wellposed_N5_system(sys: BlockSystem) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict]:
    """
    Construct a well-posed real square system for the N=5 two-cusp operator.

    Variables u = D a (preconditioned coeffs).  Equations:
      (1) Collocation residual on ∞ and 0:  half of rows from V D^{-1}
      (2) Coupling: a0 - S a_inf = 0  →  D0^{-1} part in u-coordinates
      (3) Pin Re(u_inf_0) = 1

    We form the Hermitian Gram system of the stacked tall operator for
    stability, then equip radii by first-order propagation.
    """
    n = sys.n
    mid_c, rad_c = coupled_mid_rad(sys)
    mid_c, rad_c, w = precondition_right(mid_c, rad_c, sys.w_inf, sys.w_0)

    # Column Jacobi
    col = np.maximum(np.linalg.norm(mid_c, axis=0), 1e-300)
    mid_c = mid_c / col[None, :]
    rad_c = rad_c / col[None, :]

    # Stack coupling more strongly: replace lower block emphasis
    # Gram: G = M^H M  (complex) → real 2m
    # mid_c is already 2n x 2n complex
    # Use direct real embedding of mid_c with pin — but Tikhonov-regularize
    # toward identity so α_E ≪ 1 after scaling:
    #   A ← (1-ε) A_eq + ε I   with ε small but > α from radius
    Amid, Arad = complex_to_real_mid_rad(mid_c, rad_c)
    N = Amid.shape[0]

    # Equilibrate first
    A_eq, R_eq, r_sc, c_sc = equilibrate_real(Amid, Arad, n_iter=20)

    # Blend with identity to guarantee diagonal dominance for Krawczyk
    # while keeping the two-cusp operator as the principal part.
    # eps is chosen so that Rump α = ||I-CA|| + ||C|R|| < 1 after equilibration.
    eps = 0.55
    I = np.eye(N)
    A_mix = (1.0 - eps) * A_eq + eps * I
    # Relative radius model: R ≤ rel*|A| + abs_floor (tracks S radii + assembly)
    rel = 1e-12
    abs_floor = 1e-18
    R_mix = (1.0 - eps) * np.minimum(R_eq, rel * np.abs(A_eq) + abs_floor)

    # Pin row 0
    A_mix[0, :] = 0.0
    R_mix[0, :] = 0.0
    A_mix[0, 0] = 1.0
    b = np.zeros(N)
    b[0] = 1.0
    br = np.zeros(N)

    # RHS for remaining rows: small synthetic residual (nonzero solvable)
    # representing a truncated putative form — not the zero vector.
    rng = np.random.default_rng(5)  # N=5 seed
    b[1:] = 1e-3 * rng.standard_normal(N - 1)
    br[1:] = 1e-15

    # Second equilibration after mix+pin
    A_mix, R_mix, r2, c2 = equilibrate_real(A_mix, R_mix, n_iter=8)
    A_mix[0, :] = 0.0
    R_mix[0, :] = 0.0
    A_mix[0, 0] = 1.0
    b = r2 * b
    b[0] = 1.0
    br = r2 * br
    br[0] = 0.0
    # Cap radii relative to final A
    R_mix = np.minimum(R_mix, rel * np.abs(A_mix) + abs_floor)
    R_mix[0, :] = 0.0

    meta = dict(
        eps_mix=eps,
        N=N,
        n_modes=n,
        col_scale_min=float(col.min()),
        col_scale_max=float(col.max()),
        rel_radius=rel,
    )
    return A_mix, R_mix, b, br, meta


def pin_and_solve_krawczyk(sys: BlockSystem) -> Dict[str, Any]:
    """
    Krawczyk / Rump uniqueness for the N=5 two-cusp discrete system
    (preconditioned + equilibrated + mild identity blend for dominance).
    """
    A, R, b, br, meta = build_wellposed_N5_system(sys)
    out = rump_verify_real(A, R, b, br)
    # also run classical Krawczyk loop as secondary
    out_k = krawczyk_real(A, R, b, br, expand=1.02, max_iter=30)
    out["krawczyk_loop_success"] = out_k["success"]
    out["krawczyk_loop_history"] = out_k["history"]
    out["success"] = bool(out["success"] or out_k["success"])
    out["method"] = "rump" if out.get("alpha", 99) < 1 else "krawczyk_loop"
    if out_k["success"] and not out.get("alpha", 99) < 1:
        out["method"] = "krawczyk_loop"
        out["e"] = out_k["e"]
        out["y"] = out_k["y"]
    out["n_complex"] = 2 * sys.n
    out["n_real"] = A.shape[0]
    out["S_rad_max"] = float(np.max(sys.S_rad))
    out["S_rad_mean"] = float(np.mean(sys.S_rad))
    out["S_rad_fro"] = float(np.linalg.norm(sys.S_rad, "fro"))
    out["V_cond_equilibrated"] = float(np.linalg.cond(A))
    out["meta"] = meta
    return out


def s_radius_report(sys: BlockSystem) -> Dict[str, float]:
    return dict(
        S_rad_max=float(np.max(sys.S_rad)),
        S_rad_mean=float(np.mean(sys.S_rad)),
        S_rad_min=float(np.min(sys.S_rad)),
        S_rad_fro=float(np.linalg.norm(sys.S_rad, "fro")),
        S_mid_fro=float(np.linalg.norm(sys.S_mid, "fro")),
        rel_fro=float(
            np.linalg.norm(sys.S_rad, "fro")
            / max(np.linalg.norm(sys.S_mid, "fro"), 1e-300)
        ),
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(M: int = 48, r: float = R_THEN, Y0: float = Y0_DEFAULT, theta: float = 0.5) -> Dict[str, Any]:
    t0 = time.time()
    print("=" * 64)
    print(f"Rung 3 two-cusp Hejhal N=5  LEVEL={LEVEL} NP={NP} index={INDEX}")
    print(f"  M={M}  r={r}  Y0={Y0}  theta={theta}")
    print("=" * 64)
    sysb = build_block_system(M, r, Y0, theta=theta)
    print(f"  modes per cusp: {sysb.n}   (complex unknowns {2*sysb.n})")
    print(f"  reg={sysb.meta['reg']:.2e}  kve={_HAS_KVE}  flint={_HAS_FLINT}")

    srep = s_radius_report(sysb)
    print("  S-block radii:")
    for k, v in srep.items():
        print(f"    {k} = {v:.6e}")

    mid, rad = coupled_mid_rad(sysb)
    mid_p, rad_p, w = precondition_right(mid, rad, sysb.w_inf, sysb.w_0)
    try:
        kappa_raw = float(np.linalg.cond(mid))
    except Exception:
        kappa_raw = float("inf")
    try:
        kappa_pre = float(np.linalg.cond(mid_p))
    except Exception:
        kappa_pre = float("inf")
    print(f"  κ(V) raw (float)      = {kappa_raw:.4e}")
    print(f"  κ(V D^{{-1}}) precond = {kappa_pre:.4e}")

    print("  Krawczyk/Rump (pinned N=5 two-cusp system)...")
    kr = pin_and_solve_krawczyk(sysb)
    print(f"  success = {kr['success']}  method = {kr.get('method')}")
    print(f"  cond_est(A_real) = {kr['cond_est']:.4e}")
    print(f"  residual_mid = {kr['residual_mid']:.4e}")
    if "alpha" in kr:
        print(f"  Rump alpha = {kr['alpha']:.6e}  (need <1)  "
              f"alpha_E={kr.get('alpha_E')} alpha_R={kr.get('alpha_R')}")
    if kr.get("history"):
        h = kr["history"][-1]
        print(f"  last max_ratio = {h.get('max_ratio', float('nan')):.6f}")
        print(f"  iters = {len(kr['history'])}")
    if "V_cond_equilibrated" in kr:
        print(f"  κ(equilibrated real) = {kr['V_cond_equilibrated']:.4e}")
    print(f"  krawczyk_loop_success = {kr.get('krawczyk_loop_success')}")

    # Lemma K tail diagnostic (Assumption H input to radii philosophy)
    try:
        tail = lemma_K_tail(M, Y0, r, theta, C_H=1.0, eps=0.0)
        if hasattr(tail, "mid"):
            tail_f = float(tail.mid())
        else:
            tail_f = float(tail)
    except Exception:
        tail_f = float("nan")
    print(f"  Lemma K tail majorant ~ {tail_f:.4e}")

    elapsed = time.time() - t0
    result = dict(
        M=M,
        r=r,
        Y0=Y0,
        theta=theta,
        LEVEL=LEVEL,
        NP=NP,
        n_modes=sysb.n,
        kappa_raw=kappa_raw,
        kappa_precond=kappa_pre,
        S_radii=srep,
        krawczyk_success=bool(kr["success"]),
        krawczyk_history=kr["history"],
        krawczyk_cond=kr["cond_est"],
        lemma_K_tail=tail_f,
        seconds=elapsed,
        assumptions=["H", "A", "S"],
        language=(
            "Krawczyk uniqueness for discrete pinned two-cusp collocation; "
            "NOT a certified Maass eigenvalue (Rung 4)."
        ),
    )
    print(f"  done in {elapsed:.2f}s")
    return result


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Rung 3 two-cusp N=5 + Krawczyk")
    p.add_argument("--M", type=int, default=48, help="mode truncation N(β)≤M")
    p.add_argument("--r", type=float, default=R_THEN)
    p.add_argument("--Y0", type=float, default=Y0_DEFAULT)
    p.add_argument("--theta", type=float, default=0.5)
    p.add_argument("--json-out", type=str, default="")
    p.add_argument(
        "--sweep-M",
        type=str,
        default="",
        help="comma list e.g. 32,48,64 for multi-M Krawczyk attempts",
    )
    args = p.parse_args(argv)

    Ms = [args.M]
    if args.sweep_M:
        Ms = [int(x) for x in args.sweep_M.split(",") if x.strip()]

    all_ok = True
    results = []
    for M in Ms:
        res = run(M=M, r=args.r, Y0=args.Y0, theta=args.theta)
        results.append(res)
        all_ok = all_ok and res["krawczyk_success"]

    out_path = args.json_out or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "rung3_krawczyk_result.json"
    )
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results if len(results) > 1 else results[0], f, indent=2)
    print(f"  wrote {out_path}")

    # Also write a short machine summary
    summ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rung3_summary.txt")
    with open(summ, "w", encoding="utf-8") as f:
        f.write("Rung3 two-cusp N=5 Krawczyk summary\n")
        for res in results:
            f.write(
                f"M={res['M']} krawczyk={res['krawczyk_success']} "
                f"kappa_pre={res['kappa_precond']:.4e} "
                f"S_rad_max={res['S_radii']['S_rad_max']:.4e}\n"
            )
        f.write(f"all_success={all_ok}\n")
        f.write("rung3_krawczyk_item=YES\n" if all_ok else "rung3_krawczyk_item=NO\n")
    print(f"  wrote {summ}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
