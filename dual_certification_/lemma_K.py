#!/usr/bin/env python3
"""
Lemma K — uniform K-Bessel tail with Hecke control (Arb-enclosable).

Under Assumption H(K,θ): |a_β| ≤ C_H(ε) N(β)^{θ+ε} for every ε>0,
with N(β)=|β|², prove an explicit majorant for

    S_{M,∞}(r,Y0) = ∑_{N(β)>M} |a_β|² |K_{ir}(2π|β|Y0)|²

and enclose it as an Arb ball when python-flint is available.

This module is self-contained: no FEM, no gluing, no N𝔭, no reference cell.

Mathematical summary (matches theorem_DK.tex §4)
------------------------------------------------
1. Pointwise UPPER bound (elementary, from integral representation):
       |K_{ir}(y)| ≤ √(π/(2y)) · e^{-y} · C_K(r),   y>0, r≥0,
   with sharp C_K(r) = 1 proved by
       |K_{ir}| ≤ K_0 ≤ K_{1/2} = √(π/(2y)) e^{-y}.
   Classical looser constant C_K^{cl}(r) = exp(π r / 2) also works.
   NO r-independent lower bound is claimed (false for large r at fixed y).

2. Substituting y = 2π √n Y0, n = N(β):
       |K_{ir}(2π√n Y0)|² ≤ (1/(4 √n Y0)) exp(-4π √n Y0) C_K(r)².

3. Hecke control |a_β|² ≤ C_H² n^{2(θ+ε)} and grouping by norm:
       S ≤ [C_H² C_K(r)² /(4 Y0)] ∑_{n>M} r₂(n) n^{2θ+2ε−1/2} e^{−4π√n Y0}.

4. Multiplicity (theorem default): r₂(n) ≤ 6 d(n) ≤ 6 n.
   Optional looser: r₂(n) ≤ 6 n directly (mode "crude").
   After absorbing r₂ ≤ 6n: α = 2θ+2ε+1/2, prefactor ×6.

5. Integral comparison (monotonicity threshold N_mono stated):
       ∑_{n≥N0} n^α exp(-c √n) ≤ f(N0) + ∫_{N0}^∞ f
   with closed form via upper incomplete gamma, N0 ≥ N_mono.

API
---
    C_K(r)                         -> Arb/float enclosure of C_K(r)
    lemma_K_tail(M,Y0,r,theta,...) -> Arb enclosure of the majorant bound
    tail_majorant(...)             -> alias of lemma_K_tail (AgentReady name)
    C1_C2_constants(...)           -> named C1, C2 matching the paper
    benchmark_majorant(M_list,...) -> print analytic majorant vs truncated sum
    luke_upper_validation(...)     -> 20-pair Luke-type UPPER self-test

Usage:
    python lemma_K.py --test --bench
    python lemma_K.py --M 400 --Y0 0.8 --r 6.62212 --theta 0.5
"""
from __future__ import annotations

import argparse
import math
import random
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

# ---------------------------------------------------------------------------
# Optional Arb backend (python-flint)
# ---------------------------------------------------------------------------
HAS_FLINT = False
HAS_ACB = False
_arb = None  # type: ignore
_acb = None  # type: ignore
_ctx = None  # type: ignore
try:
    from flint import arb as _arb  # type: ignore
    from flint import ctx as _ctx  # type: ignore

    HAS_FLINT = True
    try:
        from flint import acb as _acb  # type: ignore

        HAS_ACB = True
    except Exception:  # pragma: no cover
        _acb = None  # type: ignore
except Exception:  # pragma: no cover
    _arb = None  # type: ignore
    _acb = None  # type: ignore
    _ctx = None  # type: ignore


# Then's published first spectral parameter for PSL(2,Z[i]) (engineering target).
THEN_R1 = 6.62212

# Working precision for certifying Arb paths (bits). Diagnostics may use lower.
DEFAULT_ARB_PREC = 128

Number = Union[float, Any]  # float or flint.arb


def set_arb_precision(prec: int = DEFAULT_ARB_PREC) -> Optional[int]:
    """Set Arb working precision in bits. Returns previous prec, or None if no flint."""
    if not HAS_FLINT or _ctx is None:
        return None
    prev = int(_ctx.prec)
    _ctx.prec = int(prec)
    return prev


def _is_arb(x: Any) -> bool:
    return HAS_FLINT and type(x).__name__ == "arb"


def to_ball(x: Union[float, int, str, Any], radius: float = 0.0) -> Number:
    """Wrap a real number as an Arb ball if flint is available, else float."""
    if HAS_FLINT:
        if _is_arb(x):
            return x
        a = _arb(x)
        if radius > 0:
            a = a + _arb(0, radius)
        return a
    return float(x)


def arb_pi() -> Number:
    """π as Arb ball at current precision (or float fallback)."""
    if HAS_FLINT:
        return _arb.pi()
    return math.pi


def ball_mid(x: Number) -> float:
    if _is_arb(x):
        try:
            return float(x.mid())
        except Exception:
            return float(x)
    return float(x)


def ball_rad(x: Number) -> float:
    if _is_arb(x):
        try:
            return float(x.rad())
        except Exception:
            return 0.0
    return 0.0


def ball_str(x: Number, digits: int = 12) -> str:
    if _is_arb(x):
        return str(x)
    return f"{float(x):.{digits}g}  [non-certifying float]"


# ---------------------------------------------------------------------------
# C_K(r) — explicit K-Bessel prefactor (UPPER bound only)
# ---------------------------------------------------------------------------
def C_K(
    r: Union[float, Any],
    poly: bool = False,
    classical: bool = False,
) -> Number:
    """
    Explicit majorant prefactor C_K(r) for the UPPER bound

        |K_{ir}(y)| ≤ √(π/(2y)) e^{-y} C_K(r),   y>0, r≥0.

    Default (sharp, elementary proof in theorem_DK.tex Lemma Kpoint):
        C_K(r) = 1
    via |K_{ir}| ≤ K_0 ≤ K_{1/2} = √(π/(2y)) e^{-y}.

    classical=True (looser literature constant):
        C_K^{cl}(r) = exp(π r / 2) ≥ 1.

    poly=True (optional UPPER envelope for y≥1 only):
        sharp path:   1 + 1/8
        classical:    exp(π r / 2) · (1 + (1+r²)/4)
    No r-independent LOWER bound is provided or claimed.

    Returns an Arb ball when python-flint is present (certifying), else float.
    """
    r = to_ball(r)
    if classical:
        if HAS_FLINT:
            half_pi_r = (_arb(math.pi) * r) / 2
            base = half_pi_r.exp()
            if poly:
                return base * (1 + (1 + r * r) / 4)
            return base
        r_f = float(r)
        base = math.exp(math.pi * r_f / 2.0)
        if poly:
            return base * (1.0 + (1.0 + r_f * r_f) / 4.0)
        return base

    # Sharp C_K = 1
    if poly:
        # y≥1 envelope from K_0 ≤ √(π/(2y)) e^{-y} (1 + 1/(8y)) ≤ (9/8) · main term
        return to_ball(1.0 + 1.0 / 8.0)
    return to_ball(1.0)


def k_bessel_pointwise_majorant(
    y: float,
    r: float,
    poly: bool = False,
    classical: bool = False,
) -> Number:
    """
    Majorant for |K_{ir}(y)| (UPPER bound only):
        √(π/(2y)) * exp(-y) * C_K(r).

    Never a lower bound.  Referee forbids claiming
    |K_ir| ≥ √(π/(2y)) e^{-y} independent of r.
    """
    if y <= 0:
        raise ValueError("y must be positive")
    ck = C_K(r, poly=poly, classical=classical)
    if HAS_FLINT:
        yb = to_ball(y)
        pref = (arb_pi() / (to_ball(2.0) * yb)).sqrt() * (-yb).exp()
        return pref * ck
    pref = math.sqrt(math.pi / (2.0 * y)) * math.exp(-y)
    return pref * float(ck)


def k0_arb(y: float) -> Number:
    """K_0(y) via Arb bessel_k when available, else scipy/float diagnostic."""
    if y <= 0:
        raise ValueError("y must be positive")
    if HAS_FLINT:
        return to_ball(y).bessel_k(to_ball(0))
    try:
        from scipy.special import k0 as sc_k0  # type: ignore

        return float(sc_k0(y))
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("K_0 evaluation requires python-flint or scipy") from exc


def kir_abs_arb(r: float, y: float) -> Optional[Number]:
    """
    |K_{ir}(y)| via Acb when available.

    Returns None if complex-order K_ν is unavailable (no acb backend).
    """
    if y <= 0:
        raise ValueError("y must be positive")
    if not (HAS_FLINT and HAS_ACB and _acb is not None):
        return None
    # K_ν(z) with ν = i r, z = y > 0
    nu = _acb(0, r)
    z = _acb(y)
    return abs(z.bessel_k(nu))


def luke_upper_validation(
    n_samples: int = 20,
    seed: int = 20260712,
    r_range: Tuple[float, float] = (0.1, 12.0),
    y_range: Tuple[float, float] = (0.5, 8.0),
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Self-test: Luke-type UPPER bounds for random (r,y) pairs.

    Validated claims (UPPER only):
      (U1) |K_ir(y)| ≤ √(π/(2y)) e^{-y}           (sharp C_K=1), when K_ir available
      (U2) K_0(y)     ≤ √(π/(2y)) e^{-y}
      (U3) for y≥1:   K_0(y) ≤ √(π/(2y)) e^{-y} (1 + 1/(8y))   (Luke poly refinement)
      (U4) K_0(y) ≤ K_{1/2}(y) = √(π/(2y)) e^{-y}

    Explicitly NOT claimed / counterexample recorded:
      (L⊥) r-independent lower bound |K_ir| ≥ √(π/(2y)) e^{-y}  is FALSE.

    If K_ir is unavailable in the backend, (U1) is skipped and (U2)–(U4) on the
    K_0 path still run (the certified tail only needs the K_0 ≤ K_{1/2} majorant).
    """
    rng = random.Random(seed)
    set_arb_precision(DEFAULT_ARB_PREC)

    results: List[dict] = []
    n_kir_ok = 0
    n_k0_ok = 0
    n_poly_ok = 0
    n_kir_tested = 0
    n_poly_tested = 0
    kir_available = HAS_FLINT and HAS_ACB

    if verbose:
        print("Luke-type UPPER validation (never lower bound)")
        print(f"  n_samples={n_samples}, seed={seed}")
        print(f"  HAS_FLINT={HAS_FLINT}, HAS_ACB={HAS_ACB} (K_ir path={'on' if kir_available else 'off'})")
        print(f"  arb prec bits = {_ctx.prec if _ctx is not None else 'n/a'}")

    for i in range(n_samples):
        r = rng.uniform(*r_range)
        y = rng.uniform(*y_range)
        maj = k_bessel_pointwise_majorant(y, r, poly=False, classical=False)
        maj_poly = k_bessel_pointwise_majorant(y, r, poly=True, classical=False)
        k0 = k0_arb(y)

        # (U2)(U4): K_0 ≤ majorant with C_K=1
        # Ball check: maj - K_0 has nonnegative lower bound (or mid comparison with rad slack).
        if HAS_FLINT and _is_arb(k0) and _is_arb(maj):
            gap_k0 = maj - k0
            k0_ok = ball_mid(gap_k0) >= -ball_rad(gap_k0)  # 0 ∈ (-∞, gap] upper-ok
            # stronger: lower endpoint of gap ≥ 0  <=>  mid - rad ≥ 0
            k0_ok = (ball_mid(gap_k0) - ball_rad(gap_k0)) >= -1e-30
        else:
            k0_ok = float(ball_mid(k0)) <= float(ball_mid(maj)) * (1.0 + 1e-9)
        if k0_ok:
            n_k0_ok += 1

        poly_ok = True
        if y >= 1.0:
            n_poly_tested += 1
            if HAS_FLINT and _is_arb(k0) and _is_arb(maj_poly):
                gap_p = maj_poly - k0
                poly_ok = (ball_mid(gap_p) - ball_rad(gap_p)) >= -1e-30
            else:
                poly_ok = float(ball_mid(k0)) <= float(ball_mid(maj_poly)) * (1.0 + 1e-9)
            if poly_ok:
                n_poly_ok += 1

        kir_ok = True
        kir_val: Optional[Number] = None
        if kir_available:
            n_kir_tested += 1
            kir_val = kir_abs_arb(r, y)
            assert kir_val is not None
            if HAS_FLINT and _is_arb(maj):
                gap = maj - kir_val
                kir_ok = (ball_mid(gap) - ball_rad(gap)) >= -1e-30
            else:
                kir_ok = float(ball_mid(kir_val)) <= float(ball_mid(maj)) * (1.0 + 1e-9)
            if kir_ok:
                n_kir_ok += 1

        row = {
            "r": r,
            "y": y,
            "K0": k0,
            "maj": maj,
            "k0_ok": k0_ok,
            "poly_ok": poly_ok,
            "kir": kir_val,
            "kir_ok": kir_ok,
        }
        results.append(row)
        if verbose and i < 5:
            kir_s = ball_str(kir_val) if kir_val is not None else "n/a"
            print(
                f"  [{i}] r={r:.4f} y={y:.4f} |Kir|={kir_s}  "
                f"K0={ball_str(k0)}  maj={ball_str(maj)}  "
                f"U_K0={k0_ok} U_Kir={kir_ok}"
            )

    # Counterexample to false lower bound: large r, moderate y ⇒ |Kir| ≪ majorant
    r_bad, y_bad = 6.62212, 1.0
    maj_bad = k_bessel_pointwise_majorant(y_bad, r_bad)
    kir_bad = kir_abs_arb(r_bad, y_bad) if kir_available else None
    lower_bound_false = False
    if kir_bad is not None:
        # If |Kir| < 0.1 * majorant, the r-free lower bound fails badly
        lower_bound_false = ball_mid(kir_bad) < 0.1 * ball_mid(maj_bad)

    if verbose:
        print(
            f"  K0≤maj: {n_k0_ok}/{n_samples}  "
            f"poly y≥1: {n_poly_ok}/{n_poly_tested}  "
            f"|Kir|≤maj: {n_kir_ok}/{n_kir_tested if n_kir_tested else 'skipped'}"
        )
        if kir_bad is not None:
            print(
                f"  False-lower-bound check at r={r_bad}, y={y_bad}: "
                f"|Kir|={ball_str(kir_bad)}  maj={ball_str(maj_bad)}  "
                f"ratio={ball_mid(kir_bad)/ball_mid(maj_bad):.6g}  "
                f"(r-free lower bound fails: {lower_bound_false})"
            )
        elif not kir_available:
            print(
                "  K_ir unavailable in backend: validated Luke UPPER via K_0 path only "
                "(matches certified chain |K_ir|≤K_0≤K_{1/2})."
            )

    all_pass = (
        n_k0_ok == n_samples
        and (n_poly_tested == 0 or n_poly_ok == n_poly_tested)
        and (n_kir_tested == 0 or n_kir_ok == n_kir_tested)
        and (not kir_available or lower_bound_false)
    )
    if verbose:
        print(f"  Luke UPPER validation: {'PASS' if all_pass else 'FAIL'}")

    return {
        "pass": all_pass,
        "n_samples": n_samples,
        "n_k0_ok": n_k0_ok,
        "n_poly_ok": n_poly_ok,
        "n_poly_tested": n_poly_tested,
        "n_kir_ok": n_kir_ok,
        "n_kir_tested": n_kir_tested,
        "kir_available": kir_available,
        "lower_bound_false": lower_bound_false,
        "results": results,
        "backend": (
            "python-flint Arb+Acb"
            if kir_available
            else ("python-flint Arb (K_0 only)" if HAS_FLINT else "float/scipy")
        ),
    }


# ---------------------------------------------------------------------------
# Multiplicity r₂(n)
# ---------------------------------------------------------------------------
def d_of_n(n: int) -> int:
    """Number of positive divisors of n (exact, integer)."""
    if n <= 0:
        return 0
    c = 0
    i = 1
    while i * i < n:
        if n % i == 0:
            c += 2
        i += 1
    if i * i == n:
        c += 1
    return c


def r2_exact_gaussian(n: int) -> int:
    """
    Exact r₂(n) = #{(a,b)∈ℤ² : a²+b² = n} = 4(χ₄*1)(n) = 4(d1−d3).
    Used only for diagnostics / truncated majorant sums, not required for
    the certified crude majorant.
    """
    if n < 0:
        return 0
    if n == 0:
        return 1
    d1 = d3 = 0
    i = 1
    while i * i < n:
        if n % i == 0:
            for div in (i, n // i):
                if div % 4 == 1:
                    d1 += 1
                elif div % 4 == 3:
                    d3 += 1
        i += 1
    if i * i == n:
        if i % 4 == 1:
            d1 += 1
        elif i % 4 == 3:
            d3 += 1
    return 4 * (d1 - d3)


def r2_majorant(n: int, mode: str = "div") -> Number:
    """
    Explicit upper bound on r₂_K(n) valid for both Z[i] and Z[ω].

    mode:
      'div'    : 6 d(n)  (theorem default; r₂ ≤ 4 d(n) for Z[i]; 6 d for Z[ω])
      'crude'  : 6 n     (optional looser bound)
      'exact_i': exact r₂ for Z[i] (diagnostic only; not a majorant for Z[ω])
    """
    if n <= 0:
        return to_ball(0)
    if mode == "crude":
        return to_ball(6 * n)
    if mode == "div":
        return to_ball(6 * d_of_n(n))
    if mode == "exact_i":
        return to_ball(r2_exact_gaussian(n))
    raise ValueError(f"unknown r2 mode: {mode}")


# ---------------------------------------------------------------------------
# Incomplete-gamma integral comparison
# ---------------------------------------------------------------------------
def incomplete_gamma_upper(s: Number, z: Number) -> Number:
    """
    Upper incomplete gamma Γ(s,z) = ∫_z^∞ t^{s-1} e^{-t} dt.

    Arb: z.gamma_upper(s).  Float fallback: scipy if present, else asymptotic.
    """
    if HAS_FLINT:
        zs = to_ball(z)
        ss = to_ball(s)
        return zs.gamma_upper(ss)

    s_f, z_f = float(s), float(z)
    try:
        from scipy.special import gammaincc, gamma as sc_gamma  # type: ignore

        return float(sc_gamma(s_f) * gammaincc(s_f, z_f))
    except Exception:
        if z_f <= 0:
            return float("inf")
        term = 1.0
        total = 1.0
        for k in range(1, 40):
            term *= (s_f - k) / z_f
            total += term
            if abs(term) < 1e-16 * abs(total):
                break
        return (z_f ** (s_f - 1.0)) * math.exp(-z_f) * total


def integral_x_alpha_exp(c: Number, alpha: Number, x0: Number) -> Number:
    """
    Exact identity:
        ∫_{x0}^∞ x^α exp(-c √x) dx = 2 c^{-(2α+2)} Γ(2α+2, c √x0),
    valid for c>0, x0>0.
    """
    if HAS_FLINT:
        c = to_ball(c)
        alpha = to_ball(alpha)
        x0 = to_ball(x0)
        s = 2 * alpha + 2
        z = c * x0.sqrt()
        return 2 * (c ** (-(2 * alpha + 2))) * incomplete_gamma_upper(s, z)

    c_f, a_f, x0_f = float(c), float(alpha), float(x0)
    s = 2 * a_f + 2
    z = c_f * math.sqrt(x0_f)
    return 2.0 * (c_f ** (-(2 * a_f + 2))) * float(incomplete_gamma_upper(s, z))


def monotonicity_threshold(alpha: float, c: float) -> int:
    """
    N_mono for f(x)=x^α e^{-c√x}: decreasing on [N_mono, ∞).
    For α>0: N_mono = ceil((2α/c)²)+1; for α≤0: N_mono = 1.
    """
    if alpha > 0 and c > 0:
        return int(math.ceil((2.0 * alpha / c) ** 2)) + 1
    return 1


def series_majorant_integral(
    N0: int,
    alpha: Number,
    c: Number,
    r2_mode: str = "div",
) -> Number:
    """
    Bound ∑_{n≥N0} n^α exp(-c √n) by integral comparison:

        ∑_{n≥N0} f(n) ≤ f(N0) + ∫_{N0}^∞ f(x) dx

    after raising N0 to the monotonicity threshold when needed.
    Finite head [N0, N_mono) is summed in Arb when available.
    """
    if N0 < 1:
        N0 = 1

    a_f = ball_mid(alpha)
    c_f = ball_mid(c)
    N_mono = max(N0, monotonicity_threshold(a_f, c_f))

    if HAS_FLINT:
        alpha_b = to_ball(alpha)
        c_b = to_ball(c)
        head = to_ball(0.0)
        for n in range(N0, N_mono):
            nb = to_ball(n)
            term = (nb ** alpha_b) * (-c_b * nb.sqrt()).exp()
            head = head + term
        n0 = to_ball(N_mono)
        fN = (n0 ** alpha_b) * (-c_b * n0.sqrt()).exp()
        integ = integral_x_alpha_exp(c_b, alpha_b, n0)
        return head + fN + integ

    head = 0.0
    for n in range(N0, N_mono):
        nf = float(n)
        head += (nf ** a_f) * math.exp(-c_f * math.sqrt(nf))
    fN = (float(N_mono) ** a_f) * math.exp(-c_f * math.sqrt(float(N_mono)))
    integ = float(integral_x_alpha_exp(c_f, a_f, float(N_mono)))
    return float(head) + fN + integ


# ---------------------------------------------------------------------------
# Main lemma_K_tail
# ---------------------------------------------------------------------------
def lemma_K_tail(
    M: int,
    Y0: float,
    r: float,
    theta: float,
    C_H: float = 1.0,
    eps: float = 0.0,
    r2_mode: str = "div",
    poly_CK: bool = False,
    classical_CK: bool = False,
    use_integral: bool = True,
    sum_cap: Optional[int] = None,
) -> Number:
    """
    Arb enclosure of the analytic majorant for S_{M,∞}(r,Y0).

    Parameters
    ----------
    M : int
        Truncation in the norm: sum over n = N(β) > M.
    Y0 : float
        Collar height (≥ 0.5 recommended; proof assumes Y0 > 0).
    r : float
        Spectral parameter, λ = r² + 1, r ≥ 0.
    theta : float
        Hecke exponent in Assumption H (e.g. 0.5, 7/64, 0).
    C_H : float
        Hecke constant C_H(ε) from Assumption H.
    eps : float
        The ε in N(β)^{θ+ε}.
    r2_mode : {'div','crude','exact_i'}
        Multiplicity majorant.  Default 'div' (r₂ ≤ 6 d(n) ≤ 6n).
    poly_CK : bool
        Use polynomial-refined C_K envelope (y≥1).
    classical_CK : bool
        Use C_K^{cl}=exp(π r/2) instead of sharp C_K=1.
    use_integral : bool
        If True, enclose the infinite tail by incomplete-gamma comparison.
    sum_cap : optional int
        When use_integral is False, sum n = M+1 .. sum_cap (diagnostic).

    Returns
    -------
    Arb ball (certifying) or float (non-certifying) upper bound on S_{M,∞}.
    """
    if M < 0:
        raise ValueError("M must be ≥ 0")
    if Y0 <= 0:
        raise ValueError("Y0 must be positive")
    if r < 0:
        raise ValueError("r must be ≥ 0")
    if theta < 0:
        raise ValueError("theta must be ≥ 0")

    ck = C_K(r, poly=poly_CK, classical=classical_CK)
    if HAS_FLINT:
        pref = (to_ball(C_H) ** 2) * (ck ** 2) / (to_ball(4.0) * to_ball(Y0))
        beta_power = to_ball(2.0 * theta + 2.0 * eps - 0.5)
        c = to_ball(4.0) * arb_pi() * to_ball(Y0)  # exp(-c √n)
    else:
        pref = (C_H ** 2) * (float(ck) ** 2) / (4.0 * Y0)
        beta_power = 2.0 * theta + 2.0 * eps - 0.5
        c = 4.0 * math.pi * Y0

    N0 = M + 1
    if N0 < 1:
        N0 = 1

    if not use_integral:
        cap = sum_cap if sum_cap is not None else max(N0 + 5000, 10 * N0)
        ssum: Number = to_ball(0.0) if HAS_FLINT else 0.0
        for n in range(N0, cap + 1):
            r2m = r2_majorant(n, mode=r2_mode)
            if HAS_FLINT:
                nb = to_ball(n)
                term_core = (nb ** beta_power) * (-c * nb.sqrt()).exp()
                ssum = ssum + r2m * term_core
            else:
                term_core = (float(n) ** float(beta_power)) * math.exp(
                    -float(c) * math.sqrt(float(n))
                )
                ssum = ssum + float(r2m) * term_core
        return pref * ssum

    # Certified route: r₂(n) ≤ 6 n  (also covers div mode since 6d(n)≤6n)
    #   r₂ n^{β} ≤ 6 n^{β+1}  ⇒  bound 6 · ∑ n^{β+1} e^{-c√n}.
    if r2_mode in ("crude", "div"):
        alpha = beta_power + to_ball(1.0) if HAS_FLINT else (float(beta_power) + 1.0)
        mult = 6.0
        series = series_majorant_integral(N0, alpha, c)
        return pref * to_ball(mult) * series

    if r2_mode == "exact_i":
        head_end = min(N0 + 20000, max(N0 + 1000, 5 * N0))
        ssum = to_ball(0.0) if HAS_FLINT else 0.0
        for n in range(N0, head_end + 1):
            r2e = r2_exact_gaussian(n)
            if r2e == 0:
                continue
            if HAS_FLINT:
                nb = to_ball(n)
                term = to_ball(r2e) * (nb ** beta_power) * (-c * nb.sqrt()).exp()
                ssum = ssum + term
            else:
                term = r2e * (float(n) ** float(beta_power)) * math.exp(
                    -float(c) * math.sqrt(float(n))
                )
                ssum = ssum + term
        alpha = beta_power + to_ball(1.0) if HAS_FLINT else (float(beta_power) + 1.0)
        tail = to_ball(6.0) * series_majorant_integral(head_end + 1, alpha, c)
        return pref * (ssum + tail)

    raise ValueError(f"unknown r2_mode: {r2_mode}")


def tail_majorant(
    M: int,
    Y0: float,
    r: float,
    theta: float,
    C_H: float = 1.0,
    eps: float = 0.0,
    r2_mode: str = "div",
    poly_CK: bool = False,
    classical_CK: bool = False,
    use_integral: bool = True,
    sum_cap: Optional[int] = None,
) -> Number:
    """
    AgentReady API alias for lemma_K_tail.

    Returns an Arb ball enclosing the analytic tail majorant ε(M,Y0,r,θ)
    under Assumption H with parameters (C_H, θ, ε).  (C_ε in the paper is
    tracked here as the pair (C_H, eps) with exponent θ+ε.)
    """
    return lemma_K_tail(
        M,
        Y0,
        r,
        theta,
        C_H=C_H,
        eps=eps,
        r2_mode=r2_mode,
        poly_CK=poly_CK,
        classical_CK=classical_CK,
        use_integral=use_integral,
        sum_cap=sum_cap,
    )


def truncated_majorant_sum(
    M: int,
    Y0: float,
    r: float,
    theta: float,
    C_H: float = 1.0,
    eps: float = 0.0,
    Nmax: int = 20000,
    r2_mode: str = "exact_i",
    poly_CK: bool = False,
    classical_CK: bool = False,
) -> Number:
    """
    Direct truncated sum of the *analytic* majorant (not Then's a_β):

        pref · ∑_{M < n ≤ Nmax} r₂(n) n^{2θ+2ε−1/2} e^{−4π√n Y0}

    Default r2_mode='exact_i' uses the exact representation function r₂ for
    ℤ[i] (AgentReady Rung 0 benchmark requirement).  Used as a diagnostic
    tightness check against lemma_K_tail / tail_majorant (integral comparison).
    """
    ck = C_K(r, poly=poly_CK, classical=classical_CK)
    if HAS_FLINT:
        pref = (to_ball(C_H) ** 2) * (ck ** 2) / (to_ball(4.0) * to_ball(Y0))
        beta_power = to_ball(2.0 * theta + 2.0 * eps - 0.5)
        c = to_ball(4.0) * arb_pi() * to_ball(Y0)
    else:
        pref = (C_H ** 2) * (float(ck) ** 2) / (4.0 * Y0)
        beta_power = 2.0 * theta + 2.0 * eps - 0.5
        c = 4.0 * math.pi * Y0

    ssum: Number = to_ball(0.0) if HAS_FLINT else 0.0
    for n in range(M + 1, Nmax + 1):
        if r2_mode == "exact_i":
            r2v = r2_exact_gaussian(n)
        else:
            r2v = int(ball_mid(r2_majorant(n, mode=r2_mode)))
        if r2v == 0:
            continue
        if HAS_FLINT:
            nb = to_ball(n)
            term = to_ball(r2v) * (nb ** beta_power) * (-c * nb.sqrt()).exp()
            ssum = ssum + term
        else:
            term = r2v * (float(n) ** float(beta_power)) * math.exp(
                -float(c) * math.sqrt(float(n))
            )
            ssum = ssum + term
    return pref * ssum


# ---------------------------------------------------------------------------
# C1, C2 — named constants matching theorem_DK.tex §5
# ---------------------------------------------------------------------------
def C1_C2_constants(
    Y: float = 1.25,
    Y0: float = 0.8,
    r: float = THEN_R1,
    theta: float = 0.5,
    C_H: float = 1.0,
    eps: float = 0.0,
    M: int = 400,
    field: str = "i",
    D_Y: Optional[float] = None,
    h_min: Optional[float] = None,
    y_min: Optional[float] = None,
    classical_CK: bool = False,
    sharp_geom: bool = True,
    U_norm: float = 1.0,
    trace_factor: Optional[float] = None,
    bdry_sqrt_lam: bool = True,
) -> Dict[str, Number]:
    """
    Evaluate named constants for Theorem D(K).

    Parameters
    ----------
    field : {'i', 'omega'}
    classical_CK : bool
        If True, use C_K^{cl}=exp(π r/2) in A_cusp (inflates; literature only).
        Default False: sharp C_K=1 from Lemma K (integral rep + K_0 ≤ K_{1/2}).
    sharp_geom : bool
        If True (default): y_min = 1/√2 (Gaussian Humbert floor), D_Y = Euclidean
        core-box diam √(dx²+dy²+(Y−y_min)²), h_min from geometry, trace_factor=1
        so C_tr = |F|/|T| = 3/h_min (morningbrief §1.1–1.2).
        If False: conservative y_min=1/2, D_Y=10, h_min=0.25, factor=4/3 → C_tr=4/h_min.
    U_norm : float
        Explicit ‖U‖₂ of the periodized trial (default 1). Relaxes
        η₀ = min( 1/(4(1+A_bdry)²),  U_norm²/(4 C1²) )
        so the spectral pigeonhole is not forced to unit L² when ‖U‖ is known.
    bdry_sqrt_lam : bool
        If True (default, theorem as written): A_bdry includes (1+√λ).
        If False: drop (1+√λ) — **audit / hypothetical only** (morningbrief §1.6);
        not the production default until the proof is rewritten without double count.
    """
    if field in ("i", "Zi", "Z[i]", "gaussian"):
        T_abs = 0.5
        # Conservative floor y≥1/2; sharp uses Humbert floor y≥1/√2 on |z|²≤1/2.
        y_min_cons = 0.5
        y_min_sharp = 1.0 / math.sqrt(2.0)  # ≈0.7071 (morningbrief §1.1)
        # Picard core tet altitude scale for |F|/|T|=3/h
        h_min_sharp = 0.45
        box_dx, box_dy = 1.0, 0.5  # torus cell extents in (x1,x2)
    elif field in ("omega", "w", "Zomega", "Z[omega]", "eisenstein"):
        T_abs = math.sqrt(3.0) / 6.0
        y_min_cons = math.sqrt(2.0 / 3.0)  # ≈0.816
        y_min_sharp = math.sqrt(2.0 / 3.0)
        h_min_sharp = 0.50
        box_dx, box_dy = 0.5, math.sqrt(3.0) / 3.0  # P3 strip extents
    else:
        raise ValueError(f"unknown field: {field}")

    if y_min is None:
        y_min = y_min_sharp if sharp_geom else y_min_cons

    if sharp_geom:
        if D_Y is None:
            # morningbrief §1.2: Euclidean diameter of computational core box
            #   D_Y = √(dx² + dy² + (Y − y_min)²)
            # (replaces crude D_Y=3; still ≤3 for Y≲1.5, y_min≳0.5)
            dy_h = max(float(Y) - float(y_min), 0.0)
            D_Y = math.sqrt(box_dx ** 2 + box_dy ** 2 + dy_h ** 2)
        if h_min is None:
            h_min = h_min_sharp
        if trace_factor is None:
            # |F|/|T| = 3/h for right tet of altitude h (user sharp form)
            trace_factor = 1.0
            C_trace = 3.0 / h_min
        else:
            C_trace = trace_factor * (3.0 / h_min)
    else:
        if D_Y is None:
            D_Y = 10.0
        if h_min is None:
            h_min = 0.25
        if trace_factor is None:
            trace_factor = 4.0 / 3.0  # Young-expanded form in paper
        C_trace = trace_factor * (3.0 / h_min)  # = 4/h_min when factor=4/3

    lam = r * r + 1.0
    sqrt_lam = math.sqrt(lam)

    C_Poincare = D_Y / math.pi
    C_Sob = 1.0 + C_Poincare

    # Metric comparison (Lemma M) — explicit y_min, Y
    C_met_0 = y_min ** (-1.5)
    C_met_1 = y_min ** (-0.5)
    C_met_0_inv = Y ** 1.5
    C_met_1_inv = Y ** 0.5
    C_met = max(C_met_0, C_met_1, C_met_0_inv, C_met_1_inv)

    A_met = C_met
    bdry_ell = (1.0 + sqrt_lam) if bdry_sqrt_lam else 1.0
    A_bdry = C_trace * C_Sob * A_met * bdry_ell
    A_ell = 1.0 + lam
    A_res = A_ell * (1.0 + C_Poincare * A_met)

    ck = C_K(r, classical=classical_CK)
    if HAS_FLINT:
        beta_power = to_ball(2.0 * theta + 2.0 * eps - 0.5)
        c = to_ball(4.0) * arb_pi() * to_ball(Y0)
        alpha = beta_power + to_ball(1.0)
    else:
        beta_power = 2.0 * theta + 2.0 * eps - 0.5
        c = 4.0 * math.pi * Y0
        alpha = beta_power + 1.0
    N0 = max(M + 1, 1)
    series = series_majorant_integral(N0, alpha, c)
    sum_maj = to_ball(6.0) * series
    if HAS_FLINT:
        A_cusp = (to_ball(C_H) * ck / (to_ball(2.0) * to_ball(Y0).sqrt())) * sum_maj.sqrt()
    else:
        A_cusp = (C_H * float(ball_mid(ck)) / (2.0 * math.sqrt(Y0))) * math.sqrt(
            float(ball_mid(sum_maj))
        )

    A_cut = 2.0 * T_abs * (1.0 / Y) * (1.0 + C_Sob) * A_met

    U2 = max(float(U_norm), 1e-30) ** 2
    if HAS_FLINT:
        one = to_ball(1.0)
        C1 = to_ball(2.0) * to_ball(A_bdry) * to_ball(A_res) * (one + A_cusp)
        C2 = to_ball(A_cut)
        c1m = float(ball_mid(C1))
        # η₀ = min( 1/(4(1+A_bdry)²) ,  U_norm² / (4 C1²) )
        # second term: C1 √η ≤ U_norm/2  ⇒  η ≤ U_norm²/(4 C1²)
        eta0 = min(
            1.0 / (4.0 * (1.0 + A_bdry) ** 2),
            U2 / (4.0 * c1m ** 2),
        )
        eta0_b = to_ball(eta0)
    else:
        C1 = 2.0 * A_bdry * A_res * (1.0 + float(ball_mid(A_cusp)))
        C2 = A_cut
        eta0_b = min(
            1.0 / (4.0 * (1.0 + A_bdry) ** 2),
            U2 / (4.0 * C1 ** 2),
        )

    return {
        "field": field,
        "T_abs": to_ball(T_abs),
        "y_min": to_ball(y_min),
        "Y": to_ball(Y),
        "Y0": to_ball(Y0),
        "r": to_ball(r),
        "lambda": to_ball(lam),
        "C_K": ck,
        "C_trace": to_ball(C_trace),
        "C_Poincare": to_ball(C_Poincare) if HAS_FLINT else C_Poincare,
        "C_Sob": to_ball(C_Sob) if HAS_FLINT else C_Sob,
        "C_met": to_ball(C_met),
        "A_met": to_ball(A_met),
        "A_bdry": to_ball(A_bdry),
        "A_ell": to_ball(A_ell),
        "A_res": to_ball(A_res),
        "A_cusp": A_cusp,
        "A_cut": to_ball(A_cut),
        "C1": C1,
        "C2": C2,
        "eta0": eta0_b if HAS_FLINT else eta0_b,
        "U_norm": to_ball(U_norm),
        "D_Y": to_ball(D_Y),
        "h_min": to_ball(h_min),
        "sharp_geom": sharp_geom,
        "bdry_sqrt_lam": bdry_sqrt_lam,
        "backend": "python-flint Arb" if HAS_FLINT else "float (NON-CERTIFYING)",
        "classical_CK": classical_CK,
    }


def print_C1_C2(field: str = "i", **kwargs: Any) -> Dict[str, Number]:
    """Pretty-print C1_C2_constants for one field."""
    d = C1_C2_constants(field=field, **kwargs)
    print(f"--- C1/C2 constants  field={field}  backend={d['backend']} ---")
    for key in (
        "T_abs",
        "y_min",
        "C_K",
        "C_trace",
        "C_Poincare",
        "C_Sob",
        "C_met",
        "A_met",
        "A_bdry",
        "A_ell",
        "A_res",
        "A_cusp",
        "A_cut",
        "C1",
        "C2",
        "eta0",
    ):
        print(f"  {key:12s} = {ball_str(d[key])}")
    return d


# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------
def benchmark_majorant(
    M_list: Sequence[int] = (100, 200, 400),
    Y0: float = 0.8,
    r: float = THEN_R1,
    theta: float = 0.5,
    C_H: float = 1.0,
    eps: float = 0.0,
    classical_CK: bool = False,
    Nmax: int = 20000,
) -> List[dict]:
    """
    Compare the infinite-tail majorant enclosure (tail_majorant / lemma_K_tail)
    against a truncated majorant sum with exact r₂ for ℤ[i], n ≤ Nmax.
    """
    set_arb_precision(DEFAULT_ARB_PREC)
    rows = []
    backend = "python-flint Arb (certifying)" if HAS_FLINT else "float (NON-CERTIFYING)"
    print("=" * 72)
    print("Lemma K majorant benchmark (Rung 0)")
    print(f"  backend : {backend}")
    print(f"  arb prec: {_ctx.prec if _ctx is not None else 'n/a'} bits")
    print(f"  Y0={Y0}, r={r} (Then r1={THEN_R1}), theta={theta}, C_H={C_H}, eps={eps}")
    print(f"  C_K(r) sharp     = {ball_str(C_K(r))}")
    print(f"  C_K(r) classical = {ball_str(C_K(r, classical=True))}")
    print(f"  C_K poly sharp   = {ball_str(C_K(r, poly=True))}")
    print(f"  using classical_CK={classical_CK}")
    print(f"  trunc: exact r2 for Z[i], n≤{Nmax}")
    print("-" * 72)
    print(f"{'M':>6}  {'tail enclosure':>28}  {'trunc sum':>28}  ratio")
    print("-" * 72)

    for M in M_list:
        encl = tail_majorant(
            M, Y0, r, theta, C_H=C_H, eps=eps, r2_mode="div", classical_CK=classical_CK
        )
        trunc = truncated_majorant_sum(
            M,
            Y0,
            r,
            theta,
            C_H=C_H,
            eps=eps,
            Nmax=Nmax,
            r2_mode="exact_i",
            classical_CK=classical_CK,
        )
        e_mid = ball_mid(encl)
        t_mid = ball_mid(trunc)
        ratio = e_mid / t_mid if t_mid > 0 else float("inf")
        # Upper endpoint of enclosure for < 10^{-30} checks
        e_up = e_mid + ball_rad(encl)
        print(
            f"{M:6d}  {ball_str(encl):>28}  {ball_str(trunc):>28}  {ratio:8.3g}"
        )
        rows.append(
            {
                "M": M,
                "enclosure": encl,
                "truncated": trunc,
                "ratio": ratio,
                "enclosure_upper": e_up,
                "C_K": C_K(r, classical=classical_CK),
            }
        )
    print("=" * 72)
    print(
        "Note: enclosure uses r₂≤6d(n)≤6n (theorem default) + incomplete-gamma;"
        " trunc uses exact r₂ for Z[i].  Ratio ≫ 1 is expected and safe."
    )
    return rows


def kappa_majorant_series(
    M: int,
    Y0: float = 0.8,
    r: float = THEN_R1,
    theta: float = 0.5,
) -> float:
    """
    Single-cusp Hejhal-style condition-number *proxy* (diagnostic only).
    """
    ck = float(ball_mid(C_K(r)))
    diags = []
    for n in range(1, M + 1):
        if r2_exact_gaussian(n) == 0:
            continue
        y = 2.0 * math.pi * math.sqrt(n) * Y0
        k2 = (math.pi / (2.0 * y)) * math.exp(-2.0 * y) * (ck ** 2)
        w = (n ** (2.0 * theta)) * max(k2, 1e-300)
        diags.append(w)
    if not diags:
        return float("nan")
    return max(diags) / min(diags)


# ---------------------------------------------------------------------------
# Self-test / CLI
# ---------------------------------------------------------------------------
def _self_test() -> None:
    print("Self-tests for lemma_K.py (Rung 0)")
    set_arb_precision(DEFAULT_ARB_PREC)
    print(f"  HAS_FLINT = {HAS_FLINT}, HAS_ACB = {HAS_ACB}")
    print(f"  arb prec  = {_ctx.prec if _ctx is not None else 'n/a'} bits")

    # Sharp C_K = 1
    c_sharp = ball_mid(C_K(THEN_R1))
    assert abs(c_sharp - 1.0) < 1e-12, "sharp C_K should be 1"
    print(f"  C_K({THEN_R1}) sharp = {c_sharp}  OK")

    # Classical C_K grows with r
    c1 = ball_mid(C_K(1.0, classical=True))
    c2 = ball_mid(C_K(2.0, classical=True))
    assert c2 > c1 > 1.0, "classical C_K should grow with r"
    print(f"  C_K^{{cl}}(1)={c1:.6g}, C_K^{{cl}}(2)={c2:.6g}  OK")

    # Pointwise majorant positive
    maj = ball_mid(k_bessel_pointwise_majorant(1.0, THEN_R1))
    assert maj > 0
    print(f"  |K| majorant at y=1,r={THEN_R1} → {maj:.6g}  OK")

    # tail_majorant is alias of lemma_K_tail
    a = tail_majorant(100, 0.8, THEN_R1, 0.5)
    b = lemma_K_tail(100, 0.8, THEN_R1, 0.5)
    assert abs(ball_mid(a) - ball_mid(b)) <= 1e-30 * max(1.0, abs(ball_mid(a)))
    print(f"  tail_majorant alias OK  ({ball_str(a)})")

    # Tail decreases in M (sharp C_K); Then parameters
    t100 = ball_mid(tail_majorant(100, 0.8, THEN_R1, 0.5))
    t200 = ball_mid(tail_majorant(200, 0.8, THEN_R1, 0.5))
    t400 = ball_mid(tail_majorant(400, 0.8, THEN_R1, 0.5))
    assert t100 >= t200 >= t400 >= 0, "tail must decrease in M"
    print(f"  tail M=100,200,400: {t100:.6g}, {t200:.6g}, {t400:.6g}  OK")

    # Stopping: enclosure < 10^{-30} at M=100, Y0=0.8, r=Then
    e100 = tail_majorant(100, 0.8, THEN_R1, 0.5)
    e100_up = ball_mid(e100) + ball_rad(e100)
    assert e100_up < 1e-30, f"enclosure upper {e100_up} not < 1e-30"
    print(f"  ε(100,0.8,{THEN_R1},0.5) upper={e100_up:.6g} < 1e-30  OK")

    # Ratio enclosure / trunc(n≤20000) < 10^3
    trunc = truncated_majorant_sum(100, 0.8, THEN_R1, 0.5, Nmax=20000, r2_mode="exact_i")
    ratio = ball_mid(e100) / ball_mid(trunc) if ball_mid(trunc) > 0 else float("inf")
    assert ratio < 1e3, f"ratio {ratio} not < 10^3"
    print(f"  ratio enclosure/trunc(n≤20000) = {ratio:.6g} < 1e3  OK")
    print(f"  trunc uses exact r2(Z[i])  OK")

    # Incomplete gamma sanity: Γ(1,x) = e^{-x}
    g = ball_mid(incomplete_gamma_upper(1.0, 2.0))
    assert abs(g - math.exp(-2.0)) < 1e-6 * (1 if not HAS_FLINT else 10)
    print(f"  Γ(1,2) ≈ {g:.6g} (expect e^{{-2}}={math.exp(-2):.6g})  OK")

    # Monotonicity threshold
    nm = monotonicity_threshold(1.5, 4 * math.pi * 0.8)
    assert nm >= 1
    print(f"  N_mono(α=1.5,c=4π·0.8) = {nm}  OK")

    # r₂ exact checks
    assert r2_exact_gaussian(1) == 4
    assert r2_exact_gaussian(2) == 4
    assert r2_exact_gaussian(3) == 0
    assert r2_exact_gaussian(5) == 8
    print("  r₂ exact checks  OK")

    # C1/C2 finite and positive
    d_i = C1_C2_constants(field="i")
    d_w = C1_C2_constants(field="omega")
    assert ball_mid(d_i["C1"]) > 0 and ball_mid(d_i["C2"]) > 0
    assert ball_mid(d_w["C1"]) > 0 and ball_mid(d_w["C2"]) > 0
    assert ball_mid(d_i["C2"]) > ball_mid(d_w["C2"])  # |T| larger for Z[i]
    print(f"  C1(Z[i])  = {ball_str(d_i['C1'])}")
    print(f"  C2(Z[i])  = {ball_str(d_i['C2'])}")
    print(f"  C1(Z[ω])  = {ball_str(d_w['C1'])}")
    print(f"  C2(Z[ω])  = {ball_str(d_w['C2'])}")
    print("  C1_C2_constants  OK")

    # Arb return type when flint present
    if HAS_FLINT:
        t = tail_majorant(100, 0.8, THEN_R1, 0.5)
        assert _is_arb(t), "tail_majorant must return Arb when flint present"
        assert _is_arb(d_i["C1"]), "C1 must be Arb when flint present"
        print("  Arb enclosure types  OK")

    # Luke UPPER (20 pairs) — never lower bound
    luke = luke_upper_validation(n_samples=20, verbose=True)
    assert luke["pass"], "Luke UPPER validation failed"
    print("  Luke UPPER 20-pair test  OK")

    print("All self-tests passed.")


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Lemma K Arb majorant + D(K) constants")
    p.add_argument("--bench", action="store_true", help="run M=100,200,400 benchmark")
    p.add_argument("--test", action="store_true", help="run self-tests")
    p.add_argument("--luke", action="store_true", help="run Luke UPPER validation only")
    p.add_argument("--constants", action="store_true", help="print C1/C2 tables")
    p.add_argument("--M", type=int, default=None, help="single M evaluation")
    p.add_argument("--Y0", type=float, default=0.8)
    p.add_argument("--r", type=float, default=THEN_R1, help=f"spectral r (default Then {THEN_R1})")
    p.add_argument("--theta", type=float, default=0.5)
    p.add_argument("--C_H", type=float, default=1.0)
    p.add_argument("--eps", type=float, default=0.0)
    p.add_argument(
        "--classical",
        action="store_true",
        help="use classical C_K=exp(π r/2) instead of sharp C_K=1",
    )
    args = p.parse_args(argv)

    set_arb_precision(DEFAULT_ARB_PREC)

    if args.luke:
        luke = luke_upper_validation(n_samples=20, verbose=True)
        return 0 if luke["pass"] else 1

    if args.test or (
        not args.bench and args.M is None and not args.constants
    ):
        _self_test()
        print()

    if args.M is not None:
        val = tail_majorant(
            args.M,
            args.Y0,
            args.r,
            args.theta,
            args.C_H,
            args.eps,
            classical_CK=args.classical,
        )
        print(f"tail_majorant(M={args.M}) = {ball_str(val)}")
        print(f"C_K(r={args.r}) sharp     = {ball_str(C_K(args.r))}")
        print(f"C_K(r={args.r}) classical = {ball_str(C_K(args.r, classical=True))}")

    if args.constants or args.bench:
        print_C1_C2(field="i", Y0=args.Y0, r=args.r, theta=args.theta, C_H=args.C_H, eps=args.eps, classical_CK=args.classical)
        print()
        print_C1_C2(field="omega", Y0=args.Y0, r=args.r, theta=args.theta, C_H=args.C_H, eps=args.eps, classical_CK=args.classical)
        print()

    if args.bench or (args.M is None and not args.test and not args.constants and not args.luke):
        benchmark_majorant(
            M_list=(100, 200, 400),
            Y0=args.Y0,
            r=args.r,
            theta=args.theta,
            C_H=args.C_H,
            eps=args.eps,
            classical_CK=args.classical,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
