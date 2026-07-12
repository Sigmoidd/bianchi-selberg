"""Stage 4 — interpolation / pencil constants (field-agnostic core).

Lemma I1 (lower_bound_theory.md): CR interpolant constant
  kappa_sc = sqrt(1/pi^2 + 1/15)   (self-contained, arbitrary tets)

Optional sharper: kappa_1 = sqrt(1/pi^2 + 1/120)  (CP22 / CZZ20).

Window pencil (Lemma P / Thm G1 shape), with |T|-general beta:
  beta(s) = (1-s) / (|T| Y^2)
  kappa_c(s) = 1 / ((1+s) Y^2)

G1 coefficients follow independent_exclusion/cr_prototype.py (Picard),
with beta and |T| adapted for the Eisenstein cusp section.
"""
from __future__ import annotations

import math

import numpy as np

PI = math.pi
KAPPA_SC = math.sqrt(1.0 / PI ** 2 + 1.0 / 15.0)
KAPPA_1 = math.sqrt(1.0 / PI ** 2 + 1.0 / 120.0)  # optional sharper


def beta_s(s: float, area_T: float, Y: float) -> float:
    return (1.0 - s) / (area_T * Y * Y)


def kappa_c(s: float, Y: float) -> float:
    return 1.0 / ((1.0 + s) * Y * Y)


def g1_coeffs(
    s_lo: float,
    s_hi: float,
    area_T: float,
    Y: float,
    cst: dict,
    *,
    rho: float = 55.0,
    nu_star: float = 1.05,
    theta: float = 0.7,
    theta2: float = 0.5,
    thetap: float = 0.5,
    alpha: float = 0.5,
) -> dict:
    """Theorem G1 float coefficients for one s-window (Picard formula, |T|-gen).

    cst must provide: gamma, tau, alpha_h, S_Q, S_M, S_S, V_S
    (from cr_omega.compute_g1_constants).
    """
    s_lo = float(s_lo)
    s_hi = float(s_hi)
    s0 = 0.5 * (s_lo + s_hi)
    lam_p = 1.0 - s_lo ** 2
    beta_p = beta_s(s_lo, area_T, Y)
    kap0 = kappa_c(s0, Y)
    dkap = max(
        abs(kappa_c(s_lo, Y) - kap0),
        abs(kappa_c(s_hi, Y) - kap0),
    )

    omega = 2.0 * rho * (1.0 / theta2 - 1.0) * cst["V_S"]
    c_Q = 1.0 - (omega + nu_star * lam_p) * cst["S_Q"]
    c_S = 1.0 - (omega + nu_star * lam_p) * cst["S_S"]
    lam_t = nu_star * lam_p * (1.0 + cst["S_M"]) + omega * cst["S_M"]
    beta_t = nu_star * beta_p + 2.0 * rho * (1.0 / theta2 - 1.0) * dkap ** 2
    rho_t = rho * (1.0 - theta2)

    sigma_h = cst["alpha_h"] + kap0 * cst["tau"]
    c_e = c_Q - rho_t * (1.0 / theta - 1.0) * sigma_h ** 2
    d_e = (
        lam_t * (1.0 + 1.0 / alpha) * cst["gamma"] ** 2
        + beta_t * (1.0 + 1.0 / thetap) * cst["tau"] ** 2
    )

    return dict(
        s_lo=s_lo,
        s_hi=s_hi,
        lam_p=lam_p,
        kap0=kap0,
        beta_p=beta_p,
        c_Q=c_Q,
        c_S=c_S,
        lam_t=lam_t,
        beta_t=beta_t,
        rho_t=rho_t,
        sigma_h=sigma_h,
        c_e=c_e,
        d_e=d_e,
        theta=theta,
        theta2=theta2,
        thetap=thetap,
        alpha=alpha,
        rho=rho,
        nu_star=nu_star,
    )


# Back-compat alias used by older call sites
def window_coeffs(s_lo, s_hi, area_T, Y, **kwargs):
    """Legacy crude coefficients; prefer g1_coeffs with real mesh constants."""
    gamma = kwargs.pop("gamma", 0.2)
    sigma_h = kwargs.pop("sigma_h", 0.1)
    S_Q = kwargs.pop("S_Q", 0.05)
    S_M = kwargs.pop("S_M", 0.4)
    V_S = kwargs.pop("V_S", 0.02)
    S_S = kwargs.pop("S_S", S_Q)
    tau = kwargs.pop("tau", 0.05)
    alpha_h = kwargs.pop("alpha_h", sigma_h)
    cst = dict(
        gamma=gamma, tau=tau, alpha_h=alpha_h,
        S_Q=S_Q, S_M=S_M, S_S=S_S, V_S=V_S,
    )
    # map old param names
    rho = kwargs.pop("rho_t", kwargs.pop("rho", 55.0))
    # old code passed rho_t as the free parameter; interpret as rho
    if "nu" in kwargs:
        kwargs["nu_star"] = kwargs.pop("nu")
    return g1_coeffs(s_lo, s_hi, area_T, Y, cst, rho=rho, **{
        k: v for k, v in kwargs.items()
        if k in ("nu_star", "theta", "theta2", "thetap", "alpha")
    })


def tet_h_and_gamma(X, tets, kappa=KAPPA_SC):
    """Per-tet diameter h_T and crude gamma ~ kappa * h / sqrt(w) proxy.

    Prefer cr_omega.compute_g1_constants for the proper weighted gamma.
    """
    nt = len(tets)
    hT = np.empty(nt)
    gamma = 0.0
    for e, tet in enumerate(tets):
        P = X[list(tet)]
        d2 = 0.0
        for a in range(4):
            for b in range(a + 1, 4):
                d2 = max(d2, float(np.sum((P[a] - P[b]) ** 2)))
        hT[e] = math.sqrt(d2)
        ymin = max(float(P[:, 2].min()), 1e-12)
        gamma = max(gamma, kappa * hT[e] / math.sqrt(ymin))
    return hT, gamma
