#!/usr/bin/env python3
"""
Diagnostic: log κ_proxy(V(r)) vs log M for the single-cusp analytic majorant.

This is NOT a two-cusp conditioning estimate and does NOT use FEM/gluing.
It plots the pure Fourier-mode diagonal majorant ratio

    κ_proxy(M) = max_{n≤M, r₂(n)>0} w(n)  /  min_{n≤M, r₂(n)>0} w(n),

where w(n) = n^{2θ} · [√(π/(2y)) e^{-y} C_K(r)]² , y = 2π √n Y0.

Polynomial growth in M is acceptable for interval Hejhal; exponential growth
would flag that interval widths will explode.

Outputs (written next to this script):
  - kappa_vs_M.csv
  - kappa_vs_M.png  (if matplotlib is available)
"""
from __future__ import annotations

import csv
import math
import os
import sys
from typing import List, Tuple

# Local import
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lemma_K import C_K, ball_mid, kappa_majorant_series, r2_exact_gaussian  # noqa: E402


def kappa_algebraic(
    M: int,
    Y0: float = 0.8,
    r: float = 6.6,
    theta: float = 0.5,
) -> float:
    """
    De-exponentialized diagonal ratio: strip the universal e^{-2y} factor so
    only the algebraic weight n^{2θ} / y remains.  Polynomial growth here is
    the optimistic signal; the raw kappa_majorant_series includes K-decay and
    is always roughly exp(Θ(√M)).
    """
    ck = float(ball_mid(C_K(r)))
    diags = []
    for n in range(1, M + 1):
        if r2_exact_gaussian(n) == 0:
            continue
        y = 2.0 * math.pi * math.sqrt(n) * Y0
        # |K|² majorant without e^{-2y}: (π/(2y)) C_K²
        k2_alg = (math.pi / (2.0 * y)) * (ck ** 2)
        w = (n ** (2.0 * theta)) * k2_alg
        diags.append(w)
    if not diags:
        return float("nan")
    return max(diags) / min(diags)


def kappa_and_series(
    M_values: List[int],
    Y0: float = 0.8,
    r: float = 6.6,
    theta: float = 0.5,
) -> List[Tuple[int, float, float, float]]:
    """Return list of (M, kappa_raw, log10 kappa_raw, kappa_alg)."""
    rows = []
    for M in M_values:
        kap = kappa_majorant_series(M, Y0=Y0, r=r, theta=theta)
        kalg = kappa_algebraic(M, Y0=Y0, r=r, theta=theta)
        logk = math.log10(kap) if kap > 0 and math.isfinite(kap) else float("nan")
        rows.append((M, kap, logk, kalg))
    return rows


def fit_loglog(rows: List[Tuple[int, float, float]]) -> Tuple[float, float]:
    """Least-squares log κ ~ a + b log M  (returns a, b)."""
    xs, ys = [], []
    for M, kap, _ in rows:
        if kap > 0 and math.isfinite(kap):
            xs.append(math.log(float(M)))
            ys.append(math.log(kap))
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


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    Y0, r, theta = 0.8, 6.6, 0.5
    M_values = [50, 100, 150, 200, 300, 400, 600, 800]

    print("Single-cusp κ_proxy diagnostic")
    print(f"  Y0={Y0}, r={r}, theta={theta}")
    print(f"  C_K(r) = {ball_mid(C_K(r)):.6g}")
    rows = kappa_and_series(M_values, Y0=Y0, r=r, theta=theta)

    csv_path = os.path.join(here, "kappa_vs_M.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            ["M", "kappa_raw", "log10_kappa_raw", "kappa_algebraic", "log_M"]
        )
        for M, kap, logk, kalg in rows:
            w.writerow(
                [
                    M,
                    f"{kap:.16e}",
                    f"{logk:.16e}",
                    f"{kalg:.16e}",
                    f"{math.log(M):.16e}",
                ]
            )
            print(
                f"  M={M:4d}  κ_raw={kap:.6e}  log10κ={logk:.4f}  "
                f"κ_alg={kalg:.6e}"
            )

    # Fit algebraic proxy (the meaningful polynomial-vs-exponential test)
    rows_alg = [(M, kalg, math.log10(kalg) if kalg > 0 else float("nan"))
                for M, _, _, kalg in rows]
    a, b = fit_loglog(rows_alg)
    print(f"  log-log fit (algebraic): log κ_alg ≈ {a:.4f} + {b:.4f} log M")
    print(f"  effective power κ_alg ~ M^{b:.4f}")
    if math.isfinite(b):
        if b < 4:
            print(
                "  Verdict (algebraic): polynomial growth — acceptable skeleton."
            )
        else:
            print("  Verdict (algebraic): steep — inspect weights.")
    print(
        "  Note: κ_raw includes K-Bessel e^{-2y} decay and grows ~exp(c√M);"
        " increase Y0 or moderate M for interval Hejhal."
    )
    print(f"  Wrote {csv_path}")

    # Optional plot
    png_path = os.path.join(here, "kappa_vs_M.png")
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        Ms = [row[0] for row in rows]
        kaps = [row[1] for row in rows]
        kalgs = [row[3] for row in rows]
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

        ax = axes[0]
        ax.semilogy(Ms, kaps, "o-", color="C0", label=r"$\kappa_{\mathrm{raw}}$")
        ax.set_xlabel(r"$M$")
        ax.set_ylabel(r"$\kappa_{\mathrm{raw}}$ (includes $e^{-2y}$)")
        ax.set_title("Raw diagonal majorant ratio\n(~exp in $\\sqrt{M}$)")
        ax.grid(True, which="both", ls=":", alpha=0.5)
        ax.legend()

        ax = axes[1]
        ax.loglog(Ms, kalgs, "s-", color="C1", label=r"$\kappa_{\mathrm{alg}}$")
        M0, k0 = Ms[0], kalgs[0]
        ax.loglog(
            Ms,
            [k0 * (M / M0) ** 2 for M in Ms],
            "--",
            color="gray",
            alpha=0.7,
            label=r"$\propto M^{2}$",
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
        ax.set_ylabel(r"$\kappa_{\mathrm{alg}}$ (no exponential)")
        ax.set_title(
            rf"Algebraic weight ratio"
            f"\n$Y_0={Y0}$, $r={r}$, $\\theta={theta}$"
        )
        ax.grid(True, which="both", ls=":", alpha=0.5)
        ax.legend()

        fig.tight_layout()
        fig.savefig(png_path, dpi=140)
        plt.close(fig)
        print(f"  Wrote {png_path}")
    except Exception as exc:
        print(f"  matplotlib plot skipped ({exc}); CSV is the archival output.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
