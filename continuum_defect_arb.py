#!/usr/bin/env python3
"""Interval continuum automorphy defect for the six-copy trial form.

This verifies a fixed, explicitly recorded coefficient vector.  It covers the
rectangular superset

    [-1/2,1/2]^2 x [1/sqrt(2), 5/4]

of the truncated Picard reference core by Arb boxes and bounds every exact
F5 coset identity F_c(p)-F_{pi_gamma(c)}(gamma p).  Using a superset is
conservative.  The finite Whittaker sum solves (Delta-(1+r^2))f=0 termwise,
so its core PDE residual is analytically zero; no sampled residual is used.

The certificate produced here is for the vector-valued/six-copy continuum
defect.  Applying the single-cusp Theorem D(K) still requires a proved
multi-copy/two-cusp extension with the corresponding constants.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from flint import acb, acb_mat, arb, ctx

from six_copy_hejhal import gaussian_modes


GLUE = {
    "T1": [0, 5, 1, 2, 3, 4],
    "R": [0, 1, 5, 4, 3, 2],
    "TiR": [0, 3, 2, 1, 5, 4],
    "S": [1, 0, 5, 3, 4, 2],
}
GENERATORS = {
    "T1": ((1, 0), (1, 0), (0, 0), (1, 0)),
    "R": ((0, 1), (0, 0), (0, 0), (0, -1)),
    # [[1,i],[0,1]] [[i,0],[0,-i]] = [[i,1],[0,-i]]
    "TiR": ((0, 1), (1, 0), (0, 0), (0, -1)),
    "S": ((0, 0), (-1, 0), (1, 0), (0, 0)),
}


def qcomplex(pair: tuple[int, int]) -> acb:
    return acb(arb(pair[0]), arb(pair[1]))


def interval_hull(a: arb, b: arb) -> arb:
    return a.union(b)


def endpoints(n: int) -> tuple[list[arb], list[arb], list[arb]]:
    n = int(n)
    x0, x1 = -arb(1) / 2, arb(1) / 2
    y0, y1 = arb(1) / arb(2).sqrt(), arb(5) / 4
    xs = [x0 + (x1 - x0) * k / n for k in range(n + 1)]
    ys = [y0 + (y1 - y0) * k / n for k in range(n + 1)]
    return xs, xs, ys


def act(name: str, x1: arb, x2: arb, y: arb) -> tuple[arb, arb, arb]:
    if name == "T1":
        return x1 + 1, x2, y
    if name == "R":
        return -x1, -x2, y
    if name == "TiR":
        return -x1, 1 - x2, y
    if name != "S":
        raise ValueError(name)

    def square_range(t: arb) -> arb:
        # Direct t*t suffers the dependency problem and can include negative
        # values when t straddles zero.  |t| has a tight nonnegative range;
        # squaring its endpoints gives a rigorous range for t^2.
        at = abs(t)
        lo, hi = at.lower(), at.upper()
        return (lo * lo).union(hi * hi)

    denom = square_range(x1) + square_range(x2) + square_range(y)
    if not bool(denom > 0):
        raise ArithmeticError(f"S pullback denominator not positive: {denom}")
    # S(z,y)=(-conj(z)/denom, y/denom).
    return -x1 / denom, x2 / denom, y / denom


def upper(x: arb) -> float:
    return math.nextafter(float(x.upper()), math.inf)


def lower(x: arb) -> float:
    return math.nextafter(float(x.lower()), -math.inf)


class TrialEvaluator:
    def __init__(self, M: int, r_text: str, coefficients: list[dict[str, str]]):
        self.M = int(M)
        self.r = arb(r_text)
        self.order = acb(0, self.r)
        self.modes_inf = gaussian_modes(M)
        self.modes_0 = gaussian_modes(5 * M)
        expected = len(self.modes_inf) + len(self.modes_0)
        if len(coefficients) != expected:
            raise ValueError(f"coefficient count {len(coefficients)} != {expected}")
        self.coeff = [acb(arb(q["real"]), arb(q["imag"])) for q in coefficients]
        self.ni = len(self.modes_inf)
        self._component_cache: dict[tuple, acb] = {}
        self._k_cache: dict[tuple, acb] = {}
        self._gradient_cache: dict[tuple, tuple[arb, arb, arb]] = {}
        self.direct_k_count = 0
        self.mean_value_k_count = 0
        self.pi = arb.pi()

    @staticmethod
    def _key_ball(x: arb) -> str:
        return str(x)

    def _frequency(self, cusp: int, mode: tuple[int, int, int]) -> tuple[arb, arb, arb]:
        a, b, nn = mode
        if cusp == 0:
            return arb(a), arb(b), arb(nn).sqrt()
        # (a+bi)(2+i)/5 = ((2a-b)+i(a+2b))/5 exactly.
        u = arb(2 * a - b) / 5
        v = arb(a + 2 * b) / 5
        mag = (u * u + v * v).sqrt()
        return u, v, mag

    def _kir(self, cusp: int, k: int, mag: arb, y: arb) -> acb:
        key = (cusp, k, self._key_ball(y))
        out = self._k_cache.get(key)
        if out is None:
            arg = 2 * self.pi * mag * y
            out = acb(arg).bessel_k(self.order)
            if out.is_finite() and out.imag.contains(0):
                self.direct_k_count += 1
            else:
                # Arb's direct complex-order interval continuation can become
                # indeterminate on moderately wide positive real balls.  Use
                # the integral-representation derivative inequality
                #   |d/dx K_{ir}(x)| <= K_1(x), x>0,
                # and monotonicity of K_1 to enclose the whole argument ball.
                xlo = arg.lower()
                if not bool(xlo > 0):
                    raise ArithmeticError(f"K argument not provably positive: {arg}")
                mid = arg.mid()
                kmid = acb(mid).bessel_k(self.order)
                k1lo = acb(xlo).bessel_k(1)
                if not kmid.is_finite() or not k1lo.is_finite():
                    raise ArithmeticError(
                        f"mean-value K fallback failed: mid={mid}, lower={xlo}"
                    )
                derivative_upper = abs(k1lo).upper()
                extra = derivative_upper * arg.rad() + abs(kmid.imag).upper()
                out = acb(arb(kmid.real, extra), arb(0, extra))
                if not out.is_finite() or not out.imag.contains(0):
                    raise ArithmeticError(f"invalid mean-value K enclosure: {out}")
                self.mean_value_k_count += 1
            self._k_cache[key] = out
        return out

    def component(self, copy: int, x1: arb, x2: arb, y: arb) -> acb:
        key = (
            int(copy), self._key_ball(x1), self._key_ball(x2), self._key_ball(y)
        )
        cached = self._component_cache.get(key)
        if cached is not None:
            return cached
        cusp = 0 if copy == 0 else 1
        modes = self.modes_inf if cusp == 0 else self.modes_0
        offset = 0 if cusp == 0 else self.ni
        tx = arb(0) if cusp == 0 else arb(copy - 1)
        total = acb(0)
        for k, mode in enumerate(modes):
            u, v, mag = self._frequency(cusp, mode)
            theta = 2 * self.pi * (u * (x1 + tx) + v * x2)
            phase = acb(0, theta).exp()
            total += self.coeff[offset + k] * y * self._kir(cusp, k, mag, y) * phase
        self._component_cache[key] = total
        return total

    def gradient_bounds(self, copy: int, y: arb) -> tuple[arb, arb, arb]:
        """Sup bounds for |d_x1 F|, |d_x2 F|, |d_y F| on a y-ball."""
        key = (int(copy), self._key_ball(y))
        cached = self._gradient_cache.get(key)
        if cached is not None:
            return cached
        cusp = 0 if copy == 0 else 1
        modes = self.modes_inf if cusp == 0 else self.modes_0
        offset = 0 if cusp == 0 else self.ni
        ylo, yhi = y.lower(), y.upper()
        if not bool(ylo > 0):
            raise ArithmeticError(f"gradient y range not positive: {y}")
        bx1 = arb(0)
        bx2 = arb(0)
        by = arb(0)
        for k, mode in enumerate(modes):
            u, v, mag = self._frequency(cusp, mode)
            omega = 2 * self.pi * mag
            xlo = omega * ylo
            # Lemma-K majorant for |K_ir| and exact real-order K1 point
            # enclosure for |d_x K_ir| <= K1.
            kmaj = (self.pi / (2 * xlo)).sqrt() * (-xlo).exp()
            k1 = abs(acb(xlo).bessel_k(1)).upper()
            ca = abs(self.coeff[offset + k]).upper()
            radial = yhi * kmaj
            bx1 += ca * (2 * self.pi * abs(u).upper()) * radial
            bx2 += ca * (2 * self.pi * abs(v).upper()) * radial
            by += ca * (kmaj + yhi * omega * k1)
        out = (bx1, bx2, by)
        self._gradient_cache[key] = out
        return out


def parse_trial(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    row = data["best"]
    if "coefficients" not in row:
        raise ValueError("trial JSON does not contain explicit coefficients")
    return data, row


def normalization_lower(ev: TrialEvaluator, ny: int = 256) -> dict[str, Any]:
    """Exact planar Fourier Gram + interval y-integration on the half torus.

    For T=[-1/2,1/2]x[0,1/2], the x1 integral is delta_{a,a'}, while

      int_0^{1/2} exp(2 pi i d x2) dx2

    equals 1/2 for d=0, 0 for nonzero even d, and i/(pi d) for odd d.
    This directly integrates the recorded (not exactly R-invariant) component
    on the actual area-1/2 orbifold section and preserves all cancellations.
    """
    ya, yb = arb("1.001"), arb("1.249")
    ye = [ya + (yb - ya) * j / ny for j in range(ny + 1)]
    norm2 = arb(0)
    positive_segments = 0
    modes = ev.modes_inf
    groups: dict[int, list[int]] = {}
    for k, (a, _b, _nn) in enumerate(modes):
        groups.setdefault(a, []).append(k)

    # Certify the smallest eigenvalue of each exact half-period Fourier Gram
    # block.  Then v*Gv >= lambda_min ||v||^2 avoids interval cancellation.
    gram_lowers: dict[int, arb] = {}
    gram_records = {}
    for a, ids in groups.items():
        G = acb_mat(len(ids), len(ids))
        for ii, k in enumerate(ids):
            b = modes[k][1]
            for jj, ell in enumerate(ids):
                bb = modes[ell][1]
                d = b - bb
                if d == 0:
                    G[ii, jj] = acb(arb(1) / 2)
                elif d % 2 == 0:
                    G[ii, jj] = acb(0)
                else:
                    G[ii, jj] = acb(0, arb(1) / (ev.pi * d))
        eigs = G.eig(algorithm="rump")
        for eig in eigs:
            if not eig.imag.contains(0):
                raise ArithmeticError(f"Gram eigenvalue not real: {eig}")
        eig_lowers = [eig.real.lower() for eig in eigs]
        lam_hull = eig_lowers[0]
        for candidate in eig_lowers[1:]:
            lam_hull = lam_hull.union(candidate)
        lam = lam_hull.lower()
        if not bool(lam > 0):
            raise ArithmeticError(f"half-torus Gram block not certified PD: a={a}, {eigs}")
        gram_lowers[a] = lam
        gram_records[str(a)] = {"size": len(ids), "lambda_min_lower": lower(lam)}

    for jy in range(ny):
        ybox = ye[jy].union(ye[jy + 1])
        radial = []
        for k, mode in enumerate(modes):
            _u, _v, mag = ev._frequency(0, mode)
            radial.append(ev.coeff[k] * ev._kir(0, k, mag, ybox))
        planar_lower = arb(0)
        for a, ids in groups.items():
            energy_lower = arb(0)
            for k in ids:
                qlo = abs(radial[k]).lower()
                if bool(qlo > 0):
                    energy_lower += qlo * qlo
            planar_lower += gram_lowers[a] * energy_lower.lower()
        integrand = planar_lower / ybox
        lo = integrand.lower()
        if bool(lo > 0):
            norm2 += (ye[jy + 1] - ye[jy]) * lo
            positive_segments += 1
    if not bool(norm2 > 0):
        raise ArithmeticError(f"failed to prove positive L2 normalization: {norm2}")
    norm_lower_arb = norm2.sqrt().lower()
    return {
        "method": "exact half-torus Fourier Gram with interval y-integration",
        "planar_region": "[-1/2,1/2] x [0,1/2] (area 1/2)",
        "y_interval": "[1.001,1.249]",
        "y_segments": ny,
        "positive_segments": positive_segments,
        "gram_blocks": gram_records,
        "norm2_lower": lower(norm2),
        "norm_lower": lower(norm_lower_arb),
        "norm2_ball": str(norm2),
        "_norm_lower_arb": norm_lower_arb,
    }


def certify_level(ev: TrialEvaluator, subdivision: int) -> dict[str, Any]:
    xs, zs, ys = endpoints(subdivision)
    max_upper = -1.0
    max_upper_arb = arb(0)
    worst = None
    per_relation = {name: 0.0 for name in GENERATORS}
    n_boxes = 0
    for i in range(subdivision):
        xb = interval_hull(xs[i], xs[i + 1])
        for j in range(subdivision):
            zb = interval_hull(zs[j], zs[j + 1])
            for k in range(subdivision):
                yb = interval_hull(ys[k], ys[k + 1])
                n_boxes += 1
                for name in GENERATORS:
                    gx, gz, gy = act(name, xb, zb, yb)
                    cx, cz, cy = xb.mid(), zb.mid(), yb.mid()
                    gcx, gcz, gcy = act(name, cx, cz, cy)

                    def deviation(box: arb, center: arb) -> arb:
                        dl = abs(box.lower() - center).upper()
                        du = abs(box.upper() - center).upper()
                        return dl if bool(dl >= du) else du

                    dx = deviation(xb, cx)
                    dz = deviation(zb, cz)
                    dy = deviation(yb, cy)
                    dgx = deviation(gx, gcx)
                    dgz = deviation(gz, gcz)
                    dgy = deviation(gy, gcy)
                    for copy in range(6):
                        other = GLUE[name][copy]
                        defect = ev.component(copy, cx, cz, cy) - ev.component(
                            other, gcx, gcz, gcy
                        )
                        b1, b2, b3 = ev.gradient_bounds(copy, yb)
                        q1, q2, q3 = ev.gradient_bounds(other, gy)
                        bound_ball = (
                            abs(defect).upper()
                            + b1 * dx + b2 * dz + b3 * dy
                            + q1 * dgx + q2 * dgz + q3 * dgy
                        )
                        bound = upper(bound_ball)
                        max_upper_arb = max_upper_arb.union(bound_ball).upper()
                        if bound > per_relation[name]:
                            per_relation[name] = bound
                        if bound > max_upper:
                            max_upper = bound
                            worst = {
                                "relation": name, "copy": copy, "other": other,
                                "box_index": [i, j, k],
                                "center_defect_ball": str(defect),
                                "taylor_bound_ball": str(bound_ball),
                            }
    return {
        "subdivision": subdivision,
        "boxes": n_boxes,
        "relations_times_copies_times_boxes": n_boxes * 24,
        "raw_delta_upper": max_upper,
        "_raw_delta_upper_arb": max_upper_arb,
        "evaluation": "centered first-order Taylor bound with analytic gradient majorants",
        "per_relation_upper": per_relation,
        "worst": worst,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", type=Path, default=Path("six_copy_hejhal_balanced_coeffs.json"))
    ap.add_argument("--subdivisions", default="1,2")
    ap.add_argument("--bits", type=int, default=160)
    ap.add_argument("--norm-y-segments", type=int, default=2048)
    ap.add_argument("--json-out", type=Path, default=Path("continuum_defect_arb_result.json"))
    ns = ap.parse_args()
    ctx.prec = max(ns.bits, 128)
    data, row = parse_trial(ns.trial)
    M = int(data["parameters"]["M"])
    r_text = row["coefficients"] and repr(float(row["r"]))
    ev = TrialEvaluator(M, r_text, row["coefficients"])
    norm = normalization_lower(ev, ny=ns.norm_y_segments)
    norm_lower_arb = norm.pop("_norm_lower_arb")
    levels = []
    for n in [int(x) for x in ns.subdivisions.split(",") if x.strip()]:
        print(f"continuum subdivision={n}", flush=True)
        rec = certify_level(ev, n)
        # Keep the quotient in Arb from the accumulated maximum through final
        # endpoint extraction; no float/repr round-trip is load-bearing.
        raw_up = rec.pop("_raw_delta_upper_arb")
        quotient = raw_up / norm_lower_arb
        rec["normalized_delta_ball"] = str(quotient)
        rec["normalized_delta_upper"] = upper(quotient)
        levels.append(rec)
    monotone = all(
        levels[k + 1]["raw_delta_upper"] <= levels[k]["raw_delta_upper"]
        for k in range(len(levels) - 1)
    )
    out = {
        "status": "Arb continuum box certificate for fixed six-copy finite trial",
        "theorem_DK_admissible": False,
        "theorem_blocker": (
            "theorem_DK.tex is single-cusp PSL2(O_K); a proved two-cusp/six-copy "
            "extension and its constants are still required"
        ),
        "trial": str(ns.trial.resolve()),
        "trial_sha256": hashlib.sha256(ns.trial.read_bytes()).hexdigest(),
        "parameters": {"M_infinity": M, "M_zero": 5 * M, "r": r_text,
                       "precision_bits": int(ctx.prec)},
        "domain_superset": "[-1/2,1/2]^2 x [1/sqrt(2),5/4]",
        "coefficients": "exact decimal point balls from the recorded trial",
        "normalization": norm,
        "levels": levels,
        "bessel_enclosures": {
            "direct_arb": ev.direct_k_count,
            "mean_value_fallback": ev.mean_value_k_count,
            "fallback_bound": "|d K_ir(x)/dx| <= K_1(x)",
        },
        "upper_bounds_monotone": monotone,
        "tau_core": 0,
        "tau_core_justification": (
            "each finite y K_ir(2pi|mu|y) character solves "
            "(Delta-(1+r^2))f=0 analytically"
        ),
        "eta_for_multi_copy_extension_upper": levels[-1]["normalized_delta_upper"],
        "hard_map_changed": False,
    }
    ns.json_out.write_text(json.dumps(out, indent=2, allow_nan=False), encoding="utf-8")
    print(ns.json_out.resolve())
    print(json.dumps({
        "monotone": monotone,
        "raw_delta_upper": levels[-1]["raw_delta_upper"],
        "norm_lower": norm["norm_lower"],
        "normalized_delta_upper": levels[-1]["normalized_delta_upper"],
        "theorem_DK_admissible": False,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
