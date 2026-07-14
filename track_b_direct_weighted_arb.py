#!/usr/bin/env python3
r"""Certified cusp-blend part of the direct Track-B gluing residual.

This program keeps the pointwise correlation in

    || B0*d0 + B1*d1 ||_2

for the exact Track-B cutoff chi_B on 1.01 <= y <= 1.20.  It also evaluates
the still sharper cancellation-preserving commutator

    (Delta chi_B)(W_B-F) - 2 <grad chi_B, grad(W_B-F)>.

Here F is the recorded six-copy finite Whittaker field and W_B is the exact
old, (-z)-even, quarter-turn-odd cusp witness repeated in all six fibers.
The embedded cusp makes this part independent of the unresolved curved
Humbert collar.  It is only one summand of the total gluing residual; the
program therefore never sets rung4_certified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from flint import acb, arb, ctx

from continuum_defect_arb import lower, parse_trial, upper
from track_b_overlap_arb import (
    DifferentialTrial,
    action_jacobian,
    parse_matrix,
    transpose_pullback,
)
from track_b_projected_mass_arb import integrate_projected_norm, projected_coefficients


ROOT = Path(__file__).resolve().parent
Y_SUPPORT = arb("1.01")
Y_PLATEAU = arb("1.20")
WIDTH = Y_PLATEAU - Y_SUPPORT
TARGET_R = arb("0.0101903405004245")


def interval(a: arb, b: arb) -> arb:
    return a.union(b)


def square_nonnegative(x: arb) -> arb:
    ax = abs(x)
    return (ax.lower() * ax.lower()).union(ax.upper() * ax.upper())


def cutoff_derivatives(y: arb) -> tuple[arb, arb, arb]:
    """Return chi', chi'', and hyperbolic Delta chi on a transition box."""
    # Endpoint balls have outward radii, so compare for disjointness rather
    # than demanding their represented sets lie strictly inside point balls.
    if bool(y.upper() < Y_SUPPORT.lower()) or bool(y.lower() > Y_PLATEAU.upper()):
        raise ValueError(f"cutoff box outside transition: {y}")
    t = (y - Y_SUPPORT) / WIDTH
    one_minus = 1 - t
    first = 30 * t * t * one_minus * one_minus / WIDTH
    second = (60 * t - 180 * t * t + 120 * t * t * t) / (WIDTH * WIDTH)
    lap = -square_nonnegative(y) * second + y * first
    return first, second, lap


def cutoff_third_lap_prime(y: arb) -> tuple[arb, arb]:
    t = (y - Y_SUPPORT) / WIDTH
    third = (60 - 360 * t + 360 * t * t) / (WIDTH ** 3)
    first, second, _lap = cutoff_derivatives(y)
    lap_prime = first - y * second - square_nonnegative(y) * third
    return third, lap_prime


def cutoff_fourth(y: arb) -> arb:
    t = (y - Y_SUPPORT) / WIDTH
    return (-360 + 720 * t) / (WIDTH ** 4)


def witness_value_gradient(
    ev: DifferentialTrial,
    odd: list[acb],
    x1: arb,
    x2: arb,
    y: arb,
) -> tuple[acb, list[acb]]:
    value = acb(0)
    grad = [acb(0), acb(0), acb(0)]
    for k, mode in enumerate(ev.modes_inf):
        u, v, mag = ev._frequency(0, mode)
        omega = 2 * ev.pi * mag
        arg = omega * y
        kval = ev._kir(0, k, mag, y)
        kprime = -(ev._kir_order(ev.order - 1, arg)
                   + ev._kir_order(ev.order + 1, arg)) / 2
        phase = acb(0, 2 * ev.pi * (u * x1 + v * x2)).exp()
        term = odd[k] * phase
        radial = y * kval
        value += term * radial
        grad[0] += term * acb(0, 2 * ev.pi * u) * radial
        grad[1] += term * acb(0, 2 * ev.pi * v) * radial
        grad[2] += term * (kval + y * omega * kprime)
    return value, grad


def box_bounds(
    ev: DifferentialTrial,
    odd: list[acb],
    x1: arb,
    x2: arb,
    y: arb,
) -> dict[str, arb]:
    wv, wg = witness_value_gradient(ev, odd, x1, x2, y)
    fv, fg = ev.vector_value_gradient(x1, x2, y)
    first, _second, lap = cutoff_derivatives(y)

    d0sq = arb(0)
    d1sq = arb(0)
    actual2 = arb(0)
    y2 = square_nonnegative(y)
    for c in range(6):
        dv = wv - fv[c]
        dg = [wg[k] - fg[c][k] for k in range(3)]
        d0sq += abs(dv).upper() ** 2
        for k in range(3):
            d1sq += y2.upper() * abs(dg[k]).upper() ** 2
        # Only the y derivative enters <grad chi,grad(W-F)>.
        residual = lap * dv - 2 * y2 * first * dg[2]
        actual2 += abs(residual).upper() ** 2

    d0 = d0sq.sqrt().upper()
    d1 = d1sq.sqrt().upper()
    value_part = abs(lap).upper() * d0
    gradient_part = 2 * y.upper() * abs(first).upper() * d1
    requested = (value_part + gradient_part).upper()
    actual = actual2.sqrt().upper()
    # This should follow analytically; retain a machine-auditable check.
    if bool(actual > requested):
        raise ArithmeticError(
            f"cancellation-preserving residual exceeds triangle bound: {actual} > {requested}"
        )
    return {
        "d0": d0,
        "d1": d1,
        "value_part": value_part,
        "gradient_part": gradient_part,
        "requested": requested,
        "actual": actual,
    }


def edges(a: arb, b: arb, n: int) -> list[arb]:
    return [a + (b - a) * k / n for k in range(n + 1)]


def hyperbolic_box_volume(
    xa: arb, xb: arb, ta: arb, tb: arb, ya: arb, yb: arb
) -> arb:
    # Integral y^-3 dy = (ya^-2-yb^-2)/2.
    return (xb - xa) * (tb - ta) * ((1 / (ya * ya) - 1 / (yb * yb)) / 2)


def integrate_level(
    ev: DifferentialTrial,
    odd: list[acb],
    nx1: int,
    nx2: int,
    ny: int,
) -> dict[str, Any]:
    xe = edges(-arb(1) / 2, arb(1) / 2, nx1)
    te = edges(arb(0), arb(1) / 2, nx2)
    ye = edges(Y_SUPPORT, Y_PLATEAU, ny)
    sums = {k: arb(0) for k in ("requested", "actual", "value_part", "gradient_part")}
    maxima = {k: arb(0) for k in ("d0", "d1", "requested", "actual")}
    worst_ratio = 0.0
    count = 0
    for i in range(nx1):
        for j in range(nx2):
            for k in range(ny):
                xb = interval(xe[i], xe[i + 1])
                tb = interval(te[j], te[j + 1])
                yb = interval(ye[k], ye[k + 1])
                q = box_bounds(ev, odd, xb, tb, yb)
                vol = hyperbolic_box_volume(
                    xe[i], xe[i + 1], te[j], te[j + 1], ye[k], ye[k + 1]
                ).upper()
                for name in sums:
                    sums[name] += vol * q[name] * q[name]
                for name in maxima:
                    maxima[name] = maxima[name].union(q[name]).upper()
                req = float(q["requested"].upper())
                act = float(q["actual"].upper())
                if req > 0:
                    worst_ratio = max(worst_ratio, act / req)
                count += 1

    norms = {name: sums[name].sqrt().upper() for name in sums}
    return {
        "subdivision": [nx1, nx2, ny],
        "boxes": count,
        "requested_pointwise_weighted_norm_upper": upper(norms["requested"]),
        "cancellation_preserving_commutator_norm_upper": upper(norms["actual"]),
        "value_contribution_norm_upper": upper(norms["value_part"]),
        "gradient_contribution_norm_upper": upper(norms["gradient_part"]),
        "maxima": {name: upper(value) for name, value in maxima.items()},
        "actual_to_triangle_box_upper_ratio_max": math.nextafter(worst_ratio, math.inf),
        "arb_balls": {name: str(value) for name, value in norms.items()},
    }


def difference_fourier_data(
    ev: DifferentialTrial, odd: list[acb]
) -> tuple[list[tuple[int, int]], list[list[acb]], dict[tuple[int, int], tuple[int, int, arb]]]:
    """Return exact frequency/5 coefficients for W_B-F in all six fibers."""
    coeffs: dict[tuple[int, int], list[acb]] = {}
    source: dict[tuple[int, int], tuple[int, int, arb]] = {}

    def row(key: tuple[int, int]) -> list[acb]:
        return coeffs.setdefault(key, [acb(0) for _ in range(6)])

    # W_B is an infinity-lattice scalar repeated in every induced component.
    for k, (a, b, _nn) in enumerate(ev.modes_inf):
        key = (5 * a, 5 * b)
        _u, _v, mag = ev._frequency(0, ev.modes_inf[k])
        source.setdefault(key, (0, k, mag))
        for c in range(6):
            row(key)[c] += odd[k]
        row(key)[0] -= ev.coeff[k]

    # Components 1,...,5 use the zero-cusp lattice and translated x1 phase.
    for k, (a, b, _nn) in enumerate(ev.modes_0):
        key = (2 * a - b, a + 2 * b)
        u, _v, mag = ev._frequency(1, ev.modes_0[k])
        source.setdefault(key, (1, k, mag))
        for c in range(1, 6):
            phase = acb(0, 2 * ev.pi * u * (c - 1)).exp()
            row(key)[c] -= ev.coeff[ev.ni + k] * phase

    # Removing exact zero rows is optional but materially reduces Gram size.
    keys = []
    rows = []
    for key in sorted(coeffs):
        if any(not (q.real.is_zero() and q.imag.is_zero()) for q in coeffs[key]):
            keys.append(key)
            rows.append(coeffs[key])
    return keys, rows, source


def exponential_integral(
    diff_num: int, denominator: int, a: arb, b: arb
) -> acb:
    """Integral of exp(2*pi*i*(diff_num/denominator)*x) from a to b."""
    if diff_num == 0:
        return acb(b - a)
    q = arb(diff_num) / denominator
    eb = acb(0, 2 * arb.pi() * q * b).exp()
    ea = acb(0, 2 * arb.pi() * q * a).exp()
    return (eb - ea) / acb(0, 2 * arb.pi() * q)


def planar_gram(
    keys: list[tuple[int, int]],
    rows: list[list[acb]],
    x1_bounds: tuple[arb, arb] = (-arb(1) / 2, arb(1) / 2),
    x2_bounds: tuple[arb, arb] = (arb(0), arb(1) / 2),
) -> list[list[acb]]:
    """Gram including the six-fiber coefficient inner product."""
    n = len(keys)
    out = [[acb(0) for _ in range(n)] for _ in range(n)]
    for i, (ui, vi) in enumerate(keys):
        for j, (uj, vj) in enumerate(keys):
            gx = exponential_integral(uj - ui, 5, *x1_bounds)
            gt = exponential_integral(vj - vi, 5, *x2_bounds)
            fiber = sum(
                (rows[i][c].conjugate() * rows[j][c] for c in range(6)),
                acb(0),
            )
            out[i][j] = gx * gt * fiber
    return out


def radial_vectors(
    ev: DifferentialTrial,
    keys: list[tuple[int, int]],
    source: dict[tuple[int, int], tuple[int, int, arb]],
    y: arb,
) -> tuple[list[acb], list[acb]]:
    values, derivatives = [], []
    for key in keys:
        cusp, k, mag = source[key]
        omega = 2 * ev.pi * mag
        arg = omega * y
        kval = ev._kir(cusp, k, mag, y)
        kprime = -(ev._kir_order(ev.order - 1, arg)
                   + ev._kir_order(ev.order + 1, arg)) / 2
        values.append(y * kval)
        derivatives.append(kval + y * omega * kprime)
    return values, derivatives


def shell_grams(
    ev: DifferentialTrial,
    keys: list[tuple[int, int]],
    H: list[list[acb]],
    source: dict[tuple[int, int], tuple[int, int, arb]],
) -> tuple[list[int], dict[str, list[list[acb]]], dict[int, tuple[int, int, arb]]]:
    """Collapse all exactly equal radial Bessel factors before intervals."""
    shells = sorted({u * u + v * v for u, v in keys})
    pos = {nn: k for k, nn in enumerate(shells)}
    n = len(shells)
    mats = {name: [[acb(0) for _ in range(n)] for _ in range(n)]
            for name in ("value", "x1", "x2")}
    shell_source: dict[int, tuple[int, int, arb]] = {}
    for key in keys:
        shell_source.setdefault(key[0] * key[0] + key[1] * key[1], source[key])
    for i, (ui, vi) in enumerate(keys):
        ii = pos[ui * ui + vi * vi]
        fu = 2 * ev.pi * arb(ui) / 5
        fv = 2 * ev.pi * arb(vi) / 5
        for j, (uj, vj) in enumerate(keys):
            jj = pos[uj * uj + vj * vj]
            gu = 2 * ev.pi * arb(uj) / 5
            gv = 2 * ev.pi * arb(vj) / 5
            mats["value"][ii][jj] += H[i][j]
            mats["x1"][ii][jj] += fu * gu * H[i][j]
            mats["x2"][ii][jj] += fv * gv * H[i][j]
    return shells, mats, shell_source


def shell_radial_vectors(
    ev: DifferentialTrial,
    shells: list[int],
    shell_source: dict[int, tuple[int, int, arb]],
    y: arb,
) -> tuple[list[acb], list[acb]]:
    values, derivatives = [], []
    for nn in shells:
        cusp, k, mag = shell_source[nn]
        omega = 2 * ev.pi * mag
        arg = omega * y
        kval = ev._kir(cusp, k, mag, y)
        kprime = -(ev._kir_order(ev.order - 1, arg)
                   + ev._kir_order(ev.order + 1, arg)) / 2
        values.append(y * kval)
        derivatives.append(kval + y * omega * kprime)
    return values, derivatives


def quadratic_energy(H: list[list[acb]], vector: list[acb], label: str) -> arb:
    ball = quadratic_energy_ball(H, vector, label)
    if bool(ball.upper() < 0):
        raise ArithmeticError(f"{label} Fourier energy certified negative: {ball}")
    return max(ball.upper(), arb(0))


def quadratic_energy_ball(
    H: list[list[acb]], vector: list[acb], label: str
) -> arb:
    total = acb(0)
    n = len(vector)
    for i in range(n):
        vi = vector[i].conjugate()
        for j in range(n):
            total += vi * H[i][j] * vector[j]
    if not total.imag.contains(0):
        raise ArithmeticError(f"{label} Fourier energy not certified real: {total}")
    return total.real


def integrate_fourier_level(
    ev: DifferentialTrial,
    frequency_count: int,
    shells: list[int],
    mats: dict[str, list[list[acb]]],
    shell_source: dict[int, tuple[int, int, arb]],
    ny: int,
) -> dict[str, Any]:
    ye = edges(Y_SUPPORT, Y_PLATEAU, ny)
    value2 = arb(0)
    gradient2 = arb(0)
    actual2 = arb(0)
    max_e0 = arb(0)
    max_e1 = arb(0)
    for k in range(ny):
        y = interval(ye[k], ye[k + 1])
        dy = ye[k + 1] - ye[k]
        first, _second, lap = cutoff_derivatives(y)
        radial, radial_y = shell_radial_vectors(ev, shells, shell_source, y)
        e0 = quadratic_energy(mats["value"], radial, "value")
        ex = quadratic_energy(mats["x1"], radial, "x1 gradient")
        et = quadratic_energy(mats["x2"], radial, "x2 gradient")
        ey = quadratic_energy(mats["value"], radial_y, "y gradient")
        egrad = ex + et + ey
        y2 = square_nonnegative(y)
        residual_radial = [
            lap * radial[i] - 2 * y2 * first * radial_y[i]
            for i in range(len(shells))
        ]
        eres = quadratic_energy(mats["value"], residual_radial, "cutoff commutator")

        # Planar energies are already integrated over the half torus.
        # Hyperbolic measure contributes y^-3 dy.
        ylo = y.lower()
        if not bool(ylo > 0):
            raise ArithmeticError(f"nonpositive y in Fourier integration: {y}")
        value_integrand = abs(lap).upper() ** 2 * e0 / (ylo ** 3)
        # d1^2 = y^2 * egrad and 2|grad chi| = 2*y*|chi'|.
        gradient_integrand = (
            4 * y.upper() ** 4 * abs(first).upper() ** 2 * egrad / (ylo ** 3)
        )
        actual_integrand = eres / (ylo ** 3)
        value2 += dy * value_integrand
        gradient2 += dy * gradient_integrand
        actual2 += dy * actual_integrand
        max_e0 = max_e0.union(e0).upper()
        max_e1 = max_e1.union(y2.upper() * egrad).upper()

    value_norm = value2.sqrt().upper()
    gradient_norm = gradient2.sqrt().upper()
    requested = (value_norm + gradient_norm).upper()  # Minkowski in L2.
    actual = actual2.sqrt().upper()
    return {
        "y_segments": ny,
        "frequency_count": frequency_count,
        "radial_shell_count": len(shells),
        "method": "exact rational-frequency half-torus Gram, exact shell collapse, Arb y boxes",
        "requested_pointwise_weighted_norm_upper": upper(requested),
        "cancellation_preserving_commutator_norm_upper": upper(actual),
        "value_contribution_norm_upper": upper(value_norm),
        "gradient_contribution_norm_upper": upper(gradient_norm),
        "max_planar_d0_squared_upper": upper(max_e0),
        "max_planar_d1_squared_upper": upper(max_e1),
        "arb_balls": {
            "requested": str(requested),
            "actual": str(actual),
            "value_part": str(value_norm),
            "gradient_part": str(gradient_norm),
        },
    }


def midpoint_diagnostic(
    ev: DifferentialTrial,
    shells: list[int],
    mats: dict[str, list[list[acb]]],
    shell_source: dict[int, tuple[int, int, arb]],
    ny: int,
) -> dict[str, Any]:
    """High-accuracy Arb point evaluations with non-certified midpoint quadrature."""
    dy = (Y_PLATEAU - Y_SUPPORT) / ny
    value2 = arb(0)
    gradient2 = arb(0)
    actual2 = arb(0)
    for k in range(ny):
        y = Y_SUPPORT + dy * (arb(k) + arb(1) / 2)
        first, _second, lap = cutoff_derivatives(y)
        radial, radial_y = shell_radial_vectors(ev, shells, shell_source, y)
        e0 = quadratic_energy(mats["value"], radial, "midpoint value")
        ex = quadratic_energy(mats["x1"], radial, "midpoint x1")
        et = quadratic_energy(mats["x2"], radial, "midpoint x2")
        ey = quadratic_energy(mats["value"], radial_y, "midpoint y")
        egrad = ex + et + ey
        y2 = y * y
        residual_radial = [
            lap * radial[i] - 2 * y2 * first * radial_y[i]
            for i in range(len(shells))
        ]
        eres = quadratic_energy(mats["value"], residual_radial, "midpoint residual")
        value2 += dy * lap * lap * e0 / (y ** 3)
        gradient2 += dy * 4 * (y ** 4) * first * first * egrad / (y ** 3)
        actual2 += dy * eres / (y ** 3)
    value = value2.sqrt()
    gradient = gradient2.sqrt()
    return {
        "load_bearing": False,
        "method": "Arb point evaluations with midpoint quadrature (convergence diagnostic only)",
        "y_segments": ny,
        "requested_estimate": upper(value + gradient),
        "commutator_estimate": upper(actual2.sqrt()),
        "value_estimate": upper(value),
        "gradient_estimate": upper(gradient),
    }


def cutoff_mass_sweep(
    ev: DifferentialTrial,
    odd: list[acb],
    shells: list[int],
    mats: dict[str, list[list[acb]]],
    shell_source: dict[int, tuple[int, int, arb]],
    plateau_texts: list[str],
    diagnostic_segments: int,
    mass_segments: int,
) -> list[dict[str, Any]]:
    """Diagnostic residual versus certified plateau-mass budget."""
    global Y_PLATEAU, WIDTH
    saved_plateau, saved_width = Y_PLATEAU, WIDTH
    records = []
    try:
        for text in plateau_texts:
            print(f"cutoff/mass sweep plateau={text}", flush=True)
            Y_PLATEAU = arb(text)
            WIDTH = Y_PLATEAU - Y_SUPPORT
            if not bool(WIDTH > 0):
                raise ValueError(f"plateau must exceed support: {text}")
            diag = midpoint_diagnostic(
                ev, shells, mats, shell_source, diagnostic_segments
            )
            y_max = Y_PLATEAU + arb("0.25")
            try:
                mass = integrate_projected_norm(
                    ev,
                    odd,
                    mass_segments,
                    str(Y_PLATEAU.mid()),
                    str(y_max.mid()),
                )
            except ArithmeticError as exc:
                records.append({
                    "support_y": str(Y_SUPPORT.mid()),
                    "plateau_y": text,
                    "witness_slab_y": [str(Y_PLATEAU.mid()), str(y_max.mid())],
                    "certified_mu_lower": None,
                    "midpoint_commutator_estimate": diag["commutator_estimate"],
                    "mass_failure": str(exc),
                    "residual_load_bearing": False,
                    "mass_load_bearing": False,
                })
                continue
            mu = arb(repr(mass["certified_witness_norm_lower"]))
            budget = mu * arb("0.05")
            comm = arb(repr(diag["commutator_estimate"]))
            records.append({
                "support_y": str(Y_SUPPORT.mid()),
                "plateau_y": text,
                "witness_slab_y": [str(Y_PLATEAU.mid()), str(y_max.mid())],
                "certified_mu_lower": lower(mu),
                "certified_width_budget_0.05_mu": lower(budget),
                "midpoint_commutator_estimate": diag["commutator_estimate"],
                "estimated_budget_ratio": math.nextafter(
                    float(comm.upper() / budget.lower()), math.inf
                ),
                "residual_load_bearing": False,
                "mass_load_bearing": True,
            })
    finally:
        Y_PLATEAU, WIDTH = saved_plateau, saved_width
    return records


def projected_core_face_diagnostic(
    ev: DifferentialTrial, odd: list[acb], n: int = 9
) -> dict[str, Any]:
    """Arb point diagnostics for using W_B as every scalar core candidate."""
    matrices = {
        "T1": [[[1, 0], [1, 0]], [[0, 0], [1, 0]]],
        "R": [[[0, 1], [0, 0]], [[0, 0], [0, -1]]],
        "TiR": [[[0, 1], [1, 0]], [[0, 0], [0, -1]]],
        "S": [[[0, 0], [-1, 0]], [[1, 0], [0, 0]]],
    }
    records = []
    sqrt6 = arb(6).sqrt()
    for name, raw in matrices.items():
        mat = parse_matrix(raw)
        vmax, gmax = arb(0), arb(0)
        count = 0
        for i in range(n):
            x1 = arb("-0.45") + arb("0.90") * (arb(i) + arb("0.37")) / n
            for j in range(n):
                x2 = arb("0.05") + arb("0.40") * (arb(j) + arb("0.37")) / n
                if name == "S":
                    y2 = 1 - x1 * x1 - x2 * x2
                    if not bool(y2 > 0):
                        continue
                    y = y2.sqrt()
                else:
                    y = arb("0.85") + arb("0.25") * (arb((i + 2 * j) % n) + arb("0.37")) / n
                p = (x1, x2, y)
                gp, jac = action_jacobian(mat, *p)
                lv, lg = witness_value_gradient(ev, odd, *p)
                rv, rg = witness_value_gradient(ev, odd, *gp)
                pulled = transpose_pullback(jac, rg)
                vd = sqrt6 * abs(lv - rv).upper()
                gd2 = sum(
                    (abs(lg[k] - pulled[k]).upper() ** 2 for k in range(3)),
                    arb(0),
                )
                gd = sqrt6 * y.upper() * gd2.sqrt().upper()
                vmax = vmax.union(vd).upper()
                gmax = gmax.union(gd).upper()
                count += 1
        records.append({
            "relation": name,
            "points": count,
            "six_vector_value_defect_max_upper": upper(vmax),
            "six_vector_gradient_defect_max_upper": upper(gmax),
        })
    return {
        "load_bearing": False,
        "method": "Arb point evaluations; deterministic face/interior sample, no continuum fill-in",
        "candidate": "W_B repeated in all six fibers on cusp and core charts",
        "cusp_blend_defect": 0,
        "relations": records,
    }


def projected_floor_midpoint_diagnostic(
    ev: DifferentialTrial,
    odd: list[acb],
    nx1: int = 8,
    nx2: int = 4,
    ns: int = 8,
    epsilon: arb = arb("0.1"),
) -> dict[str, Any]:
    """Midpoint diagnostic for the sole nonzero projected-core S transition."""
    mat = parse_matrix([[[0, 0], [-1, 0]], [[1, 0], [0, 0]]])
    dx = arb(1) / nx1
    dt = (arb(1) / 2) / nx2
    ds = epsilon / ns
    actual2 = arb(0)
    requested2 = arb(0)
    value2 = arb(0)
    gradient2 = arb(0)
    sqrt6 = arb(6).sqrt()
    max_actual = arb(0)
    for i in range(nx1):
        x1 = -arb(1) / 2 + dx * (arb(i) + arb(1) / 2)
        for j in range(nx2):
            x2 = dt * (arb(j) + arb(1) / 2)
            for k in range(ns):
                s = ds * (arb(k) + arb(1) / 2)
                r2 = s.exp()
                y2 = r2 - x1 * x1 - x2 * x2
                if not bool(y2 > 0):
                    raise ArithmeticError("nonpositive floor-shell height")
                y = y2.sqrt()
                p = (x1, x2, y)
                gp, jac = action_jacobian(mat, *p)
                lv, lg = witness_value_gradient(ev, odd, *p)
                rv, rg = witness_value_gradient(ev, odd, *gp)
                pulled = transpose_pullback(jac, rg)
                dv = lv - rv
                dg = [lg[q] - pulled[q] for q in range(3)]

                u = s / epsilon
                hp = 30 * u * u * (1 - u) * (1 - u)
                hpp = 60 * u - 180 * u * u + 120 * u * u * u
                qs = hp / epsilon
                dq = [2 * x1 * qs / r2, 2 * x2 * qs / r2, 2 * y * qs / r2]
                lapq = -4 * y2 * hpp / (r2 * epsilon * epsilon)
                cross = sum((dq[q] * dg[q] for q in range(3)), acb(0))
                residual = lapq * dv - 2 * y2 * cross

                d0 = sqrt6 * abs(dv).upper()
                dg_e = sum((abs(q).upper() ** 2 for q in dg), arb(0)).sqrt()
                d1 = sqrt6 * y * dg_e
                gradq = y * sum((q * q for q in dq), arb(0)).sqrt()
                value = abs(lapq).upper() * d0
                gradient = 2 * abs(gradq).upper() * d1
                requested = value + gradient
                actual = sqrt6 * abs(residual).upper()

                # dx1 dx2 dy/y^3 = r2/(2 y^4) dx1 dx2 ds.
                measure = r2 / (2 * y2 * y2)
                cell = dx * dt * ds * measure
                actual2 += cell * actual * actual
                requested2 += cell * requested * requested
                value2 += cell * value * value
                gradient2 += cell * gradient * gradient
                max_actual = max_actual.union(actual).upper()
    return {
        "load_bearing": False,
        "method": "Arb midpoint values in exact (x1,x2,log(r^2)) floor-shell coordinates",
        "subdivision": [nx1, nx2, ns],
        "epsilon": str(epsilon.mid()),
        "requested_norm_estimate": upper(requested2.sqrt()),
        "commutator_norm_estimate": upper(actual2.sqrt()),
        "value_norm_estimate": upper(value2.sqrt()),
        "gradient_norm_estimate": upper(gradient2.sqrt()),
        "max_point_commutator_upper": upper(max_actual),
        "target_R_upper": upper(TARGET_R),
    }


def projected_floor_width_sweep(
    ev: DifferentialTrial,
    odd: list[acb],
    widths: list[str],
    subdivision: tuple[int, int, int],
) -> list[dict[str, Any]]:
    out = []
    for text in widths:
        print(f"projected floor diagnostic epsilon={text}", flush=True)
        q = projected_floor_midpoint_diagnostic(
            ev, odd, *subdivision, epsilon=arb(text)
        )
        q["requested_to_target_ratio"] = math.nextafter(
            q["requested_norm_estimate"] / upper(TARGET_R), math.inf
        )
        q["commutator_to_target_ratio"] = math.nextafter(
            q["commutator_norm_estimate"] / upper(TARGET_R), math.inf
        )
        out.append(q)
    return out


def certified_central_lower_bound(
    ev: DifferentialTrial,
    shells: list[int],
    mats: dict[str, list[list[acb]]],
    shell_source: dict[int, tuple[int, int, arb]],
    ny: int,
    y_start: arb = arb("1.06"),
) -> dict[str, Any]:
    """Lower-bound the actual cusp commutator where all core gates equal one."""
    if not bool((y_start * y_start).log() > arb("0.1")):
        raise ArithmeticError("central slab does not clear the floor gate")
    ye = edges(y_start, Y_PLATEAU, ny)
    bounds = certified_commutator_taylor_bounds(
        ev, shells, mats, shell_source, ny, y_start, Y_PLATEAU
    )
    return {
        "certified": bool(arb(repr(bounds["norm_lower"])) > 0),
        "domain": "x1 in [-0.39,0.39], x2 in [0.11,0.39], y in [1.06,1.20]",
        "gate_reason": (
            "all four vertical distances exceed 0.1 and log(x1^2+x2^2+y^2)>0.1; "
            "therefore the candidate core product partition has one active chart"
        ),
        "y_segments": ny,
        "positive_segments": bounds["positive_segments"],
        "commutator_norm_lower": bounds["norm_lower"],
        "commutator_norm_upper": bounds["norm_upper"],
        "target_R_upper": upper(TARGET_R),
        "exceeds_entire_target": bool(
            arb(repr(bounds["norm_lower"])).lower() > TARGET_R.upper()
        ),
        "operator_norm_upper": bounds["operator_norm_upper"],
    }


def matrix_operator_norm_upper(H: list[list[acb]]) -> arb:
    """Certified Hermitian 2-norm upper via the maximum absolute row sum."""
    out = arb(0)
    for row in H:
        value = sum((abs(q).upper() for q in row), arb(0))
        out = out.union(value).upper()
    return out


def residual_and_y_derivative(
    ev: DifferentialTrial,
    shells: list[int],
    shell_source: dict[int, tuple[int, int, arb]],
    y: arb,
) -> tuple[list[acb], list[acb]]:
    first, second, lap = cutoff_derivatives(y)
    _third, lap_prime = cutoff_third_lap_prime(y)
    r, rp = shell_radial_vectors(ev, shells, shell_source, y)
    y2 = square_nonnegative(y)
    lam = 1 + ev.r * ev.r
    rpp = []
    for i, nn in enumerate(shells):
        _cusp, _idx, mag = shell_source[nn]
        omega = 2 * ev.pi * mag
        rpp.append(
            omega * omega * r[i] + (y * rp[i] - lam * r[i]) / y2
        )
    s = [lap * r[i] - 2 * y2 * first * rp[i] for i in range(len(shells))]
    sp = [
        lap_prime * r[i]
        + (lap - 4 * y * first - 2 * y2 * second) * rp[i]
        - 2 * y2 * first * rpp[i]
        for i in range(len(shells))
    ]
    return s, sp


def residual_second_y_derivative(
    ev: DifferentialTrial,
    shells: list[int],
    shell_source: dict[int, tuple[int, int, arb]],
    y: arb,
) -> list[acb]:
    first, second, lap = cutoff_derivatives(y)
    third, lap_prime = cutoff_third_lap_prime(y)
    fourth = cutoff_fourth(y)
    r, rp = shell_radial_vectors(ev, shells, shell_source, y)
    y2 = square_nonnegative(y)
    lam = 1 + ev.r * ev.r
    rpp, rppp = [], []
    for i, nn in enumerate(shells):
        _cusp, _idx, mag = shell_source[nn]
        omega = 2 * ev.pi * mag
        q2 = omega * omega * r[i] + (y * rp[i] - lam * r[i]) / y2
        q3 = (
            omega * omega * rp[i]
            + q2 / y
            - (1 + lam) * rp[i] / y2
            + 2 * lam * r[i] / (y2 * y)
        )
        rpp.append(q2)
        rppp.append(q3)

    # s=A r+B r', where A=Delta chi and B=-2 y^2 chi'.
    A = lap
    Ap = lap_prime
    App = -3 * y * third - y2 * fourth
    B = -2 * y2 * first
    Bp = -4 * y * first - 2 * y2 * second
    Bpp = -4 * first - 8 * y * second - 2 * y2 * third
    return [
        App * r[i]
        + (2 * Ap + Bpp) * rp[i]
        + (A + 2 * Bp) * rpp[i]
        + B * rppp[i]
        for i in range(len(shells))
    ]


def certified_commutator_taylor_bounds(
    ev: DifferentialTrial,
    shells: list[int],
    mats: dict[str, list[list[acb]]],
    shell_source: dict[int, tuple[int, int, arb]],
    ny: int,
    y_start: arb,
    y_end: arb,
) -> dict[str, Any]:
    """Certified L2 bounds using point energies and shell-vector radii."""
    ye = edges(y_start, y_end, ny)
    total_lo = arb(0)
    total_hi = arb(0)
    positive = 0
    op = matrix_operator_norm_upper(mats["value"])
    op_sqrt = op.sqrt().upper()
    for k in range(ny):
        ya, yb = ye[k], ye[k + 1]
        ybox = interval(ya, yb)
        ym = (ya + yb) / 2
        radius = (yb - ya) / 2
        sm, spm = residual_and_y_derivative(ev, shells, shell_source, ym)
        em = quadratic_energy_ball(mats["value"], sm, "central midpoint residual")
        em_lo = max(em.lower(), arb(0))
        em_hi = max(em.upper(), arb(0))

        # Second-order Taylor: midpoint s and s' retain the exact Fourier
        # Gram cancellation.  Only the O(radius^2) remainder uses the coarse
        # operator norm on a shell coefficient enclosure for s''.
        sppbox = residual_second_y_derivative(ev, shells, shell_source, ybox)
        spp_coeff2 = sum(
            (abs(sppbox[i]).upper() ** 2 for i in range(len(shells))),
            arb(0),
        )
        esp_mid = quadratic_energy(mats["value"], spm, "midpoint residual derivative")
        norm_sp_mid = max(esp_mid, arb(0)).sqrt().upper()
        sup_spp = op_sqrt * spp_coeff2.sqrt().upper()
        center_lo = em_lo.sqrt().lower()
        center_hi = em_hi.sqrt().upper()
        variation = radius * norm_sp_mid + radius * radius * sup_spp / 2
        norm_lo = center_lo - variation
        norm_hi = center_hi + variation
        if bool(norm_lo > 0):
            total_lo += (yb - ya) * norm_lo * norm_lo / (yb ** 3)
            positive += 1
        total_hi += (yb - ya) * norm_hi * norm_hi / (ya ** 3)
    norm_lo = total_lo.sqrt().lower()
    norm_hi = total_hi.sqrt().upper()
    return {
        "norm_lower": lower(norm_lo),
        "norm_upper": upper(norm_hi),
        "positive_segments": positive,
        "operator_norm_upper": upper(op),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", type=Path, default=ROOT / "six_copy_hejhal_balanced_coeffs.json")
    ap.add_argument("--bits", type=int, default=160)
    ap.add_argument("--y-levels", default="8,16,32")
    ap.add_argument("--box-levels", default="")
    ap.add_argument("--diagnostic-segments", type=int, default=256)
    ap.add_argument("--lower-bound-segments", type=int, default=512)
    ap.add_argument("--sweep-plateaus", default="1.20,1.30,1.40,1.50")
    ap.add_argument("--sweep-segments", type=int, default=128)
    ap.add_argument("--floor-widths", default="0.10,0.15,0.20")
    ap.add_argument("--floor-subdivision", default="8,4,8")
    ap.add_argument("--json-out", type=Path, default=ROOT / "track_b_direct_weighted_result.json")
    ns = ap.parse_args()
    ctx.prec = max(128, ns.bits)
    data, row = parse_trial(ns.trial)
    ev = DifferentialTrial(data["parameters"]["M"], str(row["r"]), row["coefficients"])
    odd = projected_coefficients(ev)["odd"]
    keys, rows, source = difference_fourier_data(ev, odd)
    print(f"building exact half-torus Gram for {len(keys)} rational frequencies", flush=True)
    H = planar_gram(keys, rows)
    shells, mats, shell_source = shell_grams(ev, keys, H, source)
    print(f"collapsed to {len(shells)} exact radial shells", flush=True)
    levels = []
    for ny in [int(q) for q in ns.y_levels.split(",") if q.strip()]:
        print(f"direct weighted Fourier y-segments={ny}", flush=True)
        levels.append(integrate_fourier_level(
            ev, len(keys), shells, mats, shell_source, ny
        ))
    diagnostic = midpoint_diagnostic(
        ev, shells, mats, shell_source, ns.diagnostic_segments
    )

    central_H = planar_gram(
        keys,
        rows,
        (arb("-0.39"), arb("0.39")),
        (arb("0.11"), arb("0.39")),
    )
    central_shells, central_mats, central_sources = shell_grams(
        ev, keys, central_H, source
    )
    central_lower = certified_central_lower_bound(
        ev,
        central_shells,
        central_mats,
        central_sources,
        ns.lower_bound_segments,
    )
    sweep = cutoff_mass_sweep(
        ev,
        odd,
        shells,
        mats,
        shell_source,
        [q.strip() for q in ns.sweep_plateaus.split(",") if q.strip()],
        ns.sweep_segments,
        ns.sweep_segments,
    )
    projected_core = projected_core_face_diagnostic(ev, odd)
    floor_dims = tuple(int(q) for q in ns.floor_subdivision.split(","))
    if len(floor_dims) != 3 or min(floor_dims) <= 0:
        raise ValueError("floor subdivision must be nx1,nx2,ns")
    projected_floor = projected_floor_width_sweep(
        ev,
        odd,
        [q.strip() for q in ns.floor_widths.split(",") if q.strip()],
        floor_dims,
    )

    box_levels = []
    for item in [q for q in ns.box_levels.split(";") if q.strip()]:
        dims = tuple(int(q) for q in item.split(","))
        if len(dims) != 3 or min(dims) <= 0:
            raise ValueError(f"bad box subdivision {item!r}")
        print(f"independent coarse box subdivision={dims}", flush=True)
        box_levels.append(integrate_level(ev, odd, *dims))

    final = levels[-1]
    requested = arb(repr(final["requested_pointwise_weighted_norm_upper"]))
    actual = arb(repr(final["cancellation_preserving_commutator_norm_upper"]))
    requested_below = bool(requested.upper() < TARGET_R.lower())
    actual_below = bool(actual.upper() < TARGET_R.lower())
    result = {
        "status": "certified embedded-cusp contribution; total Track-B residual incomplete",
        "cusp_blend_contribution_certified": True,
        "total_direct_weighted_residual_certified": False,
        "rung4_certified": False,
        "domain": "[-1/2,1/2] x [0,1/2] x [1.01,1.20]",
        "measure": "dx1 dx2 dy / y^3",
        "cutoff": "chi_B=10t^3-15t^4+6t^5, t=(y-1.01)/0.19",
        "field_pair": "W_B=P_(Q,-)P_(-z,+)P_old F versus recorded six-copy F",
        "target_R_upper": upper(TARGET_R),
        "levels": levels,
        "midpoint_diagnostic": diagnostic,
        "certified_central_obstruction": central_lower,
        "cutoff_mass_sweep": sweep,
        "projected_core_face_diagnostic": projected_core,
        "projected_floor_direct_diagnostic": projected_floor,
        "independent_spatial_box_levels": box_levels,
        "final": final,
        "cusp_requested_bound_below_total_target": requested_below,
        "cusp_exact_commutator_bound_below_total_target": actual_below,
        "remaining_terms": [
            "core partition derivatives on complete Humbert face/edge/vertex overlaps",
            "cross terms between (1-chi_B) and the normalized core partition",
            "complete transition and elliptic-chart incidence inventory",
        ],
        "theorem_note": (
            "The requested pointwise bound is certified for the cusp blend only. "
            "The cancellation-preserving commutator is a sharper application of the "
            "same Laplacian product rule and does not change Theorem D(K)."
        ),
        "trial": str(ns.trial.resolve()),
        "trial_sha256": hashlib.sha256(ns.trial.read_bytes()).hexdigest(),
        "precision_bits": int(ctx.prec),
        "r": str(row["r"]),
        "bessel_evaluations": {
            "direct_or_mean_value_K": ev.direct_k_count + ev.mean_value_k_count,
            "shifted_order_cache_entries": len(ev._shifted_k_cache),
        },
    }
    ns.json_out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "cusp_blend_contribution_certified": True,
        "requested_upper": final["requested_pointwise_weighted_norm_upper"],
        "commutator_upper": final["cancellation_preserving_commutator_norm_upper"],
        "target": upper(TARGET_R),
        "requested_below_target": requested_below,
        "commutator_below_target": actual_below,
    }, indent=2))
    print(ns.json_out.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
