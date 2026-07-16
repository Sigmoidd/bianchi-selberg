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
from concurrent.futures import ProcessPoolExecutor
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
from track_b_floor_taylor import (
    BesselAudit,
    TaylorJet,
    TaylorModel,
    begin_arithmetic_audit,
    end_arithmetic_audit,
    polynomial_l2_squared_upper,
    radial_taylor_jet,
)


ROOT = Path(__file__).resolve().parent
Y_SUPPORT = arb("1.01")
Y_PLATEAU = arb("1.20")
WIDTH = Y_PLATEAU - Y_SUPPORT
TARGET_R = arb("0.0101903405004245")
FLOOR_WEIGHT_FORMULA_ID = "track-b-floor-quintic-logrho-width-0.30/v1"


def canonical_hash(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def interval(a: arb, b: arb) -> arb:
    return a.union(b)


def square_nonnegative(x: arb) -> arb:
    ax = abs(x)
    return (ax.lower() * ax.lower()).union(ax.upper() * ax.upper())


def cube_monotone(x: arb) -> arb:
    """Range of x^3 using monotonicity, without repeated-ball dependency."""
    return (x.lower() ** 3).union(x.upper() ** 3)


class Jet2:
    """Three-variable complex interval jet through second order."""

    def __init__(self, v: acb, g: list[acb] | None = None,
                 h: list[list[acb]] | None = None):
        self.v = acb(v)
        self.g = g if g is not None else [acb(0) for _ in range(3)]
        self.h = h if h is not None else [
            [acb(0) for _ in range(3)] for _ in range(3)
        ]

    @staticmethod
    def coerce(q: Any) -> "Jet2":
        return q if isinstance(q, Jet2) else Jet2(acb(q))

    @staticmethod
    def variable(q: arb, index: int) -> "Jet2":
        g = [acb(0) for _ in range(3)]
        g[index] = acb(1)
        return Jet2(acb(q), g)

    def __add__(self, other: Any) -> "Jet2":
        b = self.coerce(other)
        return Jet2(
            self.v + b.v,
            [self.g[i] + b.g[i] for i in range(3)],
            [[self.h[i][j] + b.h[i][j] for j in range(3)] for i in range(3)],
        )

    __radd__ = __add__

    def __neg__(self) -> "Jet2":
        return Jet2(
            -self.v,
            [-q for q in self.g],
            [[-q for q in row] for row in self.h],
        )

    def __sub__(self, other: Any) -> "Jet2":
        return self + (-self.coerce(other))

    def __rsub__(self, other: Any) -> "Jet2":
        return self.coerce(other) - self

    def __mul__(self, other: Any) -> "Jet2":
        b = self.coerce(other)
        g = [self.g[i] * b.v + self.v * b.g[i] for i in range(3)]
        h = [[acb(0) for _ in range(3)] for _ in range(3)]
        for i in range(3):
            for j in range(3):
                h[i][j] = (
                    self.h[i][j] * b.v
                    + self.g[i] * b.g[j]
                    + self.g[j] * b.g[i]
                    + self.v * b.h[i][j]
                )
        return Jet2(self.v * b.v, g, h)

    __rmul__ = __mul__

    def unary(self, value: acb, first: acb, second: acb) -> "Jet2":
        g = [first * self.g[i] for i in range(3)]
        h = [[
            second * self.g[i] * self.g[j] + first * self.h[i][j]
            for j in range(3)
        ] for i in range(3)]
        return Jet2(value, g, h)

    def reciprocal(self) -> "Jet2":
        if self.v.imag.is_zero():
            real = self.v.real
            if not real.contains(0):
                return self.unary(
                    acb(1 / real),
                    acb(-1 / square_nonnegative(real)),
                    acb(2 / cube_monotone(real)),
                )
        return self.unary(1 / self.v, -1 / (self.v * self.v), 2 / (self.v ** 3))

    def __truediv__(self, other: Any) -> "Jet2":
        return self * self.coerce(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "Jet2":
        return self.coerce(other) / self

    def __pow__(self, n: int) -> "Jet2":
        if n == 0:
            return Jet2(acb(1))
        if n < 0:
            return (self ** (-n)).reciprocal()
        out = Jet2(acb(1))
        base = self
        k = n
        while k:
            if k & 1:
                out = out * base
            base = base * base
            k >>= 1
        return out

    def exp(self) -> "Jet2":
        if self.v.imag.is_zero():
            value = self.v.real.exp()
            return self.unary(acb(value), acb(value), acb(value))
        value = self.v.exp()
        return self.unary(value, value, value)

    def log(self) -> "Jet2":
        if self.v.imag.is_zero() and bool(self.v.real > 0):
            real = self.v.real
            return self.unary(
                acb(real.log()), acb(1 / real), acb(-1 / square_nonnegative(real))
            )
        return self.unary(self.v.log(), 1 / self.v, -1 / (self.v * self.v))

    def sqrt(self) -> "Jet2":
        if self.v.imag.is_zero() and bool(self.v.real > 0):
            value = self.v.real.sqrt()
            return self.unary(
                acb(value),
                acb(1 / (2 * value)),
                acb(-1 / (4 * cube_monotone(value))),
            )
        value = self.v.sqrt()
        return self.unary(value, 1 / (2 * value), -1 / (4 * value ** 3))

    def sin(self) -> "Jet2":
        if self.v.imag.is_zero():
            value = self.v.real.sin()
            return self.unary(acb(value), acb(self.v.real.cos()), acb(-value))
        value = self.v.sin()
        return self.unary(value, self.v.cos(), -value)

    def cos(self) -> "Jet2":
        if self.v.imag.is_zero():
            value = self.v.real.cos()
            return self.unary(acb(value), acb(-self.v.real.sin()), acb(-value))
        value = self.v.cos()
        return self.unary(value, -self.v.sin(), -value)


def compose_radial(y: Jet2, value: acb, first: acb, second: acb) -> Jet2:
    return y.unary(value, first, second)


def radial_jet(
    ev: DifferentialTrial,
    mode_index: int,
    magnitude: arb,
    y: Jet2,
) -> Jet2:
    """Second-order jet of y K_ir(2*pi*magnitude*y)."""
    if not y.v.imag.contains(0) or not bool(y.v.real > 0):
        raise ArithmeticError(f"jet height not positive real: {y.v}")
    yr = y.v.real
    omega = 2 * ev.pi * magnitude
    arg = omega * yr
    kval = ev._kir(0, mode_index, magnitude, yr)
    kprime = -(ev._kir_order(ev.order - 1, arg)
               + ev._kir_order(ev.order + 1, arg)) / 2
    radial = y.v * kval
    radial_y = kval + y.v * omega * kprime
    lam = 1 + ev.r * ev.r
    radial_yy = (
        omega * omega * radial
        + (y.v * radial_y - lam * radial) / acb(square_nonnegative(yr))
    )
    return compose_radial(y, radial, radial_y, radial_yy)


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


def witness_jet(
    ev: DifferentialTrial,
    odd: list[acb],
    x1: Jet2,
    x2: Jet2,
    y: Jet2,
) -> Jet2:
    """Second-order jet of the projected scalar Whittaker field."""
    if not y.v.imag.contains(0) or not bool(y.v.real > 0):
        raise ArithmeticError(f"jet height not positive real: {y.v}")
    total = Jet2(acb(0))
    for k, mode in enumerate(ev.modes_inf):
        u, v, mag = ev._frequency(0, mode)
        rjet = radial_jet(ev, k, mag, y)
        phase_arg = acb(0, 2 * ev.pi) * (u * x1 + v * x2)
        total += odd[k] * phase_arg.exp() * rjet
    return total


def reflection_orbits(ev: DifferentialTrial) -> list[tuple[int, int]]:
    """Infinity-mode orbits under the floor reflection (u,v)->(-u,v)."""
    cached = getattr(ev, "_floor_reflection_orbits", None)
    if cached is not None:
        return cached
    index = {(a, b): k for k, (a, b, _nn) in enumerate(ev.modes_inf)}
    visited: set[int] = set()
    out: list[tuple[int, int]] = []
    for k, (a, b, _nn) in enumerate(ev.modes_inf):
        if k in visited:
            continue
        kr = index.get((-a, b))
        if kr is None:
            raise AssertionError(f"missing floor-reflected mode {(-a, b)}")
        visited.add(k)
        visited.add(kr)
        out.append((k, kr))
    if len(visited) != len(ev.modes_inf):
        raise AssertionError("floor-reflection orbit ledger is incomplete")
    setattr(ev, "_floor_reflection_orbits", out)
    return out


def quarter_turn_orbits(
    ev: DifferentialTrial, odd: list[acb]
) -> tuple[list[tuple[int, int, int]], dict[tuple[int, int], tuple[int, int]]]:
    r"""Exact Q-orbits and the sign of A_key relative to each representative.

    A_(a,b)(x,t)=2(cos(2*pi*(a*x+b*t))-cos(2*pi*(-b*x+a*t)))
    satisfies A_(Q(a,b))=-A_(a,b).  The projected coefficients obey the
    matching sign rule exactly.
    """
    cached = getattr(ev, "_floor_quarter_turn_orbits", None)
    if cached is not None:
        return cached
    index = {(a, b): k for k, (a, b, _nn) in enumerate(ev.modes_inf)}
    visited: set[tuple[int, int]] = set()
    records: list[tuple[int, int, int]] = []
    key_map: dict[tuple[int, int], tuple[int, int]] = {}
    for key in sorted(index):
        if key in visited:
            continue
        a, b = key
        orbit = [(a, b), (-b, a), (-a, -b), (b, -a)]
        rep = min(orbit)
        a, b = rep
        ordered = [(a, b), (-b, a), (-a, -b), (b, -a)]
        orbit_id = len(records)
        krep = index[rep]
        for power, qkey in enumerate(ordered):
            if qkey not in index:
                raise AssertionError(f"missing quarter-turn mode {qkey}")
            sign = 1 if power % 2 == 0 else -1
            defect = odd[index[qkey]] - sign * odd[krep]
            # The identity is algebraic in projected_coefficients(); Arb's
            # separately rounded evaluations need only overlap at zero.
            if not (defect.real.contains(0) and defect.imag.contains(0)):
                raise AssertionError(
                    f"projected coefficient fails exact Q sign at {qkey}: {defect}"
                )
            visited.add(qkey)
            key_map[qkey] = (orbit_id, sign)
        records.append((a, b, krep))
    if len(visited) != len(ev.modes_inf):
        raise AssertionError("quarter-turn orbit ledger is incomplete")
    out = (records, key_map)
    setattr(ev, "_floor_quarter_turn_orbits", out)
    return out


def quarter_turn_angular_jet(a: int, b: int, x: Jet2, t: Jet2) -> Jet2:
    scale = 2 * arb.pi()
    first = scale * (a * x + b * t)
    second = scale * (-b * x + a * t)
    return 2 * (first.cos() - second.cos())


def witness_S_defect_param_jet(
    ev: DifferentialTrial,
    odd: list[acb],
    x: Jet2,
    t: Jet2,
    s: Jet2,
) -> Jet2:
    r"""Jet of W(x,t,y)-W(S(x,t,y)) in independent (x,t,s) coordinates.

    Here rho=exp(s), y^2=rho-x^2-t^2, and
    S(x,t,y)=(-x/rho,t/rho,y/rho).  Nontrivial frequency orbits are
    evaluated in cosine/sine form.  Thus the reflection-even cancellation
    is performed before interval absolute values are taken.
    """
    rho = s.exp()
    y2 = rho - x * x - t * t
    if not y2.v.imag.contains(0) or not bool(y2.v.real > 0):
        raise ArithmeticError(f"floor parameter jet height squared not positive: {y2.v}")
    y = y2.sqrt()
    ys = y / rho
    total = Jet2(acb(0))
    shell_angular: dict[int, tuple[Jet2, Jet2, int, arb]] = {}

    def add_shell(
        nn: int, k: int, mag: arb, left_part: Jet2, right_part: Jet2
    ) -> None:
        old = shell_angular.get(nn)
        if old is None:
            shell_angular[nn] = (left_part, right_part, k, mag)
            return
        left_old, right_old, k_old, mag_old = old
        if ev.modes_inf[k_old][2] != nn:
            raise AssertionError("incorrect radial shell representative")
        shell_angular[nn] = (
            left_old + left_part,
            right_old + right_part,
            k_old,
            mag_old,
        )

    records, key_map = quarter_turn_orbits(ev, odd)
    visited: set[int] = set()
    xp, tp = x / rho, t / rho
    for orbit_id, (a, b, k) in enumerate(records):
        if orbit_id in visited:
            continue
        reflected_key = (-a, b)
        reflected_id, reflected_sign = key_map[reflected_key]
        ar, br, kr = records[reflected_id]
        _u, _v, mag = ev._frequency(0, ev.modes_inf[k])
        if ev.modes_inf[k][2] != ev.modes_inf[kr][2]:
            raise AssertionError("reflection-paired Q-orbits have different radii")
        basis_left = quarter_turn_angular_jet(a, b, x, t)
        basis_scaled = quarter_turn_angular_jet(a, b, xp, tp)
        nn = ev.modes_inf[k][2]

        if reflected_id == orbit_id:
            add_shell(
                nn,
                k,
                mag,
                odd[k] * basis_left,
                -reflected_sign * odd[k] * basis_scaled,
            )
            visited.add(orbit_id)
            continue

        reflected_left = quarter_turn_angular_jet(-a, b, x, t)
        reflected_scaled = quarter_turn_angular_jet(-a, b, xp, tp)
        aligned_reflected_coefficient = reflected_sign * odd[kr]
        average = (odd[k] + aligned_reflected_coefficient) / 2
        difference = (odd[k] - aligned_reflected_coefficient) / 2
        add_shell(
            nn,
            k,
            mag,
            average * (basis_left + reflected_left)
            + difference * (basis_left - reflected_left),
            -average * (basis_scaled + reflected_scaled)
            + difference * (basis_scaled - reflected_scaled),
        )
        visited.add(orbit_id)
        visited.add(reflected_id)
    if len(visited) != len(records):
        raise AssertionError("dihedral floor-orbit ledger is incomplete")
    for nn in sorted(shell_angular):
        left_part, right_part, k, mag = shell_angular[nn]
        total += (
            left_part * radial_jet(ev, k, mag, y)
            + right_part * radial_jet(ev, k, mag, ys)
        )
    return total


def floor_residual_param_jet(
    ev: DifferentialTrial,
    odd: list[acb],
    x_ball: arb,
    t_ball: arb,
    s_ball: arb,
    epsilon: arb,
) -> tuple[acb, list[acb]]:
    r"""Floor commutator and its (x,t,s) gradient with shell correlation exact.

    For D=W-W o S and q=q(s), harmonicity of s=log(x^2+t^2+y^2)
    for the hyperbolic Laplacian gives

        [Delta,q]D = -4 y^2/rho * (q_ss D + q_s(x D_x+t D_t+2 D_s)).

    Only the second jet of D is needed to enclose the first derivatives of
    this already-subtracted residual.
    """
    x = Jet2.variable(x_ball, 0)
    t = Jet2.variable(t_ball, 1)
    s = Jet2.variable(s_ball, 2)
    rho = s.exp()
    y2 = rho - x * x - t * t
    if not y2.v.imag.contains(0) or not bool(y2.v.real > 0):
        raise ArithmeticError(f"floor parameter height squared not positive: {y2.v}")
    defect = witness_S_defect_param_jet(ev, odd, x, t, s)

    u = s / epsilon
    qs = 30 * u * u * (1 - u) * (1 - u) / epsilon
    qss = (60 * u - 180 * u * u + 120 * u * u * u) / (epsilon * epsilon)
    prefactor = -4 * y2 / rho

    directional = x.v * defect.g[0] + t.v * defect.g[1] + 2 * defect.g[2]
    bracket = qss.v * defect.v + qs.v * directional
    directional_grad = []
    for ell in range(3):
        directional_grad.append(
            (defect.g[0] if ell == 0 else acb(0))
            + x.v * defect.h[0][ell]
            + (defect.g[1] if ell == 1 else acb(0))
            + t.v * defect.h[1][ell]
            + 2 * defect.h[2][ell]
        )
    bracket_grad = [
        qss.g[ell] * defect.v
        + qss.v * defect.g[ell]
        + qs.g[ell] * directional
        + qs.v * directional_grad[ell]
        for ell in range(3)
    ]
    value = prefactor.v * bracket
    gradient = [
        prefactor.g[ell] * bracket + prefactor.v * bracket_grad[ell]
        for ell in range(3)
    ]
    return value, gradient


def quarter_turn_angular_taylor(
    a: int, b: int, x: TaylorJet, t: TaylorJet
) -> TaylorJet:
    scale = 2 * arb.pi()
    first = scale * (a * x + b * t)
    second = scale * (-b * x + a * t)
    return 2 * (first.cos() - second.cos())


def witness_S_defect_taylor(
    ev: DifferentialTrial,
    odd: list[acb],
    x: TaylorJet,
    t: TaylorJet,
    s: TaylorJet,
    audit: BesselAudit,
    rho: TaylorJet | None = None,
    y2: TaylorJet | None = None,
) -> TaylorJet:
    r"""Validated model of D=W_B o S-W_B after all exact orbit reductions."""
    degree = x.value.degree
    rho = s.exp() if rho is None else rho
    y2 = rho - x * x - t * t if y2 is None else y2
    y = y2.sqrt()
    ys = y / rho
    xp, tp = x / rho, t / rho
    records, key_map = quarter_turn_orbits(ev, odd)
    visited: set[int] = set()
    shell_angular: dict[int, tuple[TaylorJet, TaylorJet, int, arb]] = {}

    def add_shell(
        nn: int,
        k: int,
        magnitude: arb,
        left_part: TaylorJet,
        right_part: TaylorJet,
    ) -> None:
        old = shell_angular.get(nn)
        if old is None:
            shell_angular[nn] = (left_part, right_part, k, magnitude)
            return
        left_old, right_old, k_old, magnitude_old = old
        shell_angular[nn] = (
            left_old + left_part,
            right_old + right_part,
            k_old,
            magnitude_old,
        )

    for orbit_id, (a, b, k) in enumerate(records):
        if orbit_id in visited:
            continue
        reflected_id, reflected_sign = key_map[(-a, b)]
        _ar, _br, kr = records[reflected_id]
        if ev.modes_inf[k][2] != ev.modes_inf[kr][2]:
            raise AssertionError("reflection-paired Taylor orbits have different shells")
        _u, _v, magnitude = ev._frequency(0, ev.modes_inf[k])
        nn = ev.modes_inf[k][2]
        basis_left = quarter_turn_angular_taylor(a, b, x, t)
        basis_scaled = quarter_turn_angular_taylor(a, b, xp, tp)
        if reflected_id == orbit_id:
            add_shell(
                nn,
                k,
                magnitude,
                odd[k] * basis_left,
                -reflected_sign * odd[k] * basis_scaled,
            )
            visited.add(orbit_id)
            continue

        reflected_left = quarter_turn_angular_taylor(-a, b, x, t)
        reflected_scaled = quarter_turn_angular_taylor(-a, b, xp, tp)
        aligned = reflected_sign * odd[kr]
        average = (odd[k] + aligned) / 2
        difference = (odd[k] - aligned) / 2
        add_shell(
            nn,
            k,
            magnitude,
            average * (basis_left + reflected_left)
            + difference * (basis_left - reflected_left),
            -average * (basis_scaled + reflected_scaled)
            + difference * (basis_scaled - reflected_scaled),
        )
        visited.add(orbit_id)
        visited.add(reflected_id)
    if len(visited) != len(records):
        raise AssertionError("Taylor dihedral orbit ledger is incomplete")

    left_minus_right = TaylorJet.constant(degree, 0)
    for nn in sorted(shell_angular):
        left_part, right_part, k, magnitude = shell_angular[nn]
        left_minus_right += (
            left_part * radial_taylor_jet(ev, k, magnitude, y, audit)
            + right_part * radial_taylor_jet(ev, k, magnitude, ys, audit)
        )
    # The certificate convention requested for the floor is D=W o S-W.
    return -left_minus_right


def floor_cell_taylor_model(
    ev: DifferentialTrial,
    odd: list[acb],
    x_bounds: tuple[arb, arb],
    t_bounds: tuple[arb, arb],
    s_bounds: tuple[arb, arb],
    width: arb,
    degree: int | tuple[int, int, int],
    collect_arithmetic_audit: bool = False,
) -> dict[str, Any]:
    """Construct G=sqrt(6J) R_S=P+E on one normalized floor cell."""
    if isinstance(degree, int):
        if degree not in (4, 6, 8, 10, 12):
            raise ValueError("floor Taylor degree must be one of 4,6,8,10,12")
    elif len(degree) != 3 or min(degree) < 1 or max(degree) > 12:
        raise ValueError("tensor degrees must be three integers in [1,12]")
    arithmetic = begin_arithmetic_audit() if collect_arithmetic_audit else None
    try:
        xa, xb = x_bounds
        ta, tb = t_bounds
        sa, sb = s_bounds
        cx, ct, cs = (xa + xb) / 2, (ta + tb) / 2, (sa + sb) / 2
        hx, ht, hs = (xb - xa) / 2, (tb - ta) / 2, (sb - sa) / 2
        x = TaylorJet.variable(degree, cx, hx, 0)
        t = TaylorJet.variable(degree, ct, ht, 1)
        s = TaylorJet.variable(degree, cs, hs, 2)
        audit = BesselAudit()
        rho = s.exp()
        y2 = rho - x * x - t * t
        defect = witness_S_defect_taylor(ev, odd, x, t, s, audit, rho=rho, y2=y2)

        u = s.value / width
        qs = 30 * u * u * (1 - u) * (1 - u) / width
        qss = (60 * u - 180 * u * u + 120 * u * u * u) / (width * width)
        directional = (
            x.value * defect.gradient[0]
            + t.value * defect.gradient[1]
            + 2 * defect.gradient[2]
        )
        residual = (-4 * y2.value / rho.value) * (
            qss * defect.value + qs * directional
        )
        sqrt_jacobian = (s.value / 2).exp() / (arb(2).sqrt() * y2.value)
        weighted = arb(6).sqrt() * sqrt_jacobian * residual
    finally:
        arithmetic = end_arithmetic_audit() if collect_arithmetic_audit else None
    polynomial_l2_sq, legendre_count = polynomial_l2_squared_upper(
        weighted, hx, ht, hs
    )
    polynomial_l2 = polynomial_l2_sq.sqrt().upper()
    volume = 8 * hx * ht * hs
    remainder_l2 = (weighted.error * volume.sqrt()).upper()
    combined = (polynomial_l2 + remainder_l2).upper()
    return {
        "model": weighted,
        "residual_model": residual,
        "polynomial_l2_squared_upper_ball": polynomial_l2_sq,
        "polynomial_l2_upper_ball": polynomial_l2,
        "remainder_supremum_upper_ball": weighted.error.upper(),
        "remainder_l2_upper_ball": remainder_l2,
        "combined_cell_l2_upper_ball": combined,
        "polynomial_coefficient_count": len(weighted.coeff),
        "legendre_coefficient_count": legendre_count,
        "bessel_direct_count": audit.direct_count,
        "bessel_fallback_count": audit.fallback_count,
        "bessel_majorant_count": audit.majorant_count,
        "bessel_real_order_counts": dict(sorted(audit.real_order_counts.items())),
        "arithmetic_audit": None if arithmetic is None else arithmetic.summary(),
        "fourier_tail_upper_ball": arb(0),
        "center": (cx, ct, cs),
        "halfwidth": (hx, ht, hs),
    }


def projected_symmetry_algebra_certificate(
    ev: DifferentialTrial, odd: list[acb]
) -> dict[str, Any]:
    """Prove the zero core relations by exact frequency permutations."""
    index = {(a, b): k for k, (a, b, _nn) in enumerate(ev.modes_inf)}
    minus_ok = True
    quarter_ok = True
    for k, (a, b, _nn) in enumerate(ev.modes_inf):
        dm = odd[k] - odd[index[(-a, -b)]]
        dq = odd[k] + odd[index[(-b, a)]]
        minus_ok = minus_ok and dm.real.contains(0) and dm.imag.contains(0)
        quarter_ok = quarter_ok and dq.real.contains(0) and dq.imag.contains(0)
    periodic = all(isinstance(a, int) and isinstance(b, int)
                   for a, b, _nn in ev.modes_inf)
    certified = bool(periodic and minus_ok and quarter_ok)
    exact_zero = {
        name: {
            "value_defect": "0",
            "first_gradient_defect": "0",
            "proof": proof,
        }
        for name, proof in {
            "T1": "integer-frequency phase permutation; differentiate the identity",
            "R": "(-z)-even coefficient permutation (a,b)->(-a,-b)",
            "TiR": "composition of integer i-translation with (-z)-evenness",
        }.items()
    }
    return {
        "certified": certified,
        "method": "exact lattice permutations induced by the explicit projectors",
        "minus_z_even_coefficient_identities": minus_ok,
        "quarter_turn_odd_coefficient_identities": quarter_ok,
        "integer_periodicity": periodic,
        "relations": exact_zero,
        "laplace_eigen_equation": {
            "certified": True,
            "identity": "(Delta-(1+r^2)) y K_(ir)(2*pi*|m|y) exp(2*pi*i<m,x>)=0",
            "proof": "modified-Bessel differential equation, modewise; linear projectors preserve it",
        },
    }


def _serialized_floor_cell(
    index: tuple[int, int, int],
    x_bounds: tuple[arb, arb],
    t_bounds: tuple[arb, arb],
    s_bounds: tuple[arb, arb],
    degree: int | tuple[int, int, int],
    bits: int,
    cell: dict[str, Any],
) -> dict[str, Any]:
    return {
        "cell_index": list(index),
        # Keep endpoints explicitly.  python-flint ``arb.union(a, b)`` is an
        # error-ball union about zero, not a serialization of the convex
        # interval [a,b]; the legacy *_interval fields are retained only for
        # schema compatibility and must not be used to reconstruct cells.
        "x1_bounds": [str(x_bounds[0]), str(x_bounds[1])],
        "x2_bounds": [str(t_bounds[0]), str(t_bounds[1])],
        "s_bounds": [str(s_bounds[0]), str(s_bounds[1])],
        "x1_interval": str(interval(*x_bounds)),
        "x2_interval": str(interval(*t_bounds)),
        "s_interval": str(interval(*s_bounds)),
        "polynomial_degree": list(degree) if isinstance(degree, tuple) else degree,
        "polynomial_coefficient_count": cell["polynomial_coefficient_count"],
        "legendre_coefficient_count": cell["legendre_coefficient_count"],
        "polynomial_l2_squared_upper": str(cell["polynomial_l2_squared_upper_ball"]),
        "polynomial_l2_upper": str(cell["polynomial_l2_upper_ball"]),
        "remainder_supremum_upper": str(cell["remainder_supremum_upper_ball"]),
        "remainder_l2_upper": str(cell["remainder_l2_upper_ball"]),
        "combined_cell_l2_upper": str(cell["combined_cell_l2_upper_ball"]),
        "fourier_tail_upper": str(cell["fourier_tail_upper_ball"]),
        "bessel_direct_count": cell["bessel_direct_count"],
        "bessel_fallback_count": cell["bessel_fallback_count"],
        "bessel_majorant_count": cell["bessel_majorant_count"],
        "bessel_real_order_counts": cell["bessel_real_order_counts"],
        "arithmetic_audit": cell["arithmetic_audit"],
        "working_precision": bits,
    }


_FLOOR_WORKER_STATE: dict[str, Any] = {}


def _initialize_floor_worker(trial_path: str, bits: int) -> None:
    ctx.prec = bits
    data, row = parse_trial(Path(trial_path))
    ev = DifferentialTrial(data["parameters"]["M"], str(row["r"]), row["coefficients"])
    _FLOOR_WORKER_STATE["ev"] = ev
    _FLOOR_WORKER_STATE["odd"] = projected_coefficients(ev)["odd"]
    _FLOOR_WORKER_STATE["bits"] = bits


def _floor_cell_worker(payload: tuple[Any, ...]) -> dict[str, Any]:
    if len(payload) == 6:
        index, xb, tb, sb, width_text, degree = payload
        collect_arithmetic_audit = False
    else:
        index, xb, tb, sb, width_text, degree, collect_arithmetic_audit = payload
    x_bounds = (arb(xb[0]), arb(xb[1]))
    t_bounds = (arb(tb[0]), arb(tb[1]))
    s_bounds = (arb(sb[0]), arb(sb[1]))
    cell = floor_cell_taylor_model(
        _FLOOR_WORKER_STATE["ev"],
        _FLOOR_WORKER_STATE["odd"],
        x_bounds,
        t_bounds,
        s_bounds,
        arb(width_text),
        degree,
        bool(collect_arithmetic_audit),
    )
    return _serialized_floor_cell(
        index, x_bounds, t_bounds, s_bounds, degree,
        _FLOOR_WORKER_STATE["bits"], cell,
    )


def certified_allowed_floor_budget(
    mass_certificate: dict[str, Any], width_tolerance: arb
) -> arb:
    """R < mu_B*width_tolerance/2 implies full spectral width < tolerance."""
    if not mass_certificate.get("theorem_DK_projected_mass_admissible", False):
        raise ValueError("projected-mass certificate is not theorem-admissible")
    mu_lower = arb(str(
        mass_certificate["plateau_construction"]["certified_mu_B_lower"]
    ))
    if not bool(mu_lower > 0) or not bool(width_tolerance > 0):
        raise ValueError("allowed floor budget inputs are not positive")
    return (mu_lower * width_tolerance / 2)


def certify_projected_floor_taylor_grid(
    ev: DifferentialTrial,
    odd: list[acb],
    trial_path: Path,
    width: arb,
    subdivision: tuple[int, int, int],
    degree: int | tuple[int, int, int],
    bits: int,
    allowed_budget: arb,
    mass_certificate: dict[str, Any],
    partition_certificate: dict[str, Any],
    audit_path: Path | None = None,
    workers: int = 1,
    collect_arithmetic_audit: bool = False,
) -> dict[str, Any]:
    nx, nt, ns = subdivision
    xe = edges(-arb(1) / 2, arb(1) / 2, nx)
    te = edges(arb(0), arb(1) / 2, nt)
    se = edges(arb(0), width, ns)
    payloads = []
    for i in range(nx):
        for j in range(nt):
            for k in range(ns):
                payloads.append((
                    (i, j, k),
                    (str(xe[i]), str(xe[i + 1])),
                    (str(te[j]), str(te[j + 1])),
                    (str(se[k]), str(se[k + 1])),
                    str(width),
                    degree,
                    collect_arithmetic_audit,
                ))

    if workers > 1:
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=_initialize_floor_worker,
            initargs=(str(trial_path.resolve()), bits),
        ) as pool:
            records = list(pool.map(_floor_cell_worker, payloads))
    else:
        records = []
        for index, xb, tb, sb, width_text, model_degree, arithmetic_audit in payloads:
            x_bounds = (arb(xb[0]), arb(xb[1]))
            t_bounds = (arb(tb[0]), arb(tb[1]))
            s_bounds = (arb(sb[0]), arb(sb[1]))
            cell = floor_cell_taylor_model(
                ev, odd, x_bounds, t_bounds, s_bounds,
                arb(width_text), model_degree, arithmetic_audit,
            )
            records.append(_serialized_floor_cell(
                index, x_bounds, t_bounds, s_bounds, model_degree, bits, cell
            ))

    polynomial_sq = arb(0)
    remainder_sq = arb(0)
    combined_sq = arb(0)
    direct_count = 0
    fallback_count = 0
    majorant_count = 0
    worst_record = None
    worst_remainder = arb(0)
    for record in records:  # deterministic lexicographic outward summation
        polynomial_sq += arb(record["polynomial_l2_squared_upper"]).upper()
        rem = arb(record["remainder_l2_upper"]).upper()
        combined = arb(record["combined_cell_l2_upper"]).upper()
        remainder_sq += rem * rem
        combined_sq += combined * combined
        direct_count += int(record["bessel_direct_count"])
        fallback_count += int(record["bessel_fallback_count"])
        majorant_count += int(record["bessel_majorant_count"])
        rem_sup = arb(record["remainder_supremum_upper"]).upper()
        if worst_record is None or bool(rem_sup > worst_remainder):
            worst_record = record
            worst_remainder = rem_sup

    if audit_path is not None:
        audit_path.write_text(
            "".join(json.dumps(record, sort_keys=True) + "\n" for record in records),
            encoding="utf-8",
        )

    polynomial_l2 = polynomial_sq.sqrt().upper()
    remainder_l2 = remainder_sq.sqrt().upper()
    floor_l2 = combined_sq.sqrt().upper()
    budget_lower = allowed_budget.lower()
    margin = (budget_lower - floor_l2.upper()).lower()
    ratio = (floor_l2.upper() / budget_lower).upper()
    geometry = partition_certificate.get("floor_collar_incidence", {})
    geometry_ok = bool(
        geometry.get("certified", False)
        and geometry.get("floor_collar", {}).get("containment_verified", False)
        and geometry.get("open_collar_incident_codimension_one_words") == ["S"]
    )
    symmetries = projected_symmetry_algebra_certificate(ev, odd)
    remainder_ok = bool(
        floor_l2.is_finite() and remainder_l2.is_finite()
        and all(arb(r["remainder_supremum_upper"]).is_finite() for r in records)
    )
    strict = bool(floor_l2.upper() < budget_lower)
    local_certified = bool(
        strict and remainder_ok and fallback_count == 0
        and geometry_ok and symmetries["certified"]
    )
    if not bool(margin > 0):
        status = "RED"
    elif bool(ratio <= arb("0.90")):
        status = "GREEN"
    else:
        status = "YELLOW"
    global_partition = bool(partition_certificate.get("certified", False))
    global_weights = bool(partition_certificate.get("partition_constants_certified", False))
    witness_mass = bool(mass_certificate.get("theorem_DK_projected_mass_admissible", False))
    return {
        "label": "LOCAL FLOOR CERTIFICATE" if local_certified else "DIAGNOSTIC",
        "floor_width": str(width),
        "floor_weight_formula_id": FLOOR_WEIGHT_FORMULA_ID,
        "floor_geometry_incidence_hash": canonical_hash(geometry),
        "grid": list(subdivision),
        "polynomial_degree": degree,
        "arb_bits": bits,
        "polynomial_l2_upper": str(polynomial_l2),
        "remainder_l2_upper": str(remainder_l2),
        "floor_l2_upper": str(floor_l2),
        "allowed_budget_interval": str(allowed_budget),
        "allowed_budget_lower": str(budget_lower),
        "certified_margin_lower": str(margin),
        "certified_ratio_upper": str(ratio),
        "worst_cell": worst_record["cell_index"] if worst_record else None,
        "worst_cell_remainder": str(worst_remainder),
        "maximum_fourier_bessel_tail": "0",
        "bessel_direct_count": direct_count,
        "bessel_fallback_count": fallback_count,
        "bessel_majorant_count": majorant_count,
        "fourier_tail_upper": "0",
        "continuum_remainder_certified": remainder_ok,
        "floor_residual_certified": local_certified,
        "stability_check_passed": False,
        "geometry_incidence_certified": geometry_ok,
        "projected_symmetries_certified": symmetries["certified"],
        "projected_symmetry_certificate": symmetries,
        "global_partition_certified": global_partition,
        "global_weight_bounds_certified": global_weights,
        "witness_mass_certified": witness_mass,
        "rung4_integrator_comparison_certified": False,
        "rung4_certified": False,
        "status": status,
        "cell_count": len(records),
        "audit_jsonl": None if audit_path is None else str(audit_path.resolve()),
    }


def floor_residual_jet(
    ev: DifferentialTrial,
    odd: list[acb],
    x1_ball: arb,
    x2_ball: arb,
    y_ball: arb,
    epsilon: arb,
) -> tuple[acb, list[acb]]:
    """Value and original-coordinate gradient of the scalar floor commutator."""
    x = Jet2.variable(x1_ball, 0)
    t = Jet2.variable(x2_ball, 1)
    y = Jet2.variable(y_ball, 2)
    r2 = x * x + t * t + y * y
    xs = -x / r2
    ts = t / r2
    ys = y / r2
    defect = witness_jet(ev, odd, x, t, y) - witness_jet(ev, odd, xs, ts, ys)

    u = r2.log() / epsilon
    hp = 30 * u * u * (1 - u) * (1 - u)
    hpp = 60 * u - 180 * u * u + 120 * u * u * u
    qs = hp / epsilon
    dq = [2 * x * qs / r2, 2 * t * qs / r2, 2 * y * qs / r2]
    lapq = -4 * y * y * hpp / (r2 * epsilon * epsilon)

    cross = sum((dq[i].v * defect.g[i] for i in range(3)), acb(0))
    value = lapq.v * defect.v - 2 * y.v * y.v * cross
    grad = []
    y2 = y * y
    for ell in range(3):
        dcross = sum(
            (
                dq[i].g[ell] * defect.g[i]
                + dq[i].v * defect.h[i][ell]
                for i in range(3)
            ),
            acb(0),
        )
        grad.append(
            lapq.g[ell] * defect.v
            + lapq.v * defect.g[ell]
            - 2 * (y2.g[ell] * cross + y2.v * dcross)
        )
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


def projected_floor_affine_diagnostic(
    ev: DifferentialTrial,
    odd: list[acb],
    epsilon: arb,
    subdivision: tuple[int, int, int],
) -> dict[str, Any]:
    """Center-affine cell model; diagnostic only because no remainder is added."""
    nx1, nx2, ns = subdivision
    xe = edges(-arb(1) / 2, arb(1) / 2, nx1)
    te = edges(arb(0), arb(1) / 2, nx2)
    se = edges(arb(0), epsilon, ns)
    total2 = arb(0)
    affine_l2_2 = arb(0)
    sqrt6 = arb(6).sqrt()
    max_affine = arb(0)
    for i in range(nx1):
        for j in range(nx2):
            for k in range(ns):
                xa, xb = xe[i], xe[i + 1]
                ta, tb = te[j], te[j + 1]
                sa, sb = se[k], se[k + 1]
                xm, tm, sm = (xa + xb) / 2, (ta + tb) / 2, (sa + sb) / 2
                center, grad = floor_residual_param_jet(
                    ev, odd, xm, tm, sm, epsilon
                )
                affine_sup = abs(center).upper() + (
                    (xb - xa) * abs(grad[0]).upper() / 2
                    + (tb - ta) * abs(grad[1]).upper() / 2
                    + (sb - sa) * abs(grad[2]).upper() / 2
                )
                vector_sup = sqrt6 * affine_sup
                sbox = interval(sa, sb)
                xbox = interval(xa, xb)
                tbox = interval(ta, tb)
                rho = sbox.exp()
                y2 = rho - square_nonnegative(xbox) - square_nonnegative(tbox)
                if not bool(y2 > 0):
                    raise ArithmeticError(f"affine floor box height not positive: {i,j,k}")
                coord_cell = (xb - xa) * (tb - ta) * (sb - sa)
                measure = rho.upper() / (2 * y2.lower() * y2.lower())
                total2 += coord_cell * measure * vector_sup * vector_sup
                hx, ht, hs = (xb - xa) / 2, (tb - ta) / 2, (sb - sa) / 2
                # Exact unweighted integral of the squared affine polynomial
                # on a centered rectangular cell; odd cross terms vanish.
                affine_mean_square = (
                    abs(center).upper() ** 2
                    + hx * hx * abs(grad[0]).upper() ** 2 / 3
                    + ht * ht * abs(grad[1]).upper() ** 2 / 3
                    + hs * hs * abs(grad[2]).upper() ** 2 / 3
                )
                affine_l2_2 += 6 * coord_cell * measure * affine_mean_square
                max_affine = max_affine.union(vector_sup).upper()
    norm = total2.sqrt().upper()
    affine_l2 = affine_l2_2.sqrt().upper()
    return {
        "certified": False,
        "load_bearing": False,
        "method": "center value plus center gradient on each cell; second-order remainder omitted",
        "subdivision": list(subdivision),
        "epsilon": str(epsilon),
        "affine_cell_sup_norm_diagnostic": upper(norm),
        "affine_polynomial_L2_norm_diagnostic": upper(affine_l2),
        "max_affine_cell_value": upper(max_affine),
        "target_R_upper": upper(TARGET_R),
        "diagnostic_below_target": bool(norm < TARGET_R.lower()),
        "affine_L2_diagnostic_below_target": bool(affine_l2 < TARGET_R.lower()),
    }


def projected_floor_interval_upper(
    ev: DifferentialTrial,
    odd: list[acb],
    epsilon: arb,
    subdivision: tuple[int, int, int],
) -> dict[str, Any]:
    """Rigorous box upper bound on the projected-core floor residual."""
    nx1, nx2, ns = subdivision
    mat = parse_matrix([[[0, 0], [-1, 0]], [[1, 0], [0, 0]]])
    xe = edges(-arb(1) / 2, arb(1) / 2, nx1)
    te = edges(arb(0), arb(1) / 2, nx2)
    se = edges(arb(0), epsilon, ns)
    actual2 = arb(0)
    requested2 = arb(0)
    value2 = arb(0)
    gradient2 = arb(0)
    sqrt6 = arb(6).sqrt()
    max_actual = arb(0)
    for i in range(nx1):
        for j in range(nx2):
            for k in range(ns):
                x1 = interval(xe[i], xe[i + 1])
                x2 = interval(te[j], te[j + 1])
                s = interval(se[k], se[k + 1])
                r2 = s.exp()
                y2 = r2 - square_nonnegative(x1) - square_nonnegative(x2)
                if not bool(y2 > 0):
                    raise ArithmeticError(
                        f"floor interval does not prove positive height: {i,j,k} {y2}"
                    )
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
                dg_e = sum(
                    (abs(q).upper() ** 2 for q in dg), arb(0)
                ).sqrt().upper()
                d1 = sqrt6 * y.upper() * dg_e
                gradq_e2 = sum((square_nonnegative(q) for q in dq), arb(0))
                gradq = y.upper() * gradq_e2.sqrt().upper()
                value = abs(lapq).upper() * d0
                gradient = 2 * gradq * d1
                requested = value + gradient
                actual = sqrt6 * abs(residual).upper()

                coord_cell = (
                    (xe[i + 1] - xe[i])
                    * (te[j + 1] - te[j])
                    * (se[k + 1] - se[k])
                )
                measure = r2.upper() / (2 * y2.lower() * y2.lower())
                cell = coord_cell * measure
                actual2 += cell * actual * actual
                requested2 += cell * requested * requested
                value2 += cell * value * value
                gradient2 += cell * gradient * gradient
                max_actual = max_actual.union(actual).upper()
    return {
        "certified": True,
        "method": "Arb boxes in exact (x1,x2,log(r^2)) coordinates",
        "subdivision": list(subdivision),
        "epsilon": str(epsilon),
        "requested_norm_upper": upper(requested2.sqrt().upper()),
        "commutator_norm_upper": upper(actual2.sqrt().upper()),
        "value_norm_upper": upper(value2.sqrt().upper()),
        "gradient_norm_upper": upper(gradient2.sqrt().upper()),
        "max_box_commutator_upper": upper(max_actual),
        "target_R_upper": upper(TARGET_R),
        "requested_below_target": bool(requested2.sqrt().upper() < TARGET_R.lower()),
        "commutator_below_target": bool(actual2.sqrt().upper() < TARGET_R.lower()),
    }


def projected_floor_taylor_upper(
    ev: DifferentialTrial,
    odd: list[acb],
    epsilon: arb,
    subdivision: tuple[int, int, int],
) -> dict[str, Any]:
    """Rigorous midpoint-plus-gradient upper bound on the actual commutator."""
    nx1, nx2, ns = subdivision
    xe = edges(-arb(1) / 2, arb(1) / 2, nx1)
    te = edges(arb(0), arb(1) / 2, nx2)
    se = edges(arb(0), epsilon, ns)
    total2 = arb(0)
    sqrt6 = arb(6).sqrt()
    max_sup = arb(0)
    jet_crosscheck = None
    for i in range(nx1):
        for j in range(nx2):
            for k in range(ns):
                xa, xb = xe[i], xe[i + 1]
                ta, tb = te[j], te[j + 1]
                sa, sb = se[k], se[k + 1]
                xm, tm, sm = (xa + xb) / 2, (ta + tb) / 2, (sa + sb) / 2
                r2m = sm.exp()
                y2m = r2m - xm * xm - tm * tm
                ym = y2m.sqrt()
                center, _center_grad = floor_residual_param_jet(
                    ev, odd, xm, tm, sm, epsilon
                )
                if jet_crosscheck is None:
                    legacy_center, _legacy_grad = floor_residual_jet(
                        ev, odd, xm, tm, ym, epsilon
                    )
                    mat = parse_matrix([[[0, 0], [-1, 0]], [[1, 0], [0, 0]]])
                    gp, jac = action_jacobian(mat, xm, tm, ym)
                    lv, lg = witness_value_gradient(ev, odd, xm, tm, ym)
                    rv, rg = witness_value_gradient(ev, odd, *gp)
                    pulled = transpose_pullback(jac, rg)
                    dv = lv - rv
                    dg = [lg[q] - pulled[q] for q in range(3)]
                    u = sm / epsilon
                    hp = 30 * u * u * (1 - u) * (1 - u)
                    hpp = 60 * u - 180 * u * u + 120 * u * u * u
                    r2 = sm.exp()
                    qs = hp / epsilon
                    dq = [2 * xm * qs / r2, 2 * tm * qs / r2, 2 * ym * qs / r2]
                    lapq = -4 * y2m * hpp / (r2 * epsilon * epsilon)
                    direct = lapq * dv - 2 * y2m * sum(
                        (dq[q] * dg[q] for q in range(3)), acb(0)
                    )
                    difference = center - direct
                    agrees = difference.real.contains(0) and difference.imag.contains(0)
                    if not agrees:
                        raise ArithmeticError(
                            f"parameter-jet/direct floor commutator disagreement: {difference}"
                        )
                    legacy_difference = center - legacy_center
                    legacy_agrees = (
                        legacy_difference.real.contains(0)
                        and legacy_difference.imag.contains(0)
                    )
                    if not legacy_agrees:
                        raise ArithmeticError(
                            "parameter/legacy jet floor commutator disagreement: "
                            f"{legacy_difference}"
                        )
                    jet_crosscheck = {
                        "point": [str(xm), str(tm), str(sm)],
                        "parameter_vs_direct_difference_ball": str(difference),
                        "parameter_vs_legacy_jet_difference_ball": str(legacy_difference),
                        "both_contain_zero": True,
                        "load_bearing": True,
                    }

                x = interval(xa, xb)
                t = interval(ta, tb)
                s = interval(sa, sb)
                r2 = s.exp()
                y2 = r2 - square_nonnegative(x) - square_nonnegative(t)
                if not bool(y2 > 0):
                    raise ArithmeticError(f"Taylor floor box height not positive: {i,j,k}")
                _box_value, grad = floor_residual_param_jet(
                    ev, odd, x, t, s, epsilon
                )
                variation = (
                    (xb - xa) * abs(grad[0]).upper() / 2
                    + (tb - ta) * abs(grad[1]).upper() / 2
                    + (sb - sa) * abs(grad[2]).upper() / 2
                )
                scalar_sup = abs(center).upper() + variation
                vector_sup = sqrt6 * scalar_sup
                coord_cell = (xb - xa) * (tb - ta) * (sb - sa)
                measure = r2.upper() / (2 * y2.lower() * y2.lower())
                total2 += coord_cell * measure * vector_sup * vector_sup
                max_sup = max_sup.union(vector_sup).upper()
    norm = total2.sqrt().upper()
    return {
        "certified": True,
        "method": (
            "reflection-orbit Fourier pairing and second-order jets in exact "
            "(x1,x2,log(r^2)) coordinates with first-order cell Taylor fill-in"
        ),
        "subdivision": list(subdivision),
        "epsilon": str(epsilon),
        "commutator_norm_upper": upper(norm),
        "max_cell_commutator_upper": upper(max_sup),
        "target_R_upper": upper(TARGET_R),
        "below_target": bool(norm < TARGET_R.lower()),
        "jet_direct_crosscheck": jet_crosscheck,
    }


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
    ap.add_argument("--certify-floor", action="store_true")
    ap.add_argument("--diagnostic-floor-model", action="store_true")
    ap.add_argument("--floor-model", choices=("taylor",), default="taylor")
    ap.add_argument("--floor-degree", type=int, default=8)
    ap.add_argument("--floor-tensor-degree", default="")
    ap.add_argument("--floor-arithmetic-audit", action="store_true")
    ap.add_argument("--floor-audit-jsonl", type=Path, default=None)
    ap.add_argument("--check-refinement", default="")
    ap.add_argument("--check-degree", type=int, default=0)
    ap.add_argument("--floor-workers", type=int, default=1)
    ap.add_argument(
        "--mass-certificate", type=Path,
        default=ROOT / "track_b_projected_mass_arb_result.json",
    )
    ap.add_argument(
        "--partition-certificate", type=Path,
        default=ROOT / "track_b_partition_result.json",
    )
    ap.add_argument("--spectral-width-tol", default="0.1")
    ap.add_argument("--certify-floor-width", default="")
    ap.add_argument("--certify-floor-subdivision", default="8,4,8")
    ap.add_argument("--taylor-floor-width", default="")
    ap.add_argument("--taylor-floor-subdivision", default="8,4,8")
    ap.add_argument("--affine-floor-width", default="")
    ap.add_argument("--affine-floor-subdivision", default="8,4,8")
    ap.add_argument("--json-out", type=Path, default=ROOT / "track_b_direct_weighted_result.json")
    ns = ap.parse_args()
    ctx.prec = max(128, ns.bits)
    data, row = parse_trial(ns.trial)
    ev = DifferentialTrial(data["parameters"]["M"], str(row["r"]), row["coefficients"])
    odd = projected_coefficients(ev)["odd"]
    if ns.certify_floor or ns.diagnostic_floor_model:
        width_values = [q.strip() for q in ns.floor_widths.split(",") if q.strip()]
        if len(width_values) != 1:
            raise ValueError("floor model mode requires exactly one --floor-widths value")
        width = arb(width_values[0])
        if not (width - arb("0.30")).contains(0):
            raise ValueError("the current floor certificate is restricted to width 0.30")
        floor_dims = tuple(int(q) for q in ns.floor_subdivision.split(","))
        if len(floor_dims) != 3 or min(floor_dims) <= 0:
            raise ValueError("floor subdivision must be nx1,nx2,ns")
        if ns.floor_workers <= 0:
            raise ValueError("floor workers must be positive")
        mass_certificate = json.loads(ns.mass_certificate.read_text(encoding="utf-8"))
        partition_certificate = json.loads(
            ns.partition_certificate.read_text(encoding="utf-8")
        )
        allowed_budget = certified_allowed_floor_budget(
            mass_certificate, arb(ns.spectral_width_tol)
        )
        model_degree: int | tuple[int, int, int] = ns.floor_degree
        if ns.floor_tensor_degree.strip():
            tensor = tuple(int(q) for q in ns.floor_tensor_degree.split(","))
            if len(tensor) != 3:
                raise ValueError("floor tensor degree must be px,py,ps")
            model_degree = tensor
        selected = certify_projected_floor_taylor_grid(
            ev,
            odd,
            ns.trial,
            width,
            floor_dims,
            model_degree,
            int(ctx.prec),
            allowed_budget,
            mass_certificate,
            partition_certificate,
            ns.floor_audit_jsonl,
            ns.floor_workers,
            ns.floor_arithmetic_audit,
        )
        checks = []
        if ns.check_refinement.strip():
            check_dims = tuple(int(q) for q in ns.check_refinement.split(","))
            if len(check_dims) != 3 or min(check_dims) <= 0:
                raise ValueError("check refinement must be nx1,nx2,ns")
            checks.append(certify_projected_floor_taylor_grid(
                ev, odd, ns.trial, width, check_dims, model_degree,
                int(ctx.prec), allowed_budget, mass_certificate,
                partition_certificate, None, ns.floor_workers,
            ))
        if ns.check_degree:
            checks.append(certify_projected_floor_taylor_grid(
                ev, odd, ns.trial, width, floor_dims, ns.check_degree,
                int(ctx.prec), allowed_budget, mass_certificate,
                partition_certificate, None, ns.floor_workers,
            ))

        stability_records = []
        selected_upper = arb(selected["floor_l2_upper"]).upper()
        budget_lower = arb(selected["allowed_budget_lower"]).lower()
        for check in checks:
            check_upper = arb(check["floor_l2_upper"]).upper()
            maximum = max(selected_upper, check_upper)
            tolerance = arb("0.1") * (budget_lower - maximum)
            difference = abs(check_upper - selected_upper).upper()
            passed = bool(tolerance > 0 and difference <= tolerance)
            stability_records.append({
                "grid": check["grid"],
                "polynomial_degree": check["polynomial_degree"],
                "floor_l2_upper": check["floor_l2_upper"],
                "difference_upper": str(difference),
                "allowed_difference_upper": str(tolerance),
                "passed": passed,
            })
        stability_passed = bool(stability_records and all(
            q["passed"] for q in stability_records
        ))
        selected["stability_checks"] = stability_records
        selected["stability_check_passed"] = stability_passed
        selected["provisional"] = not stability_passed
        # Preserve load-bearing trial provenance in the standalone floor
        # artifact.  Earlier artifacts dropped these fields when extracting
        # the Taylor certificate from the enclosing diagnostic result.
        selected["trial"] = str(ns.trial.resolve())
        selected["trial_sha256"] = hashlib.sha256(ns.trial.read_bytes()).hexdigest()
        selected["spectral_parameter"] = str(row["r"])
        if ns.diagnostic_floor_model or not stability_passed:
            selected["floor_residual_certified"] = False
            selected["label"] = "DIAGNOSTIC"
        selected["rung4_certified"] = False
        selected["global_label"] = "GLOBAL RUNG-4 CERTIFICATE: NOT CERTIFIED"
        ns.json_out.write_text(
            json.dumps(selected, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(json.dumps(selected, indent=2, allow_nan=False))
        print(ns.json_out.resolve())
        return 0 if selected["floor_residual_certified"] else 2

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
    certified_floor = None
    if ns.certify_floor_width.strip():
        cert_dims = tuple(int(q) for q in ns.certify_floor_subdivision.split(","))
        if len(cert_dims) != 3 or min(cert_dims) <= 0:
            raise ValueError("certify floor subdivision must be nx1,nx2,ns")
        print(
            f"certified projected floor epsilon={ns.certify_floor_width} boxes={cert_dims}",
            flush=True,
        )
        certified_floor = projected_floor_interval_upper(
            ev, odd, arb(ns.certify_floor_width), cert_dims
        )
    taylor_floor = None
    if ns.taylor_floor_width.strip():
        taylor_dims = tuple(int(q) for q in ns.taylor_floor_subdivision.split(","))
        if len(taylor_dims) != 3 or min(taylor_dims) <= 0:
            raise ValueError("Taylor floor subdivision must be nx1,nx2,ns")
        print(
            f"Taylor projected floor epsilon={ns.taylor_floor_width} boxes={taylor_dims}",
            flush=True,
        )
        taylor_floor = projected_floor_taylor_upper(
            ev, odd, arb(ns.taylor_floor_width), taylor_dims
        )
    affine_floor = None
    if ns.affine_floor_width.strip():
        affine_dims = tuple(int(q) for q in ns.affine_floor_subdivision.split(","))
        if len(affine_dims) != 3 or min(affine_dims) <= 0:
            raise ValueError("affine floor subdivision must be nx1,nx2,ns")
        print(
            f"affine projected floor epsilon={ns.affine_floor_width} boxes={affine_dims}",
            flush=True,
        )
        affine_floor = projected_floor_affine_diagnostic(
            ev, odd, arb(ns.affine_floor_width), affine_dims
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
        "projected_floor_interval_certificate": certified_floor,
        "projected_floor_taylor_certificate": taylor_floor,
        "projected_floor_affine_diagnostic": affine_floor,
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
