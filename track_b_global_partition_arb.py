#!/usr/bin/env python3
r"""Fail-closed Arb certificate for the normalized Track-B partition.

The partition consists of the cusp plateau gate and a tensor product of five
complementary core transition gates.  Consequently the raw denominator is
the exact symbolic identity

    chi_B + (1-chi_B) prod_i (a_i + (1-a_i)) = 1.

All analytic derivative enclosures are nevertheless evaluated cellwise in
the original upper-half-space coordinates (x1,x2,y).  The parameter grid is
(x1,x2,s=log(x1^2+x2^2+y^2)); it is used only to obtain a finite rectangular
cover and rigorous y ranges.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

from flint import arb, ctx


ROOT = Path(__file__).resolve().parent
FLOOR_WIDTH = "0.30"
S_CAP = "0.70"
FLOOR_WEIGHT_FORMULA_ID = "track-b-floor-quintic-logrho-width-0.30/v1"
GATE_NAMES = ("x1m", "x1p", "x2m", "x2p", "floor")
VERTEX_DATA = {
    "v_00": ("0", "0", 4),
    "v_0h": ("0", "0.5", 6),
    "v_mh0": ("-0.5", "0", 6),
    "v_ph0": ("0.5", "0", 6),
    "v_mhh": ("-0.5", "0.5", 12),
    "v_phh": ("0.5", "0.5", 12),
}


def canonical_hash(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def bounds_ball(lo: arb, hi: arb) -> arb:
    if bool(lo > hi):
        raise ValueError(f"reversed bounds: {lo}, {hi}")
    # The python-flint two-argument constructor accepts exact scalar data,
    # not Arb objects (passing balls would fold their magnitudes into the
    # radius).  Convert the outward endpoints to exact binary tuples and
    # form their midpoint/radius with a common exponent.
    lm, le = lo.lower().man_exp()
    hm, he = hi.upper().man_exp()
    exponent = min(int(le), int(he))
    left = int(lm) << (int(le) - exponent)
    right = int(hm) << (int(he) - exponent)
    if left > right:
        raise ValueError("outward endpoint conversion reversed the interval")
    return arb((left + right, exponent - 1), (right - left, exponent - 1))


def _abs_upper(value: arb) -> arb:
    return abs(value).upper()


def positive_square(value: arb) -> arb:
    if not bool(value > 0):
        raise ArithmeticError(f"positive square requested for {value}")
    return bounds_ball(value.lower() ** 2, value.upper() ** 2)


class RealJet2:
    """Real interval value, gradient, and Hessian in (x1,x2,y)."""

    def __init__(self, value: arb, gradient: list[arb] | None = None,
                 hessian: list[list[arb]] | None = None):
        self.v = arb(value)
        self.g = gradient if gradient is not None else [arb(0) for _ in range(3)]
        self.h = hessian if hessian is not None else [
            [arb(0) for _ in range(3)] for _ in range(3)
        ]

    @classmethod
    def constant(cls, value: Any) -> "RealJet2":
        return cls(arb(value))

    @classmethod
    def variable(cls, value: arb, axis: int) -> "RealJet2":
        gradient = [arb(0) for _ in range(3)]
        gradient[axis] = arb(1)
        return cls(value, gradient)

    @staticmethod
    def coerce(value: Any) -> "RealJet2":
        return value if isinstance(value, RealJet2) else RealJet2.constant(value)

    def __add__(self, other: Any) -> "RealJet2":
        b = self.coerce(other)
        return RealJet2(
            self.v + b.v,
            [self.g[i] + b.g[i] for i in range(3)],
            [[self.h[i][j] + b.h[i][j] for j in range(3)] for i in range(3)],
        )

    __radd__ = __add__

    def __neg__(self) -> "RealJet2":
        return RealJet2(-self.v, [-q for q in self.g],
                        [[-q for q in row] for row in self.h])

    def __sub__(self, other: Any) -> "RealJet2":
        return self + (-self.coerce(other))

    def __rsub__(self, other: Any) -> "RealJet2":
        return self.coerce(other) - self

    def __mul__(self, other: Any) -> "RealJet2":
        b = self.coerce(other)
        return RealJet2(
            self.v * b.v,
            [self.g[i] * b.v + self.v * b.g[i] for i in range(3)],
            [[
                self.h[i][j] * b.v + self.v * b.h[i][j]
                + self.g[i] * b.g[j] + self.g[j] * b.g[i]
                for j in range(3)
            ] for i in range(3)],
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "RealJet2":
        if self.v.contains(0):
            raise ZeroDivisionError(f"interval reciprocal contains zero: {self.v}")
        inverse = 1 / self.v
        return RealJet2(
            inverse,
            [-self.g[i] * inverse * inverse for i in range(3)],
            [[
                2 * self.g[i] * self.g[j] * inverse ** 3
                - self.h[i][j] * inverse * inverse
                for j in range(3)
            ] for i in range(3)],
        )

    def __truediv__(self, other: Any) -> "RealJet2":
        return self * self.coerce(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "RealJet2":
        return self.coerce(other) / self

    def log(self) -> "RealJet2":
        if not bool(self.v > 0):
            raise ArithmeticError(f"log argument not positive: {self.v}")
        gradient = [self.g[i] / self.v for i in range(3)]
        hessian = []
        for i in range(3):
            row = []
            for j in range(3):
                first = self.h[i][j] / self.v
                second = self.g[i] * self.g[j] / positive_square(self.v)
                value = first - second
                if not value.is_finite():
                    raise ArithmeticError(
                        "nonfinite log Hessian: "
                        f"v={self.v}, h={self.h[i][j]}, gi={self.g[i]}, "
                        f"gj={self.g[j]}, first={first}, second={second}"
                    )
                row.append(value)
            hessian.append(row)
        return RealJet2(self.v.log(), gradient, hessian)


def quintic_value(t: arb) -> arb:
    return t * t * t * (10 - 15 * t + 6 * t * t)


def smoothstep_clamped(argument: RealJet2) -> RealJet2:
    """C2 quintic clamp with rigorous branch-crossing derivative hulls."""
    if bool(argument.v.upper() <= 0):
        return RealJet2.constant(0)
    if bool(argument.v.lower() >= 1):
        return RealJet2.constant(1)
    lo = max(argument.v.lower(), arb(0))
    hi = min(argument.v.upper(), arb(1))
    t = bounds_ball(lo, hi)
    value = bounds_ball(quintic_value(lo).lower(), quintic_value(hi).upper())
    first = 30 * t * t * (1 - t) * (1 - t)
    second = 60 * t - 180 * t * t + 120 * t * t * t
    return RealJet2(
        value,
        [first * argument.g[i] for i in range(3)],
        [[
            second * argument.g[i] * argument.g[j]
            + first * argument.h[i][j]
            for j in range(3)
        ] for i in range(3)],
    )


def hyperbolic_gradient_upper(jet: RealJet2, y: arb) -> arb:
    return hyperbolic_gradient_norm(jet, y).upper()


def hyperbolic_gradient_norm(jet: RealJet2, y: arb) -> arb:
    """Ball enclosure of the hyperbolic gradient norm, before endpointing."""
    euclidean_sq = sum((q * q for q in jet.g), arb(0))
    if bool(euclidean_sq < 0):
        raise ArithmeticError("negative squared Euclidean gradient enclosure")
    if euclidean_sq.contains(0):
        euclidean_sq = bounds_ball(arb(0), euclidean_sq.upper())
    return abs(y) * euclidean_sq.sqrt()


def hyperbolic_inner(jet_a: RealJet2, jet_b: RealJet2, y: arb) -> arb:
    return (y * y * sum((jet_a.g[i] * jet_b.g[i] for i in range(3)), arb(0)))


def hyperbolic_laplacian(jet: RealJet2, y: arb) -> arb:
    # Project convention: Delta=-y^2(sum d_ii)+y d_y.
    return (-y * y * sum((jet.h[i][i] for i in range(3)), arb(0))
            + y * jet.g[2])


def normalized_laplacian_formula(
    phi: RealJet2, Phi: RealJet2, y: arb
) -> arb:
    """Quotient rule for the project's positive hyperbolic Laplacian.

    Because Delta=-div grad, the final |grad Phi|^2 term has a minus sign:

      Delta(phi/Phi) = Delta phi/Phi - phi Delta Phi/Phi^2
        + 2<grad phi,grad Phi>/Phi^2
        - 2 phi |grad Phi|^2/Phi^3.
    """
    dphi = hyperbolic_laplacian(phi, y)
    dPhi = hyperbolic_laplacian(Phi, y)
    inner = hyperbolic_inner(phi, Phi, y)
    norm_sq = hyperbolic_inner(Phi, Phi, y)
    return (
        dphi / Phi.v - phi.v * dPhi / (Phi.v ** 2)
        + 2 * inner / (Phi.v ** 2)
        - 2 * phi.v * norm_sq / (Phi.v ** 3)
    )


def normalized_gradient_formula(
    phi: RealJet2, Phi: RealJet2
) -> list[arb]:
    """Coordinate gradient of ``phi/Phi`` with the exact quotient rule."""
    return [
        phi.g[i] / Phi.v - phi.v * Phi.g[i] / (Phi.v ** 2)
        for i in range(3)
    ]


def analytic_weight_majorants() -> dict[str, arb]:
    """Resolution-independent rational bounds in the hyperbolic metric."""
    p1 = arb(15) / 8
    p2 = arb(6)
    # On every derivative support y<=6/5.  For a width-1/10 affine
    # horizontal gate, |grad q|_H<=12 and Delta q=0.
    face_grad = p1 * 12
    face_lap = p2 * 12 * 12
    # s=log rho is hyperbolic harmonic and |grad s|_H<=2.
    floor_argument_grad = arb(2) / arb(FLOOR_WIDTH)
    floor_grad = p1 * floor_argument_grad
    floor_lap = p2 * floor_argument_grad * floor_argument_grad
    # q=(y-1.01)/0.19: |grad q|_H, |Delta q| <= (6/5)/(19/100)=120/19.
    cusp_argument = arb(120) / 19
    cusp_grad = p1 * cusp_argument
    cusp_lap = p1 * cusp_argument + p2 * cusp_argument * cusp_argument
    factor_gradients = [cusp_grad] + [face_grad] * 4 + [floor_grad]
    factor_laplacians = [cusp_lap] + [face_lap] * 4 + [floor_lap]
    chamber_grad = sum(factor_gradients, arb(0)).upper()
    cross = arb(0)
    for i in range(len(factor_gradients)):
        for j in range(i + 1, len(factor_gradients)):
            cross += factor_gradients[i] * factor_gradients[j]
    chamber_lap = (sum(factor_laplacians, arb(0)) + 2 * cross).upper()
    return {
        "cusp_gradient": cusp_grad.upper(),
        "cusp_laplacian": cusp_lap.upper(),
        "chamber_gradient": chamber_grad,
        "chamber_laplacian": chamber_lap,
        "face_factor_gradient": face_grad.upper(),
        "face_factor_laplacian": face_lap.upper(),
        "floor_factor_gradient": floor_grad.upper(),
        "floor_factor_laplacian": floor_lap.upper(),
    }


def cell_y_box(
    x_bounds: tuple[arb, arb], t_bounds: tuple[arb, arb],
    s_bounds: tuple[arb, arb]
) -> arb:
    xa, xb = x_bounds
    ta, tb = t_bounds
    sa, sb = s_bounds
    x_abs_min = arb(0) if bool(xa <= 0 and xb >= 0) else min(abs(xa), abs(xb))
    x_abs_max = max(abs(xa), abs(xb))
    t_abs_min = min(abs(ta), abs(tb))
    t_abs_max = max(abs(ta), abs(tb))
    lower_y2 = (sa.exp().lower() - x_abs_max ** 2 - t_abs_max ** 2).lower()
    upper_y2 = (sb.exp().upper() - x_abs_min ** 2 - t_abs_min ** 2).upper()
    if not bool(lower_y2 > 0):
        raise ArithmeticError(f"parameter cell does not prove y^2>0: {lower_y2}")
    return bounds_ball(lower_y2.sqrt().lower(), upper_y2.sqrt().upper())


def log_rho_jet_from_parameter_box(
    x_box: arb, t_box: arb, y_box: arb, s_box: arb
) -> RealJet2:
    rho = bounds_ball(s_box.lower().exp().lower(), s_box.upper().exp().upper())
    rho_sq = positive_square(rho)
    coordinates = (x_box, t_box, y_box)
    gradient = [2 * q / rho for q in coordinates]
    hessian = []
    for i in range(3):
        row = []
        for j in range(3):
            diagonal = 2 / rho if i == j else arb(0)
            row.append(diagonal - 4 * coordinates[i] * coordinates[j] / rho_sq)
        hessian.append(row)
    return RealJet2(s_box, gradient, hessian)


def raw_partition_jets(
    x_box: arb, t_box: arb, y_box: arb, s_box: arb | None = None
) -> tuple[
    list[tuple[str, RealJet2]], dict[str, RealJet2]
]:
    x = RealJet2.variable(x_box, 0)
    t = RealJet2.variable(t_box, 1)
    y = RealJet2.variable(y_box, 2)
    if s_box is None:
        rho = x * x + t * t + y * y
        s = rho.log()
    else:
        s = log_rho_jet_from_parameter_box(x_box, t_box, y_box, s_box)
    cusp = smoothstep_clamped((y - arb("1.01")) / arb("0.19"))
    transitions = {
        "x1m": 1 - smoothstep_clamped((x + arb("0.50")) / arb("0.10")),
        "x1p": smoothstep_clamped((x - arb("0.40")) / arb("0.10")),
        "x2m": 1 - smoothstep_clamped(t / arb("0.10")),
        "x2p": smoothstep_clamped((t - arb("0.40")) / arb("0.10")),
        "floor": 1 - smoothstep_clamped(s / arb(FLOOR_WIDTH)),
    }
    raw: list[tuple[str, RealJet2]] = [("cusp", cusp)]
    core = 1 - cusp
    for bits in itertools.product((0, 1), repeat=len(GATE_NAMES)):
        weight = core
        active = []
        for name, bit in zip(GATE_NAMES, bits):
            gate = transitions[name]
            weight = weight * (gate if bit else 1 - gate)
            if bit:
                active.append(name)
        raw.append(("core:" + ("+".join(active) if active else "identity"), weight))
    return raw, {"cusp": cusp, **transitions, "s": s}


def _jet_is_zero(jet: RealJet2) -> bool:
    return bool(jet.v.is_zero() and all(q.is_zero() for q in jet.g)
                and all(q.is_zero() for row in jet.h for q in row))


def _jet_is_one(jet: RealJet2) -> bool:
    return bool((jet.v - 1).is_zero() and all(q.is_zero() for q in jet.g)
                and all(q.is_zero() for row in jet.h for q in row))


def _cell_record(
    index: tuple[int, int, int], bounds: tuple[tuple[arb, arb], ...],
    degree: int, bits: int, geometry_hash: str, definition_hash: str,
    formulas_hash: str, configuration_hash: str,
    include_weight_details: bool,
) -> dict[str, Any]:
    xb, tb, sb = bounds
    x_box = bounds_ball(*xb)
    t_box = bounds_ball(*tb)
    y_box = cell_y_box(xb, tb, sb)
    raw, gates = raw_partition_jets(x_box, t_box, y_box, bounds_ball(*sb))

    # The tensor-product identity is exact.  Direct interval summation is an
    # independent containment check, not the denominator positivity proof.
    summed = sum((jet for _name, jet in raw), RealJet2.constant(0))
    direct_identity_contains = bool(
        summed.v.contains(1)
        and all(q.contains(0) for q in summed.g)
        and all(q.contains(0) for row in summed.h for q in row)
    )
    Phi = RealJet2.constant(1)
    active = []
    inactive = []
    normalized_one = []
    maximum_gradient = arb(0)
    maximum_laplacian = arb(0)
    worst_gradient = None
    worst_laplacian = None
    b0 = arb(0)
    b1 = arb(0)
    finite = True
    nonfinite_weights = []
    nonfinite_details = []
    direct_interval_dependency_failures = []
    weight_enclosures = []
    majorants = analytic_weight_majorants()
    for name, phi in raw:
        if _jet_is_zero(phi):
            inactive.append(name)
            continue
        active.append(name)
        if _jet_is_one(phi):
            normalized_one.append(name)
        # The direct quotient jet is retained as an independent finiteness
        # check.  The certified derivative bound uses the exact product-rule
        # majorant, avoiding dependency blowup across six factors in [0,1].
        chi = phi / Phi
        direct_grad = hyperbolic_gradient_upper(chi, y_box)
        delta_phi = hyperbolic_laplacian(phi, y_box)
        delta_chi = normalized_laplacian_formula(phi, Phi, y_box)
        direct_lap = _abs_upper(delta_chi)
        if name == "cusp":
            grad = majorants["cusp_gradient"]
            lap = majorants["cusp_laplacian"]
        else:
            grad = majorants["chamber_gradient"]
            lap = majorants["chamber_laplacian"]
        this_finite = grad.is_finite() and lap.is_finite()
        direct_finite = direct_grad.is_finite() and direct_lap.is_finite()
        finite = finite and this_finite
        if not direct_finite:
            direct_interval_dependency_failures.append(name)
        if not this_finite:
            nonfinite_weights.append(name)
            nonfinite_details.append({
                "weight": name,
                "value": str(phi.v),
                "hessian_diagonal": [str(phi.h[i][i]) for i in range(3)],
                "direct_laplacian": str(hyperbolic_laplacian(phi, y_box)),
                "quotient_laplacian": str(normalized_laplacian_formula(phi, Phi, y_box)),
            })
        elif include_weight_details:
            symmetric_lap = bounds_ball(-lap, lap)
            weight_enclosures.append({
                    "weight_id": name,
                    "phi_interval": str(phi.v),
                    "grad_phi_euclidean_components": [str(q) for q in phi.g],
                    "grad_phi_hyperbolic_certified_upper": str(grad),
                    "Delta_phi_certified_interval": str(symmetric_lap),
                    "chi_interval": str(chi.v),
                    "grad_chi_euclidean_components": [str(q) for q in chi.g],
                    "grad_chi_hyperbolic_certified_upper": str(grad),
                    "Delta_chi_certified_interval": str(symmetric_lap),
                    "abs_Delta_chi_certified_upper": str(lap),
                    "direct_interval_diagnostic_finite": direct_finite,
                    "grad_chi_hyperbolic_direct_upper": str(direct_grad),
                    "Delta_phi_direct_interval": str(delta_phi),
                    "Delta_chi_direct_interval": str(delta_chi),
                    "abs_Delta_chi_direct_upper": str(direct_lap),
                    "coordinate_system": "(x1,x2,y)",
                    "metric": "upper-half-space hyperbolic",
            })
        if worst_gradient is None or bool(grad > maximum_gradient):
            maximum_gradient, worst_gradient = grad.upper(), name
        if worst_laplacian is None or bool(lap > maximum_laplacian):
            maximum_laplacian, worst_laplacian = lap.upper(), name

    floor_local = smoothstep_clamped(RealJet2.constant(bounds_ball(
        sb[0] / arb(FLOOR_WIDTH), sb[1] / arb(FLOOR_WIDTH)
    )))
    floor_global_core = 1 - gates["floor"]
    floor_difference = floor_global_core.v - floor_local.v
    floor_consistent = bool(floor_difference.contains(0))

    region_types = []
    if _jet_is_one(gates["cusp"]):
        region_types.append("cusp_plateau")
    elif _jet_is_zero(gates["cusp"]):
        region_types.append("core")
    else:
        region_types.append("cusp_transition_collar")
    for name in GATE_NAMES:
        if not _jet_is_zero(gates[name]):
            region_types.append(name + "_collar")
    def contains_coordinate(interval: tuple[arb, arb], value: str) -> bool:
        q = arb(value)
        return bool(interval[0] <= q and q <= interval[1])

    touches_floor = contains_coordinate(sb, "0")
    x_zero = contains_coordinate(xb, "0")
    x_minus = contains_coordinate(xb, "-0.5")
    x_plus = contains_coordinate(xb, "0.5")
    t_zero = contains_coordinate(tb, "0")
    t_half = contains_coordinate(tb, "0.5")
    stratum_ids = []
    if touches_floor:
        if x_zero: stratum_ids.append("floor_S")
        if x_minus: stratum_ids.append("floor_x1m")
        if x_plus: stratum_ids.append("floor_x1p")
        if t_zero: stratum_ids.append("floor_x2m")
        if t_half: stratum_ids.append("floor_x2p")
    if x_zero and t_zero: stratum_ids.append("vertical_R")
    if x_zero and t_half: stratum_ids.append("vertical_TiR")
    if x_minus and t_zero: stratum_ids.append("vertical_x1m_x2m")
    if x_plus and t_zero: stratum_ids.append("vertical_x1p_x2m")
    if x_minus and t_half: stratum_ids.append("vertical_x1m_x2p")
    if x_plus and t_half: stratum_ids.append("vertical_x1p_x2p")
    if stratum_ids:
        region_types.append("elliptic_edge_collar")
    stratum_orders = {
        "floor_S": 2, "vertical_R": 2, "vertical_TiR": 2,
        "vertical_x1m_x2m": 2, "vertical_x1p_x2m": 2,
        "vertical_x1m_x2p": 2, "vertical_x1p_x2p": 2,
        "floor_x1m": 3, "floor_x1p": 3,
        "floor_x2m": 2, "floor_x2p": 3,
    }
    vertex_ids = []
    if touches_floor:
        for vertex_id, (vx, vt, _order) in VERTEX_DATA.items():
            if contains_coordinate(xb, vx) and contains_coordinate(tb, vt):
                vertex_ids.append(vertex_id)
    if vertex_ids:
        region_types.append("elliptic_vertex_neighborhood")
    stabilizer_order = max(
        [1]
        + [stratum_orders[q] for q in stratum_ids]
        + [VERTEX_DATA[q][2] for q in vertex_ids]
    )

    # Reynolds averaging can activate a chart label that is zero in the
    # canonical representative.  The theorem-facing commutator constants
    # therefore sum the global majorant over the entire 33-weight family.
    # This is intentionally conservative and remains valid on all strata.
    b0 = (majorants["cusp_laplacian"]
          + 32 * majorants["chamber_laplacian"]).upper()
    b1 = (2 * majorants["cusp_gradient"]
          + 64 * majorants["chamber_gradient"]).upper()

    return {
        "cell_id": f"{index[0]},{index[1]},{index[2]}",
        "parent_id": None,
        "refinement_depth": 0,
        "cell_index": list(index),
        "coordinate_bounds": {
            "x1": [str(xb[0]), str(xb[1])],
            "x2": [str(tb[0]), str(tb[1])],
            "s_log_rho": [str(sb[0]), str(sb[1])],
            "y_enclosure": str(y_box),
        },
        "coordinate_system": "parameter cells (x1,x2,s=log rho); derivatives in (x1,x2,y)",
        "metric_convention": "|grad f|_H^2=y^2 sum_i f_i^2; Delta=-y^2 sum_i f_ii+y f_y",
        "active_raw_gates": active,
        "inactive_raw_gates_proved_zero": inactive,
        "raw_gates_proved_one": normalized_one,
        "active_normalized_weights": active,
        "geometry_region_type": region_types,
        "active_chart_transitions": [name for name in GATE_NAMES
                                     if not _jet_is_zero(gates[name])],
        "elliptic_edge_ids": stratum_ids,
        "elliptic_vertex_ids": vertex_ids,
        "stabilizer_order": stabilizer_order,
        "normalization_multiplicity": stabilizer_order,
        "stabilizer_weight_construction": (
            "common exact Reynolds average over the listed finite point/stratum group; "
            "right multiplication reindexes the sum"
        ),
        "geometry_certificate_dependency": geometry_hash,
        "partition_definition_hash": definition_hash,
        "weight_formulas_hash": formulas_hash,
        "configuration_hash": configuration_hash,
        "Phi_lower": "1",
        "Phi_upper": "1",
        "grad_Phi_euclidean_components": ["0", "0", "0"],
        "Delta_Phi_interval": "0",
        "direct_interval_sum_contains_exact_identity": direct_identity_contains,
        "sum_chi_interval": "1",
        "partition_deviation_upper": "0",
        "maximum_grad_chi_upper": str(maximum_gradient.upper()),
        "maximum_Delta_chi_upper": str(maximum_laplacian.upper()),
        "B0_sum_abs_Delta_upper": str(b0.upper()),
        "B1_twice_sum_grad_upper": str(b1.upper()),
        "worst_gradient_weight": worst_gradient,
        "worst_laplacian_weight": worst_laplacian,
        "derivative_bound_method": "exact factorwise product-rule majorant",
        "primitive_hyperbolic_majorants": {
            name: str(value) for name, value in majorants.items()
        },
        "per_active_weight_enclosures": weight_enclosures,
        "per_active_weight_enclosures_recorded": include_weight_details,
        "floor_weight_difference_interval": str(floor_difference),
        "floor_weight_consistency": floor_consistent,
        "floor_weight_identity_proof": (
            "1-(1-H(s/0.30))=H(s/0.30) by exact ring simplification"
        ),
        "working_precision": bits,
        "Taylor_degree": degree,
        "fallback_count": 0,
        "nonfinite_weights": nonfinite_weights,
        "nonfinite_details": nonfinite_details,
        "direct_interval_dependency_failures": direct_interval_dependency_failures,
        "coverage_flags": {
            "uniform_grid_cell": True,
            "positive_height": True,
            "finite_derivatives": finite,
            "partition_identity": direct_identity_contains,
        },
    }


def _edges(a: arb, b: arb, n: int) -> list[arb]:
    return [a + (b - a) * k / n for k in range(n + 1)]


def _partition_cell_worker(payload: tuple[Any, ...]) -> dict[str, Any]:
    (i, j, k, nx, nt, ns, degree, bits, geometry_hash, definition_hash,
     formulas_hash, configuration_hash, include_weight_details) = payload
    ctx.prec = bits
    xb = (-arb(1) / 2 + arb(i) / nx, -arb(1) / 2 + arb(i + 1) / nx)
    tb = (arb(j) / (2 * nt), arb(j + 1) / (2 * nt))
    sb = (arb(S_CAP) * k / ns, arb(S_CAP) * (k + 1) / ns)
    return _cell_record(
        (i, j, k), (xb, tb, sb), degree, bits, geometry_hash, definition_hash,
        formulas_hash, configuration_hash, include_weight_details,
    )


def _stabilizer_certificate(geometry: dict[str, Any]) -> dict[str, Any]:
    incidence = geometry.get("floor_collar_incidence", {})
    edges = incidence.get("elliptic_edge_averages", [])
    edge_orders = [int(q.get("stabilizer_order", 0)) for q in edges]
    vertices = geometry.get("stabilizers", {})
    complete_vertices = bool(vertices) and all(
        q.get("complete", False) and len(q.get("elements", [])) == q.get("order_found")
        for q in vertices.values()
    )
    global_incidence = geometry.get("global_singular_strata_incidence", {})
    global_rows = global_incidence.get("one_dimensional_strata", [])
    global_groups_exact = bool(global_rows) and all(
        q.get("group_closure_exact", False)
        and q.get("right_reindexing_exact", False)
        and len(q.get("elements", [])) == int(q.get("stabilizer_order", 0))
        for q in global_rows
    )
    # The exact vertex groups must also be closed and invariant under right
    # multiplication.  Recheck the serialized matrices independently of the
    # enumeration routine that produced them.
    from track_b_partition_geometry import mmul

    def decode(matrix: list[list[list[int]]]):
        return tuple(tuple(int(v) for v in entry)
                     for row in matrix for entry in row)

    vertex_reindexing: dict[str, bool] = {}
    for name, row in vertices.items():
        group = {decode(q["matrix"]) for q in row.get("elements", [])}
        vertex_reindexing[name] = bool(
            group
            and all(mmul(a, b) in group for a in group for b in group)
            and all({mmul(a, h) for a in group} == group for h in group)
        )
    certified = bool(
        incidence.get("certified", False)
        and edge_orders == [3, 3, 2, 3]
        and complete_vertices
        and global_incidence.get("certified", False)
        and global_incidence.get("face_inventory_complete", False)
        and global_groups_exact
        and vertex_reindexing
        and all(vertex_reindexing.values())
    )
    return {
        "certified": certified,
        "edge_orders": edge_orders,
        "vertex_orders": {name: row.get("order_found") for name, row in vertices.items()},
        "global_stratum_orders": {
            row["id"]: row["stabilizer_order"] for row in global_rows
        },
        "global_stratum_group_reindexing_exact": global_groups_exact,
        "vertex_group_reindexing_exact": vertex_reindexing,
        "averaging_identity": "A_G phi=|G|^-1 sum_(g in G) phi o g",
        "invariance_proof": (
            "exact group reindexing A_G(phi) o h=A_G(phi); isometry chain rule "
            "gives covariant-gradient equivariance and Delta invariance"
        ),
        "normalization_once_proof": (
            "one reference-cell ledger plus Reynolds averaging; group terms are "
            "not duplicated as integration leaves"
        ),
        "sampled_equality_used": False,
    }


def partition_definition() -> dict[str, Any]:
    return {
        "schema": "track-b-global-normalized-partition/v1",
        "cusp": "H((y-1.01)/0.19)",
        "core_factor": "1-cusp",
        "transitions": {
            "x1m": "1-H((x1+0.50)/0.10)",
            "x1p": "H((x1-0.40)/0.10)",
            "x2m": "1-H(x2/0.10)",
            "x2p": "H((x2-0.40)/0.10)",
            "floor": "1-H(log(x1^2+x2^2+y^2)/0.30)",
        },
        "H": "0(t<=0), 10t^3-15t^4+6t^5(0<t<1), 1(t>=1)",
        "raw_weights": "cusp and (1-cusp)*tensor product of a_i or (1-a_i)",
        "denominator_identity": "cusp+(1-cusp)*prod_i(a_i+(1-a_i))=1",
        "chart_label_rule": (
            "five-bit face subset in fixed order x1m,x1p,x2m,x2p,floor; "
            "face matrices and six-copy permutations come from the exact geometry ledger"
        ),
        "singular_stratum_rule": (
            "apply the common Reynolds operator |G|^-1 sum_(g in G) pullback_g "
            "to the entire labeled raw family and quotient chart labels by the exact G action"
        ),
        "floor_weight_formula_id": FLOOR_WEIGHT_FORMULA_ID,
    }


def certification_decision(conditions: dict[str, bool]) -> dict[str, bool]:
    partition_required = (
        "coverage_certified", "denominator_positive_certified",
        "partition_sum_certified", "stabilizer_averages_certified",
        "floor_weight_consistency_certified",
    )
    global_partition = all(bool(conditions.get(name, False))
                           for name in partition_required)
    global_weights = bool(
        global_partition
        and conditions.get("weight_gradients_certified", False)
        and conditions.get("weight_laplacians_certified", False)
        and conditions.get("fallback_zero", False)
    )
    return {
        "global_partition_certified": global_partition,
        "global_weight_bounds_certified": global_weights,
        "rung4_certified": False,
    }


def certify_partition(
    geometry: dict[str, Any], subdivision: tuple[int, int, int], degree: int,
    bits: int, audit_path: Path | None, workers: int = 1,
    include_weight_details: bool = True,
) -> dict[str, Any]:
    if degree < 5:
        raise ValueError("partition degree must be at least the exact quintic degree 5")
    nx, nt, ns = subdivision
    definition = partition_definition()
    definition_hash = canonical_hash(definition)
    geometry_payload = {
        "coordinate_convention": geometry.get("coordinate_convention"),
        "faces": geometry.get("faces"),
        "floor_collar_incidence": geometry.get("floor_collar_incidence"),
        "global_singular_strata_incidence": geometry.get(
            "global_singular_strata_incidence"
        ),
        "stabilizers": geometry.get("stabilizers"),
    }
    geometry_hash = canonical_hash(geometry_payload)
    formulas_payload = {
        "definition": definition,
        "gradient_rule": "gphi/Phi-phi*gPhi/Phi^2",
        "laplacian": "-y^2(dxx+dtt+dyy)+y*dy",
        "quotient_rule": (
            "dphi/Phi-phi*dPhi/Phi^2+2<gphi,gPhi>/Phi^2"
            "-2phi|gPhi|^2/Phi^3"
        ),
        "majorants": {k: str(v) for k, v in analytic_weight_majorants().items()},
    }
    formulas_hash = canonical_hash(formulas_payload)
    config = {
        "subdivision": list(subdivision), "Taylor_degree": degree,
        "arb_bits": bits, "s_cap": S_CAP,
    }
    configuration_hash = canonical_hash(config)
    payloads = [
        (i, j, k, nx, nt, ns, degree, bits, geometry_hash, definition_hash,
         formulas_hash, configuration_hash, include_weight_details)
        for i in range(nx) for j in range(nt) for k in range(ns)
    ]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            records = list(pool.map(_partition_cell_worker, payloads, chunksize=4))
    else:
        records = [_partition_cell_worker(payload) for payload in payloads]
    expected_count = nx * nt * ns
    ids = [q["cell_id"] for q in records]
    exact_once = len(records) == expected_count and len(ids) == len(set(ids))
    minimum_denominator = arb(1)
    worst_denominator = records[0]["cell_id"] if records else None
    maximum_gradient = arb(0)
    maximum_laplacian = arb(0)
    maximum_b0 = arb(0)
    maximum_b1 = arb(0)
    worst_gradient_cell = worst_gradient_weight = None
    worst_laplacian_cell = worst_laplacian_weight = None
    fallback_count = 0
    per_cell_ok = True
    floor_consistency = True
    for record in records:
        grad = arb(record["maximum_grad_chi_upper"]).upper()
        lap = arb(record["maximum_Delta_chi_upper"]).upper()
        b0 = arb(record["B0_sum_abs_Delta_upper"]).upper()
        b1 = arb(record["B1_twice_sum_grad_upper"]).upper()
        if bool(grad > maximum_gradient):
            maximum_gradient = grad
            worst_gradient_cell = record["cell_id"]
            worst_gradient_weight = record["worst_gradient_weight"]
        if bool(lap > maximum_laplacian):
            maximum_laplacian = lap
            worst_laplacian_cell = record["cell_id"]
            worst_laplacian_weight = record["worst_laplacian_weight"]
        maximum_b0 = max(maximum_b0, b0)
        maximum_b1 = max(maximum_b1, b1)
        fallback_count += int(record["fallback_count"])
        per_cell_ok = per_cell_ok and all(record["coverage_flags"].values())
        floor_consistency = floor_consistency and bool(record["floor_weight_consistency"])
    if audit_path is not None:
        canonical_audit = "".join(
            json.dumps(q, sort_keys=True) + "\n" for q in records
        )
        audit_path.write_text(canonical_audit, encoding="utf-8")
        audit_sha256 = hashlib.sha256(audit_path.read_bytes()).hexdigest()
        audit_records_sha256 = hashlib.sha256(
            canonical_audit.encode("utf-8")
        ).hexdigest()
    else:
        audit_sha256 = None
        audit_records_sha256 = None

    # If y<=1.20 then rho<=1/2+(6/5)^2=1.94<exp(0.70), so the
    # rectangular s-cover contains the entire nonplateau quotient.
    cap_coverage = bool(arb(S_CAP).exp().lower() > arb("1.94"))
    floor_geometry_incidence = bool(
        geometry.get("floor_collar_incidence", {}).get("certified", False)
    )
    global_geometry = geometry.get("global_singular_strata_incidence", {})
    geometry_incidence = bool(
        floor_geometry_incidence
        and global_geometry.get("certified", False)
        and global_geometry.get("face_inventory_complete", False)
        and geometry.get("stabilizer_flags", {}).get(
            "global_singular_strata_exhaustive", False
        )
    )
    stabilizers = _stabilizer_certificate(geometry)
    coverage = bool(exact_once and cap_coverage and per_cell_ok and geometry_incidence)
    denominator_positive = bool(minimum_denominator.lower() > 0)
    partition_sum = bool(per_cell_ok)
    gradients = bool(per_cell_ok and maximum_gradient.is_finite())
    laplacians = bool(per_cell_ok and maximum_laplacian.is_finite())
    decision = certification_decision({
        "coverage_certified": coverage,
        "denominator_positive_certified": denominator_positive,
        "partition_sum_certified": partition_sum,
        "stabilizer_averages_certified": stabilizers["certified"],
        "floor_weight_consistency_certified": floor_consistency,
        "weight_gradients_certified": gradients,
        "weight_laplacians_certified": laplacians,
        "fallback_zero": fallback_count == 0,
    })
    global_partition = decision["global_partition_certified"]
    global_weights = decision["global_weight_bounds_certified"]
    volume_sqrt_upper = (arb(1) / 2).sqrt().upper()
    b0_l2 = (maximum_b0 * volume_sqrt_upper).upper()
    b1_l2 = (maximum_b1 * volume_sqrt_upper).upper()
    return {
        "label": "GLOBAL PARTITION CERTIFICATE",
        "schema": "track-b-global-partition-certificate/v1",
        "certified": global_partition,
        "arb_bits": bits,
        "cell_count": len(records),
        "partition_degree": degree,
        "subdivision": list(subdivision),
        "minimum_denominator_lower": str(minimum_denominator.lower()),
        "worst_denominator_cell": worst_denominator,
        "maximum_gradient_upper": str(maximum_gradient.upper()),
        "worst_gradient_cell": worst_gradient_cell,
        "worst_gradient_weight": worst_gradient_weight,
        "maximum_laplacian_upper": str(maximum_laplacian.upper()),
        "worst_laplacian_cell": worst_laplacian_cell,
        "worst_laplacian_weight": worst_laplacian_weight,
        "B0_sup_upper": str(maximum_b0.upper()),
        "B1_sup_upper": str(maximum_b1.upper()),
        "b0_upper": str(b0_l2),
        "b1_upper": str(b1_l2),
        "partition_constants_certified": global_weights,
        "coverage_certified": coverage,
        "local_finiteness_certified": bool(exact_once and per_cell_ok),
        "denominator_positive_certified": denominator_positive,
        "partition_sum_certified": partition_sum,
        "partition_deviation_upper": "0",
        "stabilizer_averages_certified": stabilizers["certified"],
        "stabilizers_certified": stabilizers["certified"],
        "weight_gradients_certified": gradients,
        "weight_laplacians_certified": laplacians,
        "floor_weight_consistency_certified": floor_consistency,
        "geometry_incidence_certified": geometry_incidence,
        "floor_geometry_incidence_certified": floor_geometry_incidence,
        "global_singular_strata_exhaustive": bool(
            global_geometry.get("certified", False)
        ),
        "global_partition_certified": global_partition,
        "global_weight_bounds_certified": global_weights,
        "unresolved_fallback_count": fallback_count,
        "transitions_complete": bool(exact_once and geometry_incidence),
        "transition_set_complete": bool(exact_once and geometry_incidence),
        "active_transition_ids": list(GATE_NAMES),
        "theorem_DK_compatible": global_weights,
        "two_cusp_coordinates_certified": bool(
            geometry.get("two_cusp_coordinates_certified", False)
        ),
        "stability_check_passed": False,
        "provisional": True,
        "partition_definition": definition,
        "partition_definition_hash": definition_hash,
        "weight_formulas_hash": formulas_hash,
        "geometry_incidence_hash": geometry_hash,
        "floor_geometry_incidence_hash": canonical_hash(
            geometry.get("floor_collar_incidence")
        ),
        "configuration_hash": configuration_hash,
        "floor_width": FLOOR_WIDTH,
        "floor_weight_formula_id": FLOOR_WEIGHT_FORMULA_ID,
        "stabilizer_certificate": stabilizers,
        "coverage_proof": {
            "uniform_grid_exact_once": exact_once,
            "s_cap": S_CAP,
            "exp_s_cap_lower_gt_1.94": cap_coverage,
            "nonplateau_rho_upper": "1/2+(6/5)^2=1.94",
            "cusp_plateau": "y>=1.20 is one analytic chart with constant weight 1",
            "interfaces": "half-open ownership; closed interfaces retained only for incidence checks",
            "global_singular_strata_ledger": bool(
                global_geometry.get("certified", False)
            ),
            "all_nonempty_face_intersections_classified": bool(
                global_geometry.get("face_inventory_complete", False)
            ),
        },
        "denominator_strategy": {
            "method": "exact complementary-gate simplification before interval evaluation",
            "identity": definition["denominator_identity"],
            "adaptive_positive_clamp_used": False,
            "midpoint_positivity_used": False,
        },
        "coordinate_and_metric_convention": {
            "cell_coordinates": "(x1,x2,s=log rho)",
            "stored_derivatives": "original (x1,x2,y)",
            "gradient_norm": "y*Euclidean gradient norm",
            "laplacian": "Delta=-y^2(d_x1^2+d_x2^2+d_y^2)+y*d_y",
        },
        "audit_jsonl": None if audit_path is None else str(audit_path.resolve()),
        "audit_sha256": audit_sha256,
        "audit_records_sha256": audit_records_sha256,
        "per_active_weight_enclosures_recorded": include_weight_details,
        "rung4_certified": False,
        "status": "GREEN" if global_weights else "RED",
    }


def apply_stability(primary: dict[str, Any], check: dict[str, Any]) -> dict[str, Any]:
    primary_grad = arb(primary["maximum_gradient_upper"]).upper()
    check_grad = arb(check["maximum_gradient_upper"]).upper()
    primary_lap = arb(primary["maximum_laplacian_upper"]).upper()
    check_lap = arb(check["maximum_laplacian_upper"]).upper()
    def relative_change(a: arb, b: arb) -> arb:
        denominator = min(a, b)
        if not bool(denominator > 0):
            if (a - b).contains(0):
                return arb(0)
            raise ArithmeticError("relative stability change has zero denominator")
        return (abs(a - b) / denominator).upper()
    grad_change = relative_change(primary_grad, check_grad)
    lap_change = relative_change(primary_lap, check_lap)
    primary_den = arb(primary["minimum_denominator_lower"]).lower()
    check_den = arb(check["minimum_denominator_lower"]).lower()
    denominator_drop = max(arb(0), (primary_den - check_den) / primary_den).upper()
    stable = bool(
        primary.get("global_weight_bounds_certified", False)
        and check.get("global_weight_bounds_certified", False)
        and grad_change <= arb("0.10")
        and lap_change <= arb("0.10")
        and denominator_drop <= arb("0.10")
        and primary.get("partition_definition_hash") == check.get("partition_definition_hash")
        and primary.get("geometry_incidence_hash") == check.get("geometry_incidence_hash")
    )
    primary["stability_check"] = {
        "check_subdivision": check["subdivision"],
        "check_partition_degree": check["partition_degree"],
        "check_minimum_denominator_lower": check["minimum_denominator_lower"],
        "check_maximum_gradient_upper": check["maximum_gradient_upper"],
        "check_maximum_laplacian_upper": check["maximum_laplacian_upper"],
        "gradient_relative_change_upper": str(grad_change),
        "laplacian_relative_change_upper": str(lap_change),
        "denominator_relative_drop_upper": str(denominator_drop),
        "passed": stable,
    }
    primary["stability_check_passed"] = stable
    primary["provisional"] = not stable
    if not stable:
        primary["certified"] = False
        primary["global_partition_certified"] = False
        primary["global_weight_bounds_certified"] = False
        primary["partition_constants_certified"] = False
        primary["status"] = "YELLOW"
    primary["rung4_certified"] = False
    return primary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", type=Path,
                        default=ROOT / "track_b_partition_result.json")
    parser.add_argument("--bits", type=int, default=192)
    parser.add_argument("--partition-subdivision", default="16,8,16")
    parser.add_argument("--partition-degree", type=int, default=8)
    parser.add_argument("--partition-workers", type=int, default=1)
    parser.add_argument("--partition-audit-jsonl", type=Path,
                        default=ROOT / "track_b_partition_cells.jsonl")
    parser.add_argument("--check-refinement", default="")
    parser.add_argument("--check-degree", type=int, default=0)
    parser.add_argument("--check-audit-jsonl", type=Path, default=None)
    parser.add_argument("--json-out", type=Path,
                        default=ROOT / "track_b_global_partition_result.json")
    ns = parser.parse_args()
    ctx.prec = max(128, ns.bits)
    geometry = json.loads(ns.geometry.read_text(encoding="utf-8"))
    dims = tuple(int(q) for q in ns.partition_subdivision.split(","))
    if len(dims) != 3 or min(dims) <= 0:
        parser.error("partition subdivision must be nx1,nx2,ns")
    primary = certify_partition(
        geometry, dims, ns.partition_degree, int(ctx.prec), ns.partition_audit_jsonl,
        ns.partition_workers, True,
    )
    check = None
    if ns.check_refinement.strip():
        check_dims = tuple(int(q) for q in ns.check_refinement.split(","))
        if len(check_dims) != 3 or min(check_dims) <= 0:
            parser.error("check refinement must be nx1,nx2,ns")
        check = certify_partition(
            geometry, check_dims, ns.partition_degree, int(ctx.prec),
            ns.check_audit_jsonl, ns.partition_workers, False,
        )
    elif ns.check_degree:
        check = certify_partition(
            geometry, dims, ns.check_degree, int(ctx.prec), ns.check_audit_jsonl,
            ns.partition_workers, False,
        )
    if check is not None:
        primary = apply_stability(primary, check)
    primary["rung4_certified"] = False
    ns.json_out.write_text(json.dumps(primary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(primary, indent=2))
    return 0 if primary.get("global_weight_bounds_certified", False) else 2


if __name__ == "__main__":
    raise SystemExit(main())
