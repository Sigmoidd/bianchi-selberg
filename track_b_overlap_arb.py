#!/usr/bin/env python3
"""Cancellation-preserving Arb overlap certificate for Track-B D(K).

The geometry is deliberately an external, strict input.  Each active
transition certifies the six-vector differences

    W(p) - rho(g)^(-1) W(g p),
    nabla W(p) - nabla[rho(g)^(-1) W(g p)],

on a finite collection of interval boxes.  The second expression includes
the exact Jacobian pullback.  Values and derivatives are summed as Acb balls
*before* absolute values are taken, so Fourier-mode cancellation is retained.

This program fails closed if the atlas is incomplete or uncertified.  It does
not invent missing Humbert overlaps, stabilizers, cusp charts, or partition
constants.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from flint import acb, arb, ctx

from continuum_defect_arb import TrialEvaluator, lower, parse_trial, upper


SCHEMA = "track-b-humbert-partition-v1"
PARTITION_LEDGER_SCHEMA = "track_b_partition/v1"
TRANSPORT = (
    "stored inverse right action pi_g: defect[c]=F_c(p)-F_{pi_g(c)}(g p); "
    "covectors pulled back by J_g^T"
)


def _q(v: Any) -> arb:
    return arb(str(v))


def _gaussian_integer(v: Any) -> tuple[int, int]:
    if not isinstance(v, list) or len(v) != 2:
        raise ValueError(f"Gaussian number must be [real,imag], got {v!r}")
    if any(isinstance(q, bool) or not isinstance(q, int) for q in v):
        raise ValueError(f"matrix entry is not an exact Gaussian integer: {v!r}")
    return int(v[0]), int(v[1])


def _z(v: Any) -> acb:
    q = _gaussian_integer(v)
    return acb(q[0], q[1])


def parse_matrix(raw: Any) -> tuple[acb, acb, acb, acb]:
    if not isinstance(raw, list) or len(raw) != 2 or any(
        not isinstance(row, list) or len(row) != 2 for row in raw
    ):
        raise ValueError("matrix must be 2x2 with [real,imag] entries")
    aq, bq = _gaussian_integer(raw[0][0]), _gaussian_integer(raw[0][1])
    cq, dq = _gaussian_integer(raw[1][0]), _gaussian_integer(raw[1][1])

    def mul(u: tuple[int, int], v: tuple[int, int]) -> tuple[int, int]:
        return u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0]

    ad, bc = mul(aq, dq), mul(bq, cq)
    det = ad[0] - bc[0], ad[1] - bc[1]
    if det != (1, 0):
        raise ValueError(f"transition matrix determinant is not 1 exactly: {det}")
    a, b, c, d = acb(*aq), acb(*bq), acb(*cq), acb(*dq)
    return a, b, c, d


def parse_perm(raw: Any) -> tuple[int, ...]:
    out = tuple(int(v) for v in raw)
    if len(out) != 6 or sorted(out) != list(range(6)):
        raise ValueError(f"not a six-copy permutation: {raw!r}")
    return out


def compose_perm(p: tuple[int, ...], q: tuple[int, ...]) -> tuple[int, ...]:
    """Index-map composition: applying q and then p gives q[p[c]]."""
    return tuple(q[p[c]] for c in range(6))


def square_range(x: arb) -> arb:
    ax = abs(x)
    return (ax.lower() * ax.lower()).union(ax.upper() * ax.upper())


def real_part(z: acb) -> arb:
    return z.real


def action_jacobian(
    mat: tuple[acb, acb, acb, acb], x: arb, t: arb, y: arb
) -> tuple[tuple[arb, arb, arb], list[list[arb]]]:
    """Rigorous PSL(2,C) action and coordinate Jacobian on an Arb box."""
    a, b, c, d = mat
    z = acb(x, t)
    q = c * z + d
    p = a * z + b
    c2 = square_range(c.real) + square_range(c.imag)
    den = square_range(q.real) + square_range(q.imag) + c2 * square_range(y)
    if not bool(den > 0):
        raise ArithmeticError(f"action denominator is not positive: {den}")
    num = p * q.conjugate() + a * c.conjugate() * y * y
    zout = num / den
    yout = y / den

    dq = [c, acb(0, 1) * c, acb(0)]
    dp = [a, acb(0, 1) * a, acb(0)]
    dden = [
        2 * real_part(dq[0] * q.conjugate()),
        2 * real_part(dq[1] * q.conjugate()),
        2 * c2 * y,
    ]
    dnum = [
        dp[0] * q.conjugate() + p * dq[0].conjugate(),
        dp[1] * q.conjugate() + p * dq[1].conjugate(),
        2 * a * c.conjugate() * y,
    ]
    dz = [(dnum[k] * den - num * dden[k]) / (den * den) for k in range(3)]
    dy = [
        -y * dden[0] / (den * den),
        -y * dden[1] / (den * den),
        (den - y * dden[2]) / (den * den),
    ]
    jac = [
        [dz[k].real for k in range(3)],
        [dz[k].imag for k in range(3)],
        dy,
    ]
    return (zout.real, zout.imag, yout), jac


def identity_matrix() -> tuple[acb, acb, acb, acb]:
    return acb(1), acb(0), acb(0), acb(1)


def transpose_pullback(jac: list[list[arb]], grad: list[acb]) -> list[acb]:
    return [sum((grad[r] * jac[r][k] for r in range(3)), acb(0)) for k in range(3)]


class DifferentialTrial(TrialEvaluator):
    """TrialEvaluator with exact first coordinate derivatives."""

    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__(*args, **kwargs)
        self._shifted_k_cache: dict[tuple[str, str], acb] = {}

    def _kir_order(self, order: acb, arg: arb) -> acb:
        key = (str(order), self._key_ball(arg))
        cached = self._shifted_k_cache.get(key)
        if cached is not None:
            return cached
        out = acb(arg).bessel_k(order)
        if not out.is_finite():
            raise ArithmeticError(
                "direct Arb K enclosure failed for derivative order; refine geometry box"
            )
        self._shifted_k_cache[key] = out
        return out

    def component_value_gradient(
        self, copy: int, x1: arb, x2: arb, y: arb
    ) -> tuple[acb, list[acb]]:
        if not bool(y > 0):
            raise ArithmeticError(f"nonpositive height box: {y}")
        cusp = 0 if copy == 0 else 1
        modes = self.modes_inf if cusp == 0 else self.modes_0
        offset = 0 if cusp == 0 else self.ni
        tx = arb(0) if cusp == 0 else arb(copy - 1)
        value = acb(0)
        grad = [acb(0), acb(0), acb(0)]
        for k, mode in enumerate(modes):
            u, v, mag = self._frequency(cusp, mode)
            omega = 2 * self.pi * mag
            arg = omega * y
            kval = self._kir(cusp, k, mag, y)
            # K'_nu(x)=-[K_(nu-1)(x)+K_(nu+1)(x)]/2.
            km = self._kir_order(self.order - 1, arg)
            kp = self._kir_order(self.order + 1, arg)
            kprime = -(km + kp) / 2
            phase = acb(0, 2 * self.pi * (u * (x1 + tx) + v * x2)).exp()
            term = self.coeff[offset + k] * phase
            radial = y * kval
            value += term * radial
            grad[0] += term * acb(0, 2 * self.pi * u) * radial
            grad[1] += term * acb(0, 2 * self.pi * v) * radial
            grad[2] += term * (kval + y * omega * kprime)
        return value, grad

    def vector_value_gradient(
        self, x1: arb, x2: arb, y: arb
    ) -> tuple[list[acb], list[list[acb]]]:
        vals, grads = [], []
        for copy in range(6):
            v, g = self.component_value_gradient(copy, x1, x2, y)
            vals.append(v)
            grads.append(g)
        return vals, grads


def average_candidate(
    ev: DifferentialTrial,
    point: tuple[arb, arb, arb],
    stabilizer: list[dict[str, Any]],
) -> tuple[list[acb], list[list[acb]]]:
    """Finite orbifold average, returned in the input-point trivialization."""
    terms = stabilizer or [
        {"matrix": [[[1, 0], [0, 0]], [[0, 0], [1, 0]]], "pi": list(range(6))}
    ]
    vals = [acb(0) for _ in range(6)]
    grads = [[acb(0) for _ in range(3)] for _ in range(6)]
    for item in terms:
        mat = parse_matrix(item["matrix"])
        pi = parse_perm(item["pi"])
        hp, jac = action_jacobian(mat, *point)
        hv, hg = ev.vector_value_gradient(*hp)
        for c in range(6):
            src = pi[c]  # rho(h)^(-1) W(hp), using stored pi_h.
            vals[c] += hv[src]
            pulled = transpose_pullback(jac, hg[src])
            for k in range(3):
                grads[c][k] += pulled[k]
    n = len(terms)
    return [v / n for v in vals], [[q / n for q in g] for g in grads]


def defect_on_box(
    ev: DifferentialTrial, box: list[list[Any]], transition: dict[str, Any]
) -> tuple[arb, arb]:
    p = tuple(_q(bounds[0]).union(_q(bounds[1])) for bounds in box)
    if len(p) != 3 or not bool(p[2] > 0):
        raise ValueError(f"invalid H3 box {box!r}")
    left_v, left_g = average_candidate(ev, p, transition.get("source_stabilizer", []))
    mat = parse_matrix(transition["matrix"])
    pi = parse_perm(transition["pi"])
    gp, jac = action_jacobian(mat, *p)
    right_v_q, right_g_q = average_candidate(
        ev, gp, transition.get("target_stabilizer", [])
    )

    value2 = arb(0)
    grad2 = arb(0)
    y = p[2]
    for c in range(6):
        rc = pi[c]
        dv = left_v[c] - right_v_q[rc]
        value2 += abs(dv).upper() ** 2
        pulled = transpose_pullback(jac, right_g_q[rc])
        for k in range(3):
            dg = left_g[c][k] - pulled[k]
            # Hyperbolic covector norm is y times Euclidean coordinate norm.
            grad2 += y.upper() ** 2 * abs(dg).upper() ** 2
    return value2.sqrt().upper(), grad2.sqrt().upper()


def subdivide(box: list[list[Any]], n: int):
    edges = []
    for lo, hi in box:
        a, b = _q(lo), _q(hi)
        edges.append([a + (b - a) * k / n for k in range(n + 1)])
    for i in range(n):
        for j in range(n):
            for k in range(n):
                yield [
                    [str(edges[0][i]), str(edges[0][i + 1])],
                    [str(edges[1][j]), str(edges[1][j + 1])],
                    [str(edges[2][k]), str(edges[2][k + 1])],
                ]


def validate_geometry(data: dict[str, Any]) -> None:
    required = {
        "schema", "certified", "coverage_certified", "local_finiteness_certified",
        "transition_set_complete", "stabilizers_complete", "two_cusp_coordinates_certified",
        "active_transitions",
    }
    missing = sorted(required - set(data))
    if missing:
        raise ValueError(f"geometry missing required fields: {missing}")
    if data["schema"] != SCHEMA:
        raise ValueError(f"unsupported geometry schema {data['schema']!r}")
    gates = [
        "certified", "coverage_certified", "local_finiteness_certified",
        "transition_set_complete", "stabilizers_complete", "two_cusp_coordinates_certified",
    ]
    failed = [k for k in gates if data.get(k) is not True]
    if failed:
        raise ValueError(f"geometry fails closed; uncertified gates: {failed}")
    ids = set()
    if not data["active_transitions"]:
        raise ValueError("geometry contains no active transitions")
    for tr in data["active_transitions"]:
        for key in ("id", "matrix", "pi", "overlap_boxes", "cusp_coordinate"):
            if key not in tr:
                raise ValueError(f"transition missing {key}: {tr}")
        if tr["id"] in ids:
            raise ValueError(f"duplicate transition id {tr['id']}")
        ids.add(tr["id"])
        parse_matrix(tr["matrix"])
        parse_perm(tr["pi"])
        if tr["cusp_coordinate"] not in ("infinity", "zero", "compact"):
            raise ValueError(f"bad cusp coordinate tag: {tr['cusp_coordinate']}")
        if not tr["overlap_boxes"]:
            raise ValueError(f"transition {tr['id']} has no overlap boxes")


def fail_closed_result(
    geometry: dict[str, Any], geometry_path: Path, trial_path: Path, reason: str
) -> dict[str, Any]:
    blockers = list(geometry.get("blockers", []))
    if reason not in blockers:
        blockers.insert(0, reason)
    return {
        "status": "fail-closed Track-B overlap certificate: geometry incomplete",
        "certified": False,
        "delta0_upper": None,
        "delta1_upper": None,
        "weighted_residual_upper": None,
        "direct_weighted_defect_L2_upper": None,
        "track_B_R_target_upper": 0.0101903405004245,
        "closes_lambda_width_below_0.1": False,
        "active_transition_ids": geometry.get("active_transition_ids", []),
        "all_transitions_covered": False,
        "common_fiber_transport_certified": False,
        "stabilizer_averaging_certified": False,
        "first_gradients_certified": False,
        "theorem_DK_compatible": False,
        "transport_convention": geometry.get("coordinate_convention", {}).get(
            "transport", TRANSPORT
        ),
        "stabilization_certified": False,
        "two_cusp_coordinates_certified": False,
        "refinement_table": [],
        "refinement_monotone": False,
        "independent_check": None,
        "theorem_compatibility": {
            "theorem": "Track-B six-copy/two-cusp D(K)",
            "tau": 0,
            "admissible": False,
            "reason": reason,
        },
        "blockers": blockers,
        "conditional_partition_bounds_not_used": geometry.get(
            "conditional_bounds_not_theorem_inputs"
        ),
        "geometry": str(geometry_path.resolve()),
        "geometry_sha256": hashlib.sha256(geometry_path.read_bytes()).hexdigest(),
        "trial": str(trial_path.resolve()),
        "trial_sha256": hashlib.sha256(trial_path.read_bytes()).hexdigest(),
        "rung4_certified": False,
    }


def midpoint_jacobian_check(geometry: dict[str, Any]) -> dict[str, Any]:
    """Independent floating finite-difference check of the analytic Jacobian."""
    worst = 0.0
    records = []
    for tr in geometry["active_transitions"]:
        box = tr["overlap_boxes"][0]
        p = [0.5 * (float(q[0]) + float(q[1])) for q in box]
        mat = parse_matrix(tr["matrix"])
        _gp, jac = action_jacobian(mat, *[arb(repr(v)) for v in p])
        jf = [[float(jac[r][k].mid()) for k in range(3)] for r in range(3)]
        h = 1e-6
        numeric = [[0.0] * 3 for _ in range(3)]
        for k in range(3):
            pp, pm = p.copy(), p.copy()
            pp[k] += h
            pm[k] -= h
            ap, _ = action_jacobian(mat, *[arb(repr(v)) for v in pp])
            am, _ = action_jacobian(mat, *[arb(repr(v)) for v in pm])
            for r in range(3):
                numeric[r][k] = (float(ap[r].mid()) - float(am[r].mid())) / (2 * h)
                worst = max(worst, abs(numeric[r][k] - jf[r][k]))
        records.append({"transition_id": tr["id"], "max_abs_error": max(
            abs(numeric[r][k] - jf[r][k]) for r in range(3) for k in range(3)
        )})
    return {
        "method": "independent centered finite-difference versus analytic PSL2(C) Jacobian",
        "load_bearing": False,
        "worst_abs_error": worst,
        "passed": worst < 1e-6,
        "records": records,
    }


def certify_level(
    ev: DifferentialTrial, geometry: dict[str, Any], refinement: int
) -> dict[str, Any]:
    d0 = arb(0)
    d1 = arb(0)
    per = []
    n_boxes = 0
    for tr in geometry["active_transitions"]:
        t0, t1 = arb(0), arb(0)
        count = 0
        for base in tr["overlap_boxes"]:
            for box in subdivide(base, refinement):
                q0, q1 = defect_on_box(ev, box, tr)
                t0 = t0.union(q0).upper()
                t1 = t1.union(q1).upper()
                count += 1
        d0 = d0.union(t0).upper()
        d1 = d1.union(t1).upper()
        n_boxes += count
        per.append({
            "transition_id": tr["id"],
            "cusp_coordinate": tr["cusp_coordinate"],
            "boxes": count,
            "delta0_upper": upper(t0),
            "delta1_upper": upper(t1),
        })
    return {
        "refinement": refinement,
        "evaluated_boxes": n_boxes,
        "delta0_upper": upper(d0),
        "delta1_upper": upper(d1),
        "per_transition": per,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", type=Path, default=Path("track_b_partition_result.json"))
    ap.add_argument("--trial", type=Path, default=Path("six_copy_hejhal_balanced_coeffs.json"))
    ap.add_argument("--bits", type=int, default=192)
    ap.add_argument("--refinements", default="1,2")
    ap.add_argument("--json-out", type=Path, default=Path("track_b_overlap_result.json"))
    ns = ap.parse_args()
    ctx.prec = max(128, ns.bits)
    geometry = json.loads(ns.geometry.read_text(encoding="utf-8"))
    if geometry.get("schema") == PARTITION_LEDGER_SCHEMA and geometry.get("certified") is not True:
        result = fail_closed_result(
            geometry,
            ns.geometry,
            ns.trial,
            "partition ledger is explicitly uncertified and supplies no complete overlap boxes/permutations",
        )
        ns.json_out.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
        print(json.dumps({
            "certified": False,
            "weighted_residual_upper": None,
            "blockers": result["blockers"],
        }, indent=2))
        print(ns.json_out.resolve())
        return 0
    validate_geometry(geometry)
    trial, row = parse_trial(ns.trial)
    ev = DifferentialTrial(trial["parameters"]["M"], str(row["r"]), row["coefficients"])
    levels = []
    for n in [int(q) for q in ns.refinements.split(",") if q.strip()]:
        print(f"overlap refinement={n}", flush=True)
        levels.append(certify_level(ev, geometry, n))
    monotone = all(
        levels[k + 1]["delta0_upper"] <= levels[k]["delta0_upper"]
        and levels[k + 1]["delta1_upper"] <= levels[k]["delta1_upper"]
        for k in range(len(levels) - 1)
    )
    final = levels[-1]
    b0 = geometry.get("b0_upper")
    b1 = geometry.get("b1_upper")
    weighted_ball = None if b0 is None or b1 is None else (
        arb(str(b0)) * arb(repr(final["delta0_upper"]))
        + arb(str(b1)) * arb(repr(final["delta1_upper"]))
    )
    weighted = None if weighted_ball is None else upper(weighted_ball)
    threshold_ball = arb("0.0101903405004245")
    threshold = upper(threshold_ball)
    weighted_below_threshold = bool(
        weighted_ball is not None and weighted_ball.upper() < threshold_ball.lower()
    )
    certified = bool(
        geometry.get("partition_constants_certified") is True
        and monotone and weighted_below_threshold
    )
    result = {
        "status": "cancellation-preserving Arb Track-B overlap certificate",
        "certified": certified,
        "delta0_upper": final["delta0_upper"],
        "delta1_upper": final["delta1_upper"],
        "direct_weighted_defect_L2_upper": weighted,
        "weighted_residual_upper": weighted,
        "track_B_R_target_upper": threshold,
        "closes_lambda_width_below_0.1": bool(certified and weighted_below_threshold),
        "active_transition_ids": [tr["id"] for tr in geometry["active_transitions"]],
        "all_transitions_covered": geometry["transition_set_complete"],
        "common_fiber_transport_certified": True,
        "stabilizer_averaging_certified": geometry["stabilizers_complete"],
        "first_gradients_certified": True,
        "theorem_DK_compatible": certified,
        "transport_convention": TRANSPORT,
        "stabilization_certified": geometry["stabilizers_complete"],
        "two_cusp_coordinates_certified": geometry["two_cusp_coordinates_certified"],
        "cusp_coordinate_tags": sorted({tr["cusp_coordinate"] for tr in geometry["active_transitions"]}),
        "refinement_table": levels,
        "refinement_monotone": monotone,
        "independent_check": midpoint_jacobian_check(geometry),
        "theorem_compatibility": {
            "theorem": "Track-B six-copy/two-cusp D(K)",
            "tau": 0,
            "six_vector_norm": True,
            "first_covariant_gradient": True,
            "cancellation_preserved_before_norm": True,
            "partition_constants_certified": geometry.get("partition_constants_certified") is True,
            "admissible": certified,
        },
        "geometry": str(ns.geometry.resolve()),
        "geometry_sha256": hashlib.sha256(ns.geometry.read_bytes()).hexdigest(),
        "trial": str(ns.trial.resolve()),
        "trial_sha256": hashlib.sha256(ns.trial.read_bytes()).hexdigest(),
        "precision_bits": int(ctx.prec),
        "bessel_enclosures": {
            "value_direct_or_mean_value": True,
            "derivatives_direct_arb_shifted_order": True,
            "absolute_modal_taylor_majorant_used": False,
        },
        "rung4_certified": False,
    }
    ns.json_out.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps({k: result[k] for k in (
        "certified", "delta0_upper", "delta1_upper",
        "direct_weighted_defect_L2_upper", "closes_lambda_width_below_0.1"
    )}, indent=2))
    print(ns.json_out.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
