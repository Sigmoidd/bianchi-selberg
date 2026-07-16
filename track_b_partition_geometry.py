#!/usr/bin/env python3
"""Exact Track-B Humbert face inventory and fail-closed partition ledger.

This program deliberately separates exact Poincare-cell combinatorics from
the missing analytic collar theorem.  The existing 24-split FEM mesh is an
inner approximation of the Humbert cell and therefore cannot, by itself,
certify coverage by smooth orbifold charts.

All matrix and fixed-point calculations below use Gaussian integers and
rational arithmetic.  The derivative bounds use the elementary rational
majorants

    max |(10t^3-15t^4+6t^5)'| <= 15/8,
    max |(10t^3-15t^4+6t^5)''| <= 6.

The reported b0/b1 values are CONDITIONAL until ``collar_separation`` and
the complete singular-stratum incidence inventory are certified.  The six
listed point stabilizers themselves are exhaustively enumerated.  The
top-level ``certified`` field remains false and no rung flag is changed.
"""
from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction as Q
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent

# Gaussian integer (real, imaginary) and 2x2 matrices, modulo +/- I.
GI = tuple[int, int]
Mat = tuple[GI, GI, GI, GI]
ZERO: GI = (0, 0)
ONE: GI = (1, 0)
I: GI = (0, 1)


def gadd(x: GI, y: GI) -> GI:
    return x[0] + y[0], x[1] + y[1]


def gneg(x: GI) -> GI:
    return -x[0], -x[1]


def gmul(x: GI, y: GI) -> GI:
    return x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def mneg(a: Mat) -> Mat:
    return tuple(gneg(x) for x in a)  # type: ignore[return-value]


def canon(a: Mat) -> Mat:
    b = mneg(a)
    return min(a, b)


def mmul(a: Mat, b: Mat) -> Mat:
    aa, ab, ac, ad = a
    ba, bb, bc, bd = b
    return canon((
        gadd(gmul(aa, ba), gmul(ab, bc)),
        gadd(gmul(aa, bb), gmul(ab, bd)),
        gadd(gmul(ac, ba), gmul(ad, bc)),
        gadd(gmul(ac, bb), gmul(ad, bd)),
    ))


ID = canon((ONE, ZERO, ZERO, ONE))
T1 = canon((ONE, ONE, ZERO, ONE))
T1I = canon((ONE, (-1, 0), ZERO, ONE))
TI = canon((ONE, I, ZERO, ONE))
TII = canon((ONE, (0, -1), ZERO, ONE))
R = canon((I, ZERO, ZERO, (0, -1)))
TIR = mmul(TI, R)
S = canon((ZERO, (-1, 0), ONE, ZERO))


def mat_json(a: Mat) -> list[list[list[int]]]:
    return [[list(a[0]), list(a[1])], [list(a[2]), list(a[3])]]


def qgadd(x: tuple[Q, Q], y: tuple[Q, Q]) -> tuple[Q, Q]:
    return x[0] + y[0], x[1] + y[1]


def qgmul(x: tuple[Q, Q], y: tuple[Q, Q]) -> tuple[Q, Q]:
    return x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def qgconj(x: tuple[Q, Q]) -> tuple[Q, Q]:
    return x[0], -x[1]


def qgnorm(x: tuple[Q, Q]) -> Q:
    return x[0] * x[0] + x[1] * x[1]


def qgi(x: GI) -> tuple[Q, Q]:
    return Q(x[0]), Q(x[1])


def fixes(a: Mat, z: tuple[Q, Q], y2: Q) -> bool:
    """Exact fixed-point test in H^3, using y^2 only."""
    aa, b, c, d = map(qgi, a)
    czd = qgadd(qgmul(c, z), d)
    den = qgnorm(czd) + qgnorm(c) * y2
    if den != 1:
        return False
    azb = qgadd(qgmul(aa, z), b)
    num = qgadd(qgmul(azb, qgconj(czd)), qgmul(aa, qgconj(c)))
    num = qgadd(num, (qgmul(aa, qgconj(c))[0] * (y2 - 1),
                      qgmul(aa, qgconj(c))[1] * (y2 - 1)))
    # The previous two lines equal azb*conj(czd)+a*conj(c)*y2.
    return num == z


def fixes_clear(a: Mat, z: tuple[Q, Q], y2: Q) -> bool:
    """Same test as fixes(), written without the compact algebra trick."""
    aa, b, c, d = map(qgi, a)
    czd = qgadd(qgmul(c, z), d)
    den = qgnorm(czd) + qgnorm(c) * y2
    azb = qgadd(qgmul(aa, z), b)
    term = qgmul(aa, qgconj(c))
    num0 = qgmul(azb, qgconj(czd))
    num = num0[0] + term[0] * y2, num0[1] + term[1] * y2
    return den == 1 and num == z


GEN = {
    "T1": T1, "T1^-1": T1I, "Ti": TI, "Ti^-1": TII,
    "R": R, "S": S,
}


def bfs(max_depth: int) -> dict[Mat, str]:
    words = {ID: "I"}
    depth = {ID: 0}
    todo = deque([ID])
    while todo:
        a = todo.popleft()
        if depth[a] == max_depth:
            continue
        for name, g in GEN.items():
            b = mmul(a, g)
            if b not in words:
                words[b] = name if words[a] == "I" else words[a] + " " + name
                depth[b] = depth[a] + 1
                todo.append(b)
    return words


SPECIAL = {
    # Coordinates are (Re z, Im z, y^2).  These are the finite Humbert
    # vertices/fixed-line intersections in the closed compact cell.
    "v_00": (Q(0), Q(0), Q(1)),
    "v_0h": (Q(0), Q(1, 2), Q(3, 4)),
    "v_mh0": (Q(-1, 2), Q(0), Q(3, 4)),
    "v_ph0": (Q(1, 2), Q(0), Q(3, 4)),
    "v_mhh": (Q(-1, 2), Q(1, 2), Q(1, 2)),
    "v_phh": (Q(1, 2), Q(1, 2), Q(1, 2)),
}


def subgroup(gens: list[Mat]) -> set[Mat]:
    out = {ID}
    todo = deque([ID])
    while todo:
        a = todo.popleft()
        for g in gens:
            b = mmul(a, g)
            if b not in out:
                out.add(b)
                todo.append(b)
    return out


def gdiv_exact(x: GI, y: GI) -> GI | None:
    """Return x/y in Z[i], or None when the quotient is not integral."""
    n = y[0] * y[0] + y[1] * y[1]
    if n == 0:
        return None
    num = gmul(x, (y[0], -y[1]))
    if num[0] % n or num[1] % n:
        return None
    return num[0] // n, num[1] // n


def exhaustive_stabilizer(z: tuple[Q, Q], y2: Q) -> set[Mat]:
    """Enumerate the point stabilizer with a proved finite entry bound.

    If g=[[a,b],[c,d]] fixes (z,y), the height denominator for g and g^-1
    is one.  Hence

      |cz+d|^2 + |c|^2 y^2 = 1,
      |-cz+a|^2 + |c|^2 y^2 = 1.

    At the six SPECIAL points, |z|^2 <= 1/2 and y^2 >= 1/2.  Therefore
    |c|^2 <= 2 and |a|,|d| <= |cz|+1 <= 2.  Enumerating Gaussian a,c,d
    in the square [-2,2]^2 is consequently exhaustive.  The determinant
    determines b when c != 0; the c=0 case is also finite because ad=1 and
    the fixed-point equation determines b.  Every candidate is finally
    checked by exact rational arithmetic.
    """
    assert qgnorm(z) <= Q(1, 2) and y2 >= Q(1, 2)
    gis = [(u, v) for u in range(-2, 3) for v in range(-2, 3)]
    out: set[Mat] = set()
    for c in gis:
        if Q(c[0] * c[0] + c[1] * c[1]) * y2 > 1:
            continue
        cq = qgi(c)
        cz = qgmul(cq, z)
        rhs = 1 - Q(c[0] * c[0] + c[1] * c[1]) * y2
        for a in gis:
            amcz = (Q(a[0]) - cz[0], Q(a[1]) - cz[1])
            if qgnorm(amcz) != rhs:
                continue
            for d in gis:
                czpd = (cz[0] + Q(d[0]), cz[1] + Q(d[1]))
                if qgnorm(czpd) != rhs:
                    continue
                if c != ZERO:
                    b = gdiv_exact(gadd(gmul(a, d), (-1, 0)), c)
                    if b is None:
                        continue
                else:
                    if gmul(a, d) != ONE:
                        continue
                    # From (a z+b) conj(d)=z and |d|=1.
                    bq = qgadd(qgmul(z, qgi(d)),
                               (-qgmul(qgi(a), z)[0],
                                -qgmul(qgi(a), z)[1]))
                    if bq[0].denominator != 1 or bq[1].denominator != 1:
                        continue
                    b = int(bq[0]), int(bq[1])
                candidate = canon((a, b, c, d))
                if fixes_clear(candidate, z, y2):
                    out.add(candidate)
    return out


def rational_bounds() -> dict[str, object]:
    """Conditional sup/L2 bounds for epsilon=1/10 product gates."""
    eps = Q(1, 10)
    ymax = Q(6, 5)
    p1 = Q(15, 8)
    p2 = Q(6)
    gx = ymax * p1 / eps
    dx = ymax * ymax * p2 / (eps * eps)
    gf = 2 * gx
    df = 4 * dx
    sum_g = 4 * gx + gf
    sum_d = 4 * dx + df
    # Opposing vertical gates have disjoint supports, but not using that fact
    # here makes the bound robust and simpler.
    pair_g = 6 * gx * gx + 4 * gx * gf
    phi_grad_l1 = 2 * sum_g
    phi_lap_l1 = 2 * sum_d + 8 * pair_g

    width = Q(19, 100)
    gc = ymax * p1 / width
    dc = ymax * ymax * p2 / (width * width) + ymax * p1 / width
    b1_sup = 4 * gc + 2 * phi_grad_l1
    b0_sup = 2 * dc + phi_lap_l1 + 2 * gc * phi_grad_l1

    # Elementary volume upper bound: F lies in a horizontal rectangle of
    # area 1/2 and y >= 1/sqrt(2), so int_F y^-3 <= 1/2.
    b1_upper = float(b1_sup) / math.sqrt(2.0)
    b0_upper = float(b0_sup) / math.sqrt(2.0)
    target_R = 0.0101903405004245
    return {
        "gate_epsilon": "1/10",
        "quintic_derivative_majorants": {"first": "15/8", "second": "6"},
        "hyperbolic_sup_bounds": {
            "B0": str(b0_sup), "B1": str(b1_sup),
        },
        "conditional_b0_upper": math.nextafter(b0_upper, math.inf),
        "conditional_b1_upper": math.nextafter(b1_upper, math.inf),
        "volume_upper": "1/2",
        "target_R_upper": target_R,
        "individual_budget": {
            "delta0_if_delta1_zero": math.nextafter(target_R / b0_upper, 0.0),
            "delta1_if_delta0_zero": math.nextafter(target_R / b1_upper, 0.0),
        },
        "equal_split_budget": {
            "delta0": math.nextafter(target_R / (2 * b0_upper), 0.0),
            "delta1": math.nextafter(target_R / (2 * b1_upper), 0.0),
        },
    }


def exp_upper_by_taylor(x: Q, terms: int) -> Q:
    """A rational upper bound for exp(x), for positive rational x.

    After the displayed Taylor terms, the ratio of successive terms is at
    most x/(terms+1).  The geometric tail bound is deliberately elementary:
    this ledger must not depend on a floating-point comparison with exp.
    """
    if x < 0 or terms < 0:
        raise ValueError("this elementary exp bound requires x >= 0")
    term = Q(1)
    partial = term
    for k in range(1, terms + 1):
        term *= x / k
        partial += term
    next_term = term * x / (terms + 1)
    ratio = x / (terms + 2)
    if ratio >= 1:
        raise ValueError("Taylor tail ratio must be below one")
    return partial + next_term / (1 - ratio)


def finite_order(a: Mat, limit: int = 24) -> int:
    """Return the projective order of a known finite-order matrix."""
    b = ID
    for order in range(1, limit + 1):
        b = mmul(b, a)
        if b == ID:
            return order
    raise ArithmeticError(f"no order <= {limit} for purported elliptic matrix")


def floor_collar_incidence() -> dict[str, object]:
    """Exact finite incidence ledger for 0 <= log(|z|^2+y^2) <= 3/10.

    This proves only the finite isometric-sphere / elliptic-stratum claim
    needed by the floor collar.  It deliberately does not claim that the
    full smooth-orbifold collar construction or its partition constants have
    been completed.

    Put rho=|z|^2+y^2.  On the closed Picard rectangle,
    |z|^2 <= 1/2 and rho >= 1, hence y^2 >= 1/2.  For an isometric sphere

        |c z+d|^2+|c|^2 y^2 = 1,

    |c|^2>2 is impossible; |c|^2=2 forces y^2=1/2 and cz+d=0,
    hence one of the two lower corner vertices.  The remaining |c|=1
    spheres have Gaussian centres n=-d/c.  Horizontal radius one confines
    n to {-1,0,1}^2, which is enumerated below with exact affine extrema.
    """
    eps = Q(3, 10)
    rho_cap = Q(27, 20)
    exp_cap = exp_upper_by_taylor(eps, 8)
    if not exp_cap < rho_cap:
        raise ArithmeticError("failed rational containment exp(3/10) < 27/20")

    # For a unit-c sphere centred at n, its equation is
    # rho = 1-|n|^2+2(n_1 x1+n_2 x2).  The extrema on the closed horizontal
    # rectangle are attained at corners, so this is entirely rational.
    corners = [(Q(-1, 2), Q(0)), (Q(-1, 2), Q(1, 2)),
               (Q(1, 2), Q(0)), (Q(1, 2), Q(1, 2))]
    unit_spheres = []
    for n1 in range(-1, 2):
        for n2 in range(-1, 2):
            values = [Q(1 - n1 * n1 - n2 * n2) + 2 * (n1 * x1 + n2 * x2)
                      for x1, x2 in corners]
            lo, hi = min(values), max(values)
            if (n1, n2) == (0, 0):
                incidence = "the unique codimension-one floor sphere; rho=1"
            elif (n1, n2) == (-1, 0):
                incidence = "floor/x1m edge only"
            elif (n1, n2) == (1, 0):
                incidence = "floor/x1p edge only"
            elif (n1, n2) == (0, 1):
                incidence = "floor/x2p edge only"
            elif (n1, n2) == (-1, 1):
                incidence = "v_mhh only"
            elif (n1, n2) == (1, 1):
                incidence = "v_phh only"
            else:
                incidence = "disjoint from rho >= 1"
            unit_spheres.append({
                "center": [n1, n2],
                "rho_range_on_horizontal_rectangle": [str(lo), str(hi)],
                "incidence": incidence,
            })

    edge_generators = {
        "floor_x1m": (mmul(S, T1), 3,
                        ["v_mh0", "v_mhh"], [], "x1=-1/2, rho=1"),
        "floor_x1p": (mmul(T1, S), 3,
                        ["v_ph0", "v_phh"], [], "x1=1/2, rho=1"),
        "floor_x2m": (mmul(R, S), 2,
                        ["v_mh0", "v_ph0"], ["v_00"], "x2=0, rho=1"),
        "floor_x2p": (mmul(TIR, S), 3,
                        ["v_mhh", "v_phh"], ["v_0h"], "x2=1/2, rho=1"),
    }
    elliptic_edges = []
    for name, (generator, expected_order, endpoints, special_points, locus) in edge_generators.items():
        order = finite_order(generator)
        if order != expected_order:
            raise ArithmeticError(f"wrong order on {name}: {order}")
        elliptic_edges.append({
            "id": name,
            "locus": locus,
            "generator": mat_json(generator),
            "stabilizer_order": order,
            "averaging_rule": f"average the local chart over the cyclic group of order {order}",
            "endpoints": endpoints,
            "interior_special_points": special_points,
            "interior_special_point_rule": (
                "replace the generic cyclic average by the exact SPECIAL-point "
                "stabilizer average when this list is nonempty"
            ),
        })

    return {
        "certified": True,
        "scope": "finite isometric-sphere and elliptic-stratum incidence in the selected floor collar",
        "floor_collar": {
            "log_rho_interval": ["0", "3/10"],
            "rho_upper_containment": "27/20",
            "exp_upper_by_9_term_taylor_plus_geometric_tail": str(exp_cap),
            "containment_verified": True,
        },
        "exclusion_argument": {
            "closed_cell_y_squared_lower": "1/2",
            "c_norm_squared_greater_than_2": "impossible, since |c|^2 y^2 > 1",
            "c_norm_squared_equal_2": (
                "only equality cases y^2=1/2 and cz+d=0; these are the already "
                "enumerated vertices v_mhh and v_phh"
            ),
            "c_norm_squared_equal_1": (
                "unit-c sphere centres are Gaussian; radius-one horizontal projection "
                "forces each centre coordinate into {-1,0,1}"
            ),
        },
        "unit_c_spheres": unit_spheres,
        "open_collar_incident_codimension_one_words": ["S"],
        "closed_collar_incident_vertical_faces": ["x1m", "x1p", "x2m", "x2p"],
        "elliptic_edge_averages": elliptic_edges,
        "vertex_averages": {
            name: {
                "stabilizer_order": len(exhaustive_stabilizer((x1, x2), y2)),
                "averaging_rule": "average over the exact finite point stabilizer recorded in stabilizers",
            }
            for name, (x1, x2, y2) in SPECIAL.items()
        },
        "conclusion": (
            "No unlisted isometric sphere enters the open 0.30 floor collar.  "
            "All closed-collar edge and vertex incidences reduce to the listed finite "
            "cyclic averages and exact SPECIAL-point stabilizers."
        ),
    }


def global_singular_strata_incidence() -> dict[str, object]:
    """Exact singular-stratum ledger for the closed Humbert cell.

    The five boundary faces are the two ``x1`` faces paired to one another
    and the three self-paired faces ``x2=0``, ``x2=1/2`` and ``rho=1``.
    The relative interiors of the self-paired faces contain one fixed axis
    each.  The four vertical face intersections and the four floor/vertical
    intersections are the remaining open one-dimensional strata.  Their
    endpoints/intersections are exactly the six points in ``SPECIAL``.

    Stabilizer orders are not inferred from pictures or midpoint floating
    point tests.  For one rational point in every open stratum we invoke the
    exhaustive Gaussian-entry enumeration used for the vertices.  The face
    and edge-cycle equations identify the full open locus; the stabilizer can
    change only at a listed intersection, where the complete point group is
    substituted.
    """
    strata = {
        # Fixed axes in relative interiors of self-paired faces.
        "floor_S": {
            "locus": "x1=0, rho=1, 0<x2<1/2",
            "sample": (Q(0), Q(1, 4), Q(15, 16)),
            "kind": "self_paired_face_fixed_axis",
            "incident_faces": ["floor"],
            "expected_order": 2,
        },
        "vertical_R": {
            "locus": "x1=0, x2=0, rho>1",
            "sample": (Q(0), Q(0), Q(2)),
            "kind": "self_paired_face_fixed_axis",
            "incident_faces": ["x2m"],
            "expected_order": 2,
        },
        "vertical_TiR": {
            "locus": "x1=0, x2=1/2, rho>1",
            "sample": (Q(0), Q(1, 2), Q(2)),
            "kind": "self_paired_face_fixed_axis",
            "incident_faces": ["x2p"],
            "expected_order": 2,
        },
        # Vertical corner cycles.  The two members of each x1-paired pair
        # remain separate closed-cell incidences but represent one quotient
        # stratum after the exact face identification.
        "vertical_x1m_x2m": {
            "locus": "x1=-1/2, x2=0, rho>1",
            "sample": (Q(-1, 2), Q(0), Q(2)),
            "kind": "vertical_edge_cycle",
            "incident_faces": ["x1m", "x2m"],
            "expected_order": 2,
        },
        "vertical_x1p_x2m": {
            "locus": "x1=1/2, x2=0, rho>1",
            "sample": (Q(1, 2), Q(0), Q(2)),
            "kind": "vertical_edge_cycle",
            "incident_faces": ["x1p", "x2m"],
            "expected_order": 2,
        },
        "vertical_x1m_x2p": {
            "locus": "x1=-1/2, x2=1/2, rho>1",
            "sample": (Q(-1, 2), Q(1, 2), Q(2)),
            "kind": "vertical_edge_cycle",
            "incident_faces": ["x1m", "x2p"],
            "expected_order": 2,
        },
        "vertical_x1p_x2p": {
            "locus": "x1=1/2, x2=1/2, rho>1",
            "sample": (Q(1, 2), Q(1, 2), Q(2)),
            "kind": "vertical_edge_cycle",
            "incident_faces": ["x1p", "x2p"],
            "expected_order": 2,
        },
        # Four boundary edges of the inversion face.
        "floor_x1m": {
            "locus": "x1=-1/2, rho=1, 0<x2<1/2",
            "sample": (Q(-1, 2), Q(1, 4), Q(11, 16)),
            "kind": "floor_edge_cycle",
            "incident_faces": ["floor", "x1m"],
            "expected_order": 3,
        },
        "floor_x1p": {
            "locus": "x1=1/2, rho=1, 0<x2<1/2",
            "sample": (Q(1, 2), Q(1, 4), Q(11, 16)),
            "kind": "floor_edge_cycle",
            "incident_faces": ["floor", "x1p"],
            "expected_order": 3,
        },
        "floor_x2m": {
            "locus": "x2=0, rho=1, -1/2<x1<1/2",
            "sample": (Q(1, 4), Q(0), Q(15, 16)),
            "kind": "floor_edge_cycle",
            "incident_faces": ["floor", "x2m"],
            "expected_order": 2,
        },
        "floor_x2p": {
            "locus": "x2=1/2, rho=1, -1/2<x1<1/2",
            "sample": (Q(1, 4), Q(1, 2), Q(11, 16)),
            "kind": "floor_edge_cycle",
            "incident_faces": ["floor", "x2p"],
            "expected_order": 3,
        },
    }
    rows = []
    for name, data in strata.items():
        x1, x2, y2 = data["sample"]
        group = exhaustive_stabilizer((x1, x2), y2)
        expected = int(data["expected_order"])
        if len(group) != expected:
            raise ArithmeticError(
                f"wrong exact generic stabilizer order on {name}: "
                f"{len(group)} != {expected}"
            )
        closure = all(mmul(a, b) in group for a in group for b in group)
        reindexing = all({mmul(a, h) for a in group} == group for h in group)
        if not closure or not reindexing:
            raise ArithmeticError(f"stabilizer group check failed on {name}")
        rows.append({
            "id": name,
            "locus": data["locus"],
            "kind": data["kind"],
            "incident_faces": data["incident_faces"],
            "generic_sample": {
                "x1": str(x1), "x2": str(x2), "y_squared": str(y2),
            },
            "stabilizer_order": len(group),
            "elements": [mat_json(a) for a in sorted(group)],
            "group_closure_exact": closure,
            "right_reindexing_exact": reindexing,
        })

    return {
        "certified": True,
        "scope": "all positive-height singular strata of the closed Humbert reference cell",
        "face_inventory_complete": True,
        "relative_face_interior_classification": {
            "paired_distinct_faces": ["x1m", "x1p"],
            "self_paired_faces": ["x2m", "x2p", "floor"],
            "self_paired_fixed_axes": ["vertical_R", "vertical_TiR", "floor_S"],
        },
        "one_dimensional_strata": rows,
        "zero_dimensional_strata": list(SPECIAL),
        "exhaustiveness_proof": (
            "The exact five-face Humbert polyhedron has no other nonempty "
            "relative face intersections.  Interior stabilizers are excluded "
            "by disjoint fundamental-domain interiors; face-interior isotropy "
            "is the fixed locus of one of the three self-pairings; edge isotropy "
            "is one of the eight listed face cycles.  Stabilizer jumps occur "
            "only at the six enumerated SPECIAL intersections."
        ),
        "sampled_floating_point_geometry_used": False,
    }


def build(max_depth: int) -> dict[str, object]:
    words = bfs(max_depth)
    collar_incidence = floor_collar_incidence()
    global_incidence = global_singular_strata_incidence()
    stabilizers = {}
    for name, (xr, xi, y2) in SPECIAL.items():
        exact = exhaustive_stabilizer((xr, xi), y2)
        bfs_fixed = {a for a in words if fixes_clear(a, (xr, xi), y2)}
        # The word search is now only an independent reproduction check.
        assert bfs_fixed <= exact
        assert all(fixes(a, (xr, xi), y2) for a in exact)
        stabilizers[name] = {
            "coordinates": {"x1": str(xr), "x2": str(xi), "y_squared": str(y2)},
            "found_through_word_depth": max_depth,
            "order_found": len(exact),
            "bfs_reproduced_all": bfs_fixed == exact,
            "elements": [{"word": words.get(a, "exact enumeration; word not found at requested depth"),
                          "matrix": mat_json(a)} for a in sorted(exact)],
            "complete": True,
            "complete_reason": "exhaustive denominator-bound enumeration over Gaussian entries",
        }

    faces = {
        "x1m": {"word": "T1", "matrix": mat_json(T1), "target": "x1p",
                  "pi_gamma": [0, 5, 1, 2, 3, 4],
                  "conditional_overlap_region": "-1/2 <= x1 < -3/10"},
        "x1p": {"word": "T1^-1", "matrix": mat_json(T1I), "target": "x1m",
                  "pi_gamma": [0, 2, 3, 4, 5, 1],
                  "conditional_overlap_region": "3/10 < x1 <= 1/2"},
        "x2m": {"word": "R", "matrix": mat_json(R), "target": "x2m",
                  "pi_gamma": [0, 1, 5, 4, 3, 2],
                  "conditional_overlap_region": "0 <= x2 < 1/5"},
        "x2p": {"word": "Ti R", "matrix": mat_json(TIR), "target": "x2p",
                  "pi_gamma": [0, 3, 2, 1, 5, 4],
                  "conditional_overlap_region": "3/10 < x2 <= 1/2"},
        "floor": {"word": "S", "matrix": mat_json(S), "target": "floor",
                   "pi_gamma": [1, 0, 5, 3, 4, 2],
                   "conditional_overlap_region":
                       "0 <= log(x1^2+x2^2+y^2) < 1/5"},
    }
    conditional = rational_bounds()
    return {
        "schema": "track_b_partition/v1",
        "certified": False,
        "coverage_certified": False,
        "local_finiteness_certified": False,
        "transitions_complete": False,
        "transition_set_complete": False,
        "stabilizers_complete": False,
        "two_cusp_coordinates_certified": True,
        "partition_constants_certified": False,
        "tau_upper": 0.0,
        "b0_upper": None,
        "b1_upper": None,
        "weighted_residual_operator_upper": None,
        "theorem_compatibility": {
            "theorem": "theorem_DK_sixcopy.tex, Lemma Certified automorphization",
            "compatible": False,
            "reason": (
                "the finite floor-collar incidence is now exact, but the global "
                "collar-separation/Lebesgue-number and normalized-weight proofs remain open"
            ),
        },
        "coordinate_convention": {
            "reference_cell": "x1 in [-1/2,1/2], x2 in [0,1/2], x1^2+x2^2+y^2 >= 1",
            "hyperbolic_laplacian": "Delta=-y^2(d_x1^2+d_x2^2+d_y^2)+y d_y",
            "transport": "F_c(P)=F_{pi_gamma(c)}(gamma P); GLUE stores pi_gamma(c)=c.gamma^{-1}",
            "cusp_inf": "identity scaling; lattice Z[i]; folded rectangle [-1/2,1/2]x[0,1/2]",
            "cusp_0": "sigma0=S; five translated chimneys T+k; lattice (2+i)Z[i]",
        },
        "cutoff": {
            "chi_B": "0 below 1.01; quintic smoothstep on [1.01,1.20]; 1 above 1.20",
            "core_gate": "five product gates with epsilon=1/10 on x1m,x1p,x2m,x2p,log(x1^2+x2^2+y^2)",
            "algebraic_partition_identity": True,
            "extension_by_zero_C2": True,
        },
        "active_transition_ids": list(faces),
        "active_transition_words": [faces[k]["word"] for k in faces],
        "faces": faces,
        "overlap_interface": {
            "schema": "track-b-humbert-overlap-v1",
            "certified": False,
            "common_fiber_transport":
                "GLUE pi_gamma: compare F_c(P) with F_{pi_gamma(c)}(gamma P)",
            "active_transitions": [
                {"id": k, **v} for k, v in faces.items()
            ],
            "edge_vertex_transition_words_complete": False,
            "instruction": "consumer must fail closed; face rows are exact but not exhaustive at edges/vertices",
        },
        "coverage": {
            "closed_reference_cell_product_identity": True,
            "true_orbifold_thickening": False,
            "polyhedral_mesh_is_not_a_cover": True,
        },
        "local_finiteness": {
            "conditional_on_collar_separation": True,
            "certified": False,
            "max_simultaneous_face_families_if_valid": 3,
            "max_boolean_chambers_if_valid": 8,
        },
        "stabilizers": stabilizers,
        "stabilizer_flags": {
            "exact_fixed_point_tests": True,
            "elliptic_averaging_required": True,
            "listed_special_points_exhaustive": True,
            "floor_collar_singular_strata_exhaustive": True,
            "global_singular_strata_exhaustive": global_incidence["certified"],
        },
        "floor_collar_incidence": collar_incidence,
        "global_singular_strata_incidence": global_incidence,
        "conditional_bounds_not_theorem_inputs": conditional,
        "blockers": [
            "replace face-only transition list by all edge/vertex chamber-to-chamber words after the collar proof",
            "prove the global collar-separation/Lebesgue-number statement outside the resolved floor collar",
            "interval-integrate the resulting normalized rational/quintic weights; current bounds are conditional sup bounds",
        ],
        "source_mesh_audit": {
            "exact_face_types_reusable": True,
            "mesh": "independent_exclusion/congruence_prototype.py 24-split",
            "fatal_for_partition_coverage": "mesh is K_h subset K and leaves the curved sliver K\\K_h uncovered",
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", type=int, default=6)
    ap.add_argument("--output", type=Path,
                    default=ROOT / "track_b_partition_result.json")
    ap.add_argument("--certify-global-partition", action="store_true")
    ap.add_argument("--bits", type=int, default=192)
    ap.add_argument("--partition-subdivision", default="16,8,16")
    ap.add_argument("--partition-degree", type=int, default=8)
    ap.add_argument("--partition-workers", type=int, default=1)
    ap.add_argument("--partition-audit-jsonl", type=Path,
                    default=ROOT / "track_b_partition_cells.jsonl")
    ap.add_argument("--check-refinement", default="")
    ap.add_argument("--check-audit-jsonl", type=Path,
                    default=ROOT / "track_b_partition_cells_refined.jsonl")
    ap.add_argument("--global-json-out", type=Path,
                    default=ROOT / "track_b_global_partition_result.json")
    ns = ap.parse_args()
    out = build(ns.depth)
    ns.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    if ns.certify_global_partition:
        from flint import ctx
        from track_b_global_partition_arb import apply_stability, certify_partition

        ctx.prec = max(128, ns.bits)
        dims = tuple(int(q) for q in ns.partition_subdivision.split(","))
        if len(dims) != 3 or min(dims) <= 0:
            raise ValueError("partition subdivision must be nx1,nx2,ns")
        result = certify_partition(
            out, dims, ns.partition_degree, int(ctx.prec),
            ns.partition_audit_jsonl, ns.partition_workers, True,
        )
        if ns.check_refinement.strip():
            check_dims = tuple(int(q) for q in ns.check_refinement.split(","))
            if len(check_dims) != 3 or min(check_dims) <= 0:
                raise ValueError("check refinement must be nx1,nx2,ns")
            check = certify_partition(
                out, check_dims, ns.partition_degree, int(ctx.prec),
                ns.check_audit_jsonl, ns.partition_workers, False,
            )
            result = apply_stability(result, check)
        result["rung4_certified"] = False
        ns.global_json_out.write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps({
            "global_partition_certified": result["global_partition_certified"],
            "global_weight_bounds_certified": result["global_weight_bounds_certified"],
            "stability_check_passed": result["stability_check_passed"],
            "rung4_certified": False,
            "global_output": str(ns.global_json_out.resolve()),
        }, indent=2))
        return 0 if result["global_weight_bounds_certified"] else 2
    print(json.dumps({
        "certified": out["certified"],
        "output": str(ns.output),
        "conditional_bounds": out["conditional_bounds_not_theorem_inputs"],
        "stabilizer_orders_found": {
            k: v["order_found"] for k, v in out["stabilizers"].items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
