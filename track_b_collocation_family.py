#!/usr/bin/env python3
"""Deterministic exact non-symmetric collocation families for Track B."""
from __future__ import annotations

import argparse
import cmath
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

from track_b_two_cusp_data import GENERATORS, GLUE, canonical_hash, cusp_modes


ROOT = Path(__file__).resolve().parent

# Hand-auditable rational points.  Heights exceed one, so every point lies in
# the embedded cusp part of the closed Humbert core; horizontal coordinates
# avoid x1=0, x2=0 and x1=+/-x2.  Denominators and height layers deliberately
# vary instead of forming a symmetry-closed grid.
_EXACT_POINTS = (
    ("1/5", "2/7", "17/15"),
    ("-2/7", "1/5", "19/16"),
    ("1/8", "-3/11", "21/19"),
    ("-1/9", "-2/5", "23/20"),
    ("2/11", "3/8", "25/21"),
    ("-3/13", "2/9", "27/23"),
    ("4/17", "-1/7", "29/24"),
    ("-2/15", "-3/17", "31/25"),
    ("3/19", "4/21", "33/26"),
    ("-4/21", "5/23", "35/27"),
    ("5/22", "-3/20", "37/29"),
    ("-5/24", "4/19", "39/31"),
    ("2/13", "-5/18", "41/32"),
    ("-3/20", "7/25", "43/33"),
    ("5/27", "2/17", "45/34"),
    ("-7/29", "-3/22", "47/35"),
    ("4/23", "-5/31", "49/36"),
    ("-5/26", "7/32", "51/37"),
    ("7/30", "-4/27", "53/38"),
    ("-6/31", "5/28", "55/39"),
)


def f(text: str) -> Fraction:
    return Fraction(text)


def exact_candidate_points(count: int = 20) -> list[dict[str, Any]]:
    if count < 1 or count > len(_EXACT_POINTS):
        raise ValueError(f"candidate point count must be in [1,{len(_EXACT_POINTS)}]")
    points = []
    for index, (x1s, x2s, ys) in enumerate(_EXACT_POINTS[:count]):
        x1, x2, y = f(x1s), f(x2s), f(ys)
        if not y > 1:
            raise ArithmeticError("candidate height must exceed one")
        if x1 == 0 or x2 == 0 or x1 == x2 or x1 == -x2:
            raise ArithmeticError("candidate lies on a forbidden symmetry axis")
        if x1 * x1 + x2 * x2 + y * y <= 1:
            raise ArithmeticError("candidate lies outside the closed Humbert core")
        points.append({
            "point_id": f"q{index:02d}",
            "x1": x1s,
            "x2": x2s,
            "y": ys,
            "height_layer": ys,
            "horizontal_denominators": [x1.denominator, x2.denominator],
            "design": "exact rational non-symmetric Humbert-core point",
        })
    return points


def candidate_point_checks(points: list[dict[str, Any]]) -> dict[str, bool]:
    coordinates = [(str(q["x1"]), str(q["x2"]), str(q["y"])) for q in points]
    ids = [str(q.get("point_id")) for q in points]
    physical = True
    axes_avoided = True
    generator_fixed_sets_avoided = True
    for record in points:
        try:
            x1, x2, y = f(record["x1"]), f(record["x2"]), f(record["y"])
            physical = physical and y > 1 and x1 * x1 + x2 * x2 + y * y > 1
            axes_avoided = axes_avoided and x1 != 0 and x2 != 0 and x1 != x2 and x1 != -x2
            point = (x1, x2, y)
            generator_fixed_sets_avoided = generator_fixed_sets_avoided and all(
                exact_h3_action(matrix, point) != point for matrix in GENERATORS.values()
            )
        except (KeyError, ValueError, ZeroDivisionError):
            physical = False
            axes_avoided = False
            generator_fixed_sets_avoided = False
    return {
        "candidate_coordinates_exact_once": len(coordinates) == len(set(coordinates)),
        "candidate_ids_exact_once": len(ids) == len(set(ids)),
        "all_points_physical": physical,
        "symmetry_axes_avoided": axes_avoided,
        "exact_generator_fixed_sets_avoided": generator_fixed_sets_avoided,
        "no_synthetic_rows": all(q.get("design") != "synthetic" for q in points),
    }


Pair = tuple[Fraction, Fraction]


def cadd(a: Pair, b: Pair) -> Pair:
    return a[0] + b[0], a[1] + b[1]


def cmul(a: Pair, b: Pair) -> Pair:
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def cconj(a: Pair) -> Pair:
    return a[0], -a[1]


def exact_h3_action(matrix: Any, point: tuple[Fraction, Fraction, Fraction]) -> tuple[Fraction, Fraction, Fraction]:
    x1, x2, y = point
    a, b, c, d = tuple((Fraction(q[0]), Fraction(q[1])) for q in matrix)
    z = (x1, x2)
    q = cadd(cmul(c, z), d)
    denominator = q[0] * q[0] + q[1] * q[1] + (c[0] * c[0] + c[1] * c[1]) * y * y
    numerator = cadd(cmul(cadd(cmul(a, z), b), cconj(q)), (
        (cmul(a, cconj(c))[0]) * y * y,
        (cmul(a, cconj(c))[1]) * y * y,
    ))
    return numerator[0] / denominator, numerator[1] / denominator, y / denominator


def frequency(cusp: str, mode: tuple[int, int, int]) -> Pair:
    a, b, _ = mode
    if cusp == "infinity":
        return Fraction(a), Fraction(b)
    return Fraction(2 * a - b, 5), Fraction(a + 2 * b, 5)


def effective_horizontal_samples(points: list[dict[str, Any]], cusp: str) -> list[Pair]:
    samples: set[Pair] = set()
    for record in points:
        point = f(record["x1"]), f(record["x2"]), f(record["y"])
        # Every source component and every exact target component occurring in
        # the complete physical row family contributes phase information.
        for copy in range(6):
            if (copy == 0) == (cusp == "infinity"):
                translate = Fraction(0 if copy == 0 else copy - 1)
                samples.add((point[0] + translate, point[1]))
        for generator, matrix in GENERATORS.items():
            transformed = exact_h3_action(matrix, point)
            for source_copy in range(6):
                target = GLUE[generator][source_copy]
                if (target == 0) == (cusp == "infinity"):
                    translate = Fraction(0 if target == 0 else target - 1)
                    samples.add((transformed[0] + translate, transformed[1]))
    return sorted(samples)


def phase_audit(cutoff: int, points: list[dict[str, Any]]) -> dict[str, Any]:
    modes_inf, modes_zero = cusp_modes(cutoff)
    per_cusp: dict[str, Any] = {}
    all_aliases = []
    maximum_coherence = 0.0
    worst_pair = None
    for cusp, modes in (("infinity", modes_inf), ("zero", modes_zero)):
        samples = effective_horizontal_samples(points, cusp)
        vectors = []
        for mode in modes:
            u, v = frequency(cusp, mode)
            vectors.append([
                cmath.exp(2j * math.pi * float(u * x1 + v * x2))
                for x1, x2 in samples
            ])
        aliases = []
        cusp_max = 0.0
        cusp_worst = None
        for i in range(len(modes)):
            ui, vi = frequency(cusp, modes[i])
            for j in range(i + 1, len(modes)):
                uj, vj = frequency(cusp, modes[j])
                exact_alias = all(
                    ((ui - uj) * x1 + (vi - vj) * x2).denominator == 1
                    for x1, x2 in samples
                )
                if exact_alias:
                    alias = {"cusp": cusp, "mode_a": list(modes[i]), "mode_b": list(modes[j])}
                    aliases.append(alias)
                    all_aliases.append(alias)
                inner = sum(a.conjugate() * b for a, b in zip(vectors[i], vectors[j]))
                coherence = abs(inner) / len(samples)
                if coherence > cusp_max:
                    cusp_max = coherence
                    cusp_worst = [list(modes[i]), list(modes[j])]
                if coherence > maximum_coherence:
                    maximum_coherence = coherence
                    worst_pair = {"cusp": cusp, "modes": cusp_worst}
        per_cusp[cusp] = {
            "sample_count": len(samples),
            "mode_count": len(modes),
            "exact_aliases": aliases,
            "maximum_phase_coherence_diagnostic": cusp_max,
            "worst_pair": cusp_worst,
        }
    return {
        "schema": "track-b-phase-alias-audit/v1",
        "fourier_cutoff": cutoff,
        "per_cusp": per_cusp,
        "exact_alias_count": len(all_aliases),
        "exact_aliases": all_aliases,
        "no_unexplained_exact_phase_aliases": len(all_aliases) == 0,
        "mu_phase_diagnostic": maximum_coherence,
        "worst_pair": worst_pair,
    }


def transformed_height_audit(points: list[dict[str, Any]]) -> dict[str, Any]:
    records = []
    signatures: dict[tuple[Fraction, ...], list[str]] = {}
    all_heights: set[Fraction] = set()
    for record in points:
        point = f(record["x1"]), f(record["x2"]), f(record["y"])
        transformed = []
        for name, matrix in GENERATORS.items():
            target = exact_h3_action(matrix, point)
            transformed.append(target[2])
            all_heights.add(target[2])
            records.append({
                "point_id": record["point_id"],
                "pairing_relation": name,
                "source_height": str(point[2]),
                "true_transformed_height": str(target[2]),
            })
        signature = (point[2], *transformed)
        signatures.setdefault(signature, []).append(record["point_id"])
    duplicate_signatures = [ids for ids in signatures.values() if len(ids) > 1]
    return {
        "schema": "track-b-transformed-height-audit/v1",
        "point_count": len(points),
        "unique_source_height_count": len({f(q["y"]) for q in points}),
        "unique_transformed_height_count": len(all_heights),
        "duplicate_complete_height_signatures": duplicate_signatures,
        "different_source_points_have_distinct_height_signatures": not duplicate_signatures,
        "records": records,
    }


def build_family(count: int, cutoff: int) -> dict[str, Any]:
    points = exact_candidate_points(count)
    point_checks = candidate_point_checks(points)
    phase = phase_audit(cutoff, points)
    height = transformed_height_audit(points)
    definition = {
        "schema": "track-b-exact-collocation-family/v1",
        "label": "NON-SYMMETRIC EXACT PHYSICAL COLLOCATION FAMILY",
        "fourier_cutoff_audited": cutoff,
        "points": points,
        "candidate_point_checks": point_checks,
        "candidate_pool_partitions": {
            "cusps": ["infinity", "zero"],
            "height_layers": [q["y"] for q in points],
            "horizontal_x1_coordinates": [q["x1"] for q in points],
            "horizontal_x2_coordinates": [q["x2"] for q in points],
            "six_copy_fibers": list(range(6)),
            "pairing_relations": list(GENERATORS),
            "row_partition_fields": [
                "source_cusp", "target_cusp", "source_height",
                "true_transformed_height", "source_copy", "target_copy",
                "group_element", "matrix_block",
            ],
        },
        "point_count": len(points),
        "row_expansion": "each point x four exact generators x six source fibers",
        "candidate_physical_row_count": 24 * len(points),
        "phase_alias_audit": phase,
        "transformed_height_audit": height,
        "geometric_conditions": {
            "all_y_gt_one": True,
            "all_inside_closed_humbert_core": True,
            "symmetry_axes_avoided": True,
            "exact_generator_fixed_sets_avoided": point_checks[
                "exact_generator_fixed_sets_avoided"
            ],
            "exact_rational_coordinates": True,
        },
    }
    definition["collocation_family_hash"] = canonical_hash(definition)
    return definition


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=20)
    parser.add_argument("--fourier-cutoff", type=int, default=5)
    parser.add_argument("--json-out", type=Path, default=ROOT / "track_b_collocation_family_result.json")
    ns = parser.parse_args()
    result = build_family(ns.count, ns.fourier_cutoff)
    ns.json_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "point_count": result["point_count"],
        "candidate_physical_row_count": result["candidate_physical_row_count"],
        "exact_alias_count": result["phase_alias_audit"]["exact_alias_count"],
        "mu_phase_diagnostic": result["phase_alias_audit"]["mu_phase_diagnostic"],
        "unique_transformed_height_count": result["transformed_height_audit"]["unique_transformed_height_count"],
        "collocation_family_hash": result["collocation_family_hash"],
    }, indent=2))
    return 0 if result["phase_alias_audit"]["no_unexplained_exact_phase_aliases"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
