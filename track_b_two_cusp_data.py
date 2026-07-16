#!/usr/bin/env python3
"""Exact cusp and six-copy group data for the Track-B physical Hejhal system."""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from typing import Any, Sequence

from flint import acb, arb


GI = tuple[int, int]
Matrix = tuple[GI, GI, GI, GI]

I: GI = (0, 1)
ZERO: GI = (0, 0)
ONE: GI = (1, 0)

GENERATORS: dict[str, Matrix] = {
    "T1": (ONE, ONE, ZERO, ONE),
    "R": (I, ZERO, ZERO, (0, -1)),
    "TiR": (I, ONE, ZERO, (0, -1)),
    "S": (ZERO, (-1, 0), ONE, ZERO),
}

# pi_gamma(c)=c.gamma^{-1} in the project's six-copy convention.
GLUE: dict[str, tuple[int, ...]] = {
    "T1": (0, 5, 1, 2, 3, 4),
    "R": (0, 1, 5, 4, 3, 2),
    "TiR": (0, 3, 2, 1, 5, 4),
    "S": (1, 0, 5, 3, 4, 2),
}

SIGMA_INFINITY: Matrix = (ONE, ZERO, ZERO, ONE)
SIGMA_ZERO: Matrix = GENERATORS["S"]
LEVEL: GI = (2, 1)


def canonical_hash(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def matrix_json(matrix: Matrix) -> list[list[list[int]]]:
    a, b, c, d = matrix
    return [[list(a), list(b)], [list(c), list(d)]]


def cusp_data() -> dict[str, Any]:
    data = {
        "schema": "track-b-exact-two-cusp-data/v1",
        "level": [2, 1],
        "index": 6,
        "scaling_matrices": {
            "infinity": matrix_json(SIGMA_INFINITY),
            "zero": matrix_json(SIGMA_ZERO),
        },
        "scaling_convention": (
            "sigma_infinity=I; sigma_zero=S; local cusp coordinate is "
            "sigma_a^{-1}P.  In PSL2, S^{-1}=S."
        ),
        "cusp_lattices": {
            "infinity": "Z[i]",
            "zero": "(2+i) Z[i]",
            "zero_dual": "(1/conj(2+i)) Z[i] = ((2+i)/5) Z[i]",
            "zero_covolume": "5",
        },
        "coset_representatives": [
            "I", "S", "S T_1", "S T_2", "S T_3", "S T_4"
        ],
        "component_cusp_classes": ["infinity", "zero", "zero", "zero", "zero", "zero"],
        "component_zero_translates": [None, 0, 1, 2, 3, 4],
        "generators": {name: matrix_json(matrix) for name, matrix in GENERATORS.items()},
        "glue": {name: list(values) for name, values in GLUE.items()},
    }
    data["cusp_data_hash"] = canonical_hash(data)
    return data


def qcomplex(value: GI) -> acb:
    return acb(arb(value[0]), arb(value[1]))


def h3_action(
    matrix: Matrix, x1: arb, x2: arb, y: arb
) -> tuple[arb, arb, arb]:
    """Validated upper-half-space action of an exact Gaussian matrix."""
    a, b, c, d = (qcomplex(q) for q in matrix)
    z = acb(x1, x2)
    q = c * z + d
    denominator = q.real * q.real + q.imag * q.imag + abs(c) ** 2 * y * y
    if not bool(denominator > 0):
        raise ArithmeticError(f"Möbius denominator is not positive: {denominator}")
    numerator = (a * z + b) * q.conjugate() + a * c.conjugate() * y * y
    zp = numerator / denominator
    yp = y / denominator
    if not zp.is_finite() or not yp.is_finite():
        raise ArithmeticError("nonfinite exact H3 action enclosure")
    return zp.real, zp.imag, yp


def sigma_zero_specialized(x1: arb, x2: arb, y: arb) -> tuple[arb, arb, arb]:
    """Specialized S^{-1}=S cusp-zero coordinate formula."""
    rho = x1 * x1 + x2 * x2 + y * y
    if not bool(rho > 0):
        raise ArithmeticError(f"sigma0 coordinate denominator not positive: {rho}")
    return -x1 / rho, x2 / rho, y / rho


def gaussian_modes(cutoff: int) -> list[tuple[int, int, int]]:
    if cutoff < 1:
        raise ValueError("Fourier cutoff must be positive")
    out = []
    radius = math.isqrt(cutoff)
    for a in range(-radius, radius + 1):
        for b in range(-radius, radius + 1):
            norm = a * a + b * b
            if 0 < norm <= cutoff:
                out.append((a, b, norm))
    return sorted(out, key=lambda row: (row[2], row[0], row[1]))


def cusp_modes(cutoff: int) -> tuple[list[tuple[int, int, int]], list[tuple[int, int, int]]]:
    # |(2+i) beta / 5|^2=N(beta)/5, so N(beta)<=5M gives the
    # same physical frequency cutoff at cusp zero.
    return gaussian_modes(cutoff), gaussian_modes(5 * cutoff)


def mode_frequency(
    cusp: str, mode: Sequence[int]
) -> tuple[arb, arb, arb]:
    a, b, norm = (int(q) for q in mode)
    if cusp == "infinity":
        return arb(a), arb(b), arb(norm).sqrt()
    if cusp != "zero":
        raise ValueError(f"unknown cusp {cusp!r}")
    # (a+bi)(2+i)/5=((2a-b)+i(a+2b))/5.
    u = arb(2 * a - b) / 5
    v = arb(a + 2 * b) / 5
    return u, v, (u * u + v * v).sqrt()


_POINTS = (
    ((1, 10), (1, 8), (9, 8)),
    ((-1, 7), (2, 11), (7, 6)),
    ((2, 9), (1, 5), (6, 5)),
    ((-2, 9), (3, 10), (11, 10)),
    ((1, 6), (4, 13), (23, 20)),
    ((-3, 14), (1, 9), (19, 16)),
)


def exact_collocation_points(count: int) -> list[dict[str, str]]:
    if count < 1 or count > len(_POINTS):
        raise ValueError(f"collocation point count must be in [1,{len(_POINTS)}]")
    rows = []
    for point in _POINTS[:count]:
        values = [Fraction(num, den) for num, den in point]
        if values[0] * values[0] + values[1] * values[1] + values[2] * values[2] <= 1:
            raise ArithmeticError("collocation point is outside the closed Humbert core")
        rows.append({
            "x1": str(values[0]), "x2": str(values[1]), "y": str(values[2])
        })
    return rows


def arb_fraction(text: str) -> arb:
    value = Fraction(text)
    return arb(value.numerator) / value.denominator

