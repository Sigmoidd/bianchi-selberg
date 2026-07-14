#!/usr/bin/env python3
r"""Exact algebra checks for the six-copy D(K) extension.

The gluing arrays store the inverse right action

    pi_g(c) = c . g^{-1}

on Gamma_0(2+i)\Gamma.  The induced representation used by the
vector-valued theorem is therefore

    (rho(g) v)_c = v_{c.g} = v_{pi_g^{-1}(c)}.

Everything in this file is finite permutation/rational arithmetic.  It does
not use floating point and it does not certify any analytic defect.
"""
from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
from typing import Iterable


GLUE: dict[str, tuple[int, ...]] = {
    "T1": (0, 5, 1, 2, 3, 4),
    "R": (0, 1, 5, 4, 3, 2),
    "TiR": (0, 3, 2, 1, 5, 4),
    "S": (1, 0, 5, 3, 4, 2),
}
N = 6


def identity() -> tuple[int, ...]:
    return tuple(range(N))


def compose(p: tuple[int, ...], q: tuple[int, ...]) -> tuple[int, ...]:
    """Return p after q: (p o q)(c)=p(q(c))."""
    return tuple(p[q[c]] for c in range(N))


def inverse(p: tuple[int, ...]) -> tuple[int, ...]:
    q = [0] * N
    for c, pc in enumerate(p):
        q[pc] = c
    return tuple(q)


def generated_group(generators: Iterable[tuple[int, ...]]) -> set[tuple[int, ...]]:
    gens = tuple(generators)
    group = {identity()}
    frontier = [identity()]
    while frontier:
        a = frontier.pop()
        for g in gens:
            for b in (compose(a, g), compose(g, a)):
                if b not in group:
                    group.add(b)
                    frontier.append(b)
    return group


def orbits(group: Iterable[tuple[int, ...]]) -> list[tuple[int, ...]]:
    group = tuple(group)
    unseen = set(range(N))
    out: list[tuple[int, ...]] = []
    while unseen:
        seed = min(unseen)
        orbit = {g[seed] for g in group}
        out.append(tuple(sorted(orbit)))
        unseen -= orbit
    return out


def projector_for_orbits(parts: list[tuple[int, ...]]) -> list[list[Fraction]]:
    p = [[Fraction(0) for _ in range(N)] for _ in range(N)]
    for orbit in parts:
        w = Fraction(1, len(orbit))
        for i in orbit:
            for j in orbit:
                p[i][j] = w
    return p


def matmul(a: list[list[Fraction]], b: list[list[Fraction]]) -> list[list[Fraction]]:
    return [
        [sum(a[i][k] * b[k][j] for k in range(N)) for j in range(N)]
        for i in range(N)
    ]


def perm_matrix_for_pullback(forward_action: tuple[int, ...]) -> list[list[Fraction]]:
    """Matrix rho with (rho v)_c = v_{forward_action[c]}."""
    return [
        [Fraction(int(j == forward_action[i])) for j in range(N)]
        for i in range(N)
    ]


def fraction_matrix_text(a: list[list[Fraction]]) -> list[list[str]]:
    return [[str(x) for x in row] for row in a]


def main() -> int:
    for name, p in GLUE.items():
        if tuple(sorted(p)) != identity():
            raise AssertionError(f"{name} is not a permutation")

    # TiR = Ti * R.  Forward right actions are anti-homomorphic:
    # f_{TiR}=f_R o f_Ti, hence f_Ti=f_R^{-1} o f_{TiR}.  Convert the stored
    # inverse actions to forward actions before composing.
    forward = {name: inverse(p) for name, p in GLUE.items()}
    ti_forward = compose(inverse(forward["R"]), forward["TiR"])

    parabolic_group = generated_group((forward["T1"], ti_forward, forward["R"]))
    translation_group = generated_group((forward["T1"], ti_forward))
    t1_r_group = generated_group((forward["T1"], forward["R"]))
    full_group = generated_group(tuple(forward.values()))
    parabolic_orbits = orbits(parabolic_group)
    translation_orbits = orbits(translation_group)
    full_orbits = orbits(full_group)
    expected = [(0,), (1, 2, 3, 4, 5)]
    if parabolic_orbits != expected:
        raise AssertionError(f"unexpected parabolic orbits {parabolic_orbits}")
    if translation_orbits != expected:
        raise AssertionError(f"unexpected translation orbits {translation_orbits}")
    if t1_r_group != parabolic_group:
        raise AssertionError("Ti/TiR unexpectedly enlarges <T1,R>")
    if forward["TiR"] not in t1_r_group:
        raise AssertionError("TiR is not reproduced by <T1,R>")
    if len(full_group) != 60 or full_orbits != [tuple(range(N))]:
        raise AssertionError(
            f"unexpected full image: order={len(full_group)}, orbits={full_orbits}"
        )

    singular_projector = projector_for_orbits(parabolic_orbits)
    if matmul(singular_projector, singular_projector) != singular_projector:
        raise AssertionError("singular projector is not idempotent")
    for g in parabolic_group:
        rho = perm_matrix_for_pullback(g)
        if matmul(rho, singular_projector) != matmul(singular_projector, rho):
            raise AssertionError("singular projector does not commute with cusp holonomy")

    # The T1 action is a 5-cycle on the finite-cusp orbit.  This identifies
    # the four nontrivial regular character channels in addition to the two
    # invariant (singular) channels.
    t1_cycles = orbits(generated_group((forward["T1"],)))
    if t1_cycles != expected:
        raise AssertionError(f"unexpected T1 cycles {t1_cycles}")

    result = {
        "status": "exact finite permutation/rational verification",
        "glue_convention": "pi_g(c)=c.g^{-1}",
        "representation_convention": "(rho(g)v)_c=v_{c.g}=v_{pi_g^{-1}(c)}",
        "forward_actions": {k: list(v) for k, v in forward.items()},
        "Ti_forward_derived_from_TiR_R_inverse": list(ti_forward),
        "translation_group_order": len(translation_group),
        "parabolic_group_order": len(parabolic_group),
        "T1_R_group_order": len(t1_r_group),
        "TiR_is_redundant_in_T1_R_group": forward["TiR"] in t1_r_group,
        "full_group_order": len(full_group),
        "full_group_orbits": [list(q) for q in full_orbits],
        "full_group_is_transitive": full_orbits == [tuple(range(N))],
        "translation_orbits": [list(q) for q in translation_orbits],
        "parabolic_orbits": [list(q) for q in parabolic_orbits],
        "singular_degree": len(parabolic_orbits),
        "singular_basis": [
            ["1", "0", "0", "0", "0", "0"],
            ["0", "1/sqrt(5)", "1/sqrt(5)", "1/sqrt(5)", "1/sqrt(5)", "1/sqrt(5)"],
        ],
        "singular_projector": fraction_matrix_text(singular_projector),
        "regular_channel_dimension": N - len(parabolic_orbits),
        "T1_orbits": [list(q) for q in t1_cycles],
        "proved_here": [
            "all four glue arrays are permutations",
            "the induced pullback matrices are unitary permutation matrices",
            "the base-cusp parabolic holonomy has exactly two orbits",
            "the singular subspace has dimension two",
            "the translation-fixed and full-parabolic-fixed subspaces coincide",
            "the remaining four dimensions are regular cusp channels",
            "TiR is already in <T1,R> and is only a cross-check transition",
            "the full four-generator image has order 60 and is transitive",
        ],
        "not_proved_here": [
            "analytic automorphy defect bounds",
            "partition-of-unity derivative constants",
            "cuspidal quasimode residual or eigenvalue enclosure",
        ],
    }
    out = Path("six_copy_DK_algebra_result.json")
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(out.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
