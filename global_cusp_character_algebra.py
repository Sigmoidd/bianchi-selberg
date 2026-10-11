#!/usr/bin/env python3
r"""Exact finite algebra for the Gamma_0(2+i) six-copy representation.

Status of everything in this file: PROVED (finite exact arithmetic over Q and
Q(zeta_5)).  No analytic statement is made here.

Conventions (identical to ``six_copy_DK_algebra.py``):

* ``GLUE[g][c] = c . g^{-1}`` (stored inverse right action);
* the *forward* action is ``fwd[g] = inverse(GLUE[g])``, ``c -> c . g``;
* the pullback representation is ``(rho(g) v)_c = v_{c.g}``, so
  ``rho(g) rho(h) = rho(gh)`` and the matrix of a forward permutation ``p`` has
  entries ``rho[c][p(c)] = 1``.

Coset labels are *derived* here, not assumed: label 0 is H_pi.I, label 1 is
H_pi.S and label 1+k (k=1..4) is H_pi.S.T_1^k.  The bottom row (c1 : c2) of a
coset representative modulo pi is the point t = c1/c2 of P^1(F_5).  The module
recomputes all four generator permutations from Gaussian matrices and asserts
that they equal the repository's arrays.
"""
from __future__ import annotations

from fractions import Fraction
import itertools
import json
from pathlib import Path
from typing import Sequence

from global_cusp_gaussian import (
    GENERATORS, QMAT, PI, PI_BAR, G, Mat, g, gmul, gneg, mmul, minv, mdet,
    random_word, is_gaussian_integer, divisible_by, I as GI_I,
)
from six_copy_DK_algebra import (
    GLUE, N, compose, generated_group, identity, inverse, orbits,
    perm_matrix_for_pullback, projector_for_orbits, matmul,
)

Perm = tuple[int, ...]
FMat = list[list[Fraction]]

FWD: dict[str, Perm] = {name: inverse(p) for name, p in GLUE.items()}


# ------------------------------------------------------------ Q(zeta_5) arithmetic

class Zeta5:
    """Element of Q(zeta_5) as coefficients on 1, z, z^2, z^3 (z^4 = -1-z-z^2-z^3)."""

    __slots__ = ("c",)

    def __init__(self, c: Sequence[int | Fraction]):
        if len(c) != 4:
            raise ValueError("need four coefficients")
        self.c = tuple(Fraction(x) for x in c)

    @staticmethod
    def zeta_pow(j: int) -> "Zeta5":
        j %= 5
        if j < 4:
            v = [0, 0, 0, 0]
            v[j] = 1
            return Zeta5(v)
        return Zeta5([-1, -1, -1, -1])

    @staticmethod
    def const(x: int | Fraction) -> "Zeta5":
        return Zeta5([x, 0, 0, 0])

    def __add__(self, other: "Zeta5") -> "Zeta5":
        return Zeta5([a + b for a, b in zip(self.c, other.c)])

    def __mul__(self, other: "Zeta5") -> "Zeta5":
        folded = [Fraction(0)] * 5
        for i, a in enumerate(self.c):
            for j, b in enumerate(other.c):
                folded[(i + j) % 5] += a * b          # x^5 = 1
        top = folded[4]                                # x^4 = -(1+x+x^2+x^3)
        return Zeta5([folded[i] - top for i in range(4)])

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Zeta5) and self.c == other.c

    def is_zero(self) -> bool:
        return all(x == 0 for x in self.c)

    def __repr__(self) -> str:
        return f"Zeta5{tuple(str(x) for x in self.c)}"


def zeta_sums() -> dict[int, bool]:
    """sum_{j=0}^{4} zeta^{kj}: equals 5 for k=0 and 0 for k=1..4 (exactly)."""
    out = {}
    for k in range(5):
        s = Zeta5.const(0)
        for j in range(5):
            s = s + Zeta5.zeta_pow(k * j)
        out[k] = (s == Zeta5.const(5)) if k == 0 else s.is_zero()
    return out


def character_vector(k: int) -> list[Zeta5]:
    """w_k: component 0 is 0; component 1+m is zeta^{k m} (m = 0..4)."""
    return [Zeta5.const(0)] + [Zeta5.zeta_pow(k * m) for m in range(5)]


def apply_rho(perm: Perm, vec: list[Zeta5]) -> list[Zeta5]:
    """(rho(p) v)_c = v_{p(c)} for a forward permutation p."""
    return [vec[perm[c]] for c in range(N)]


def scale_vec(s: Zeta5, vec: list[Zeta5]) -> list[Zeta5]:
    return [s * x for x in vec]


def vec_eq(a: list[Zeta5], b: list[Zeta5]) -> bool:
    return all(x == y for x, y in zip(a, b))


# ------------------------------------------------------------ permutation helpers

def perm_pow(p: Perm, n: int) -> Perm:
    out = identity()
    for _ in range(n % _perm_order(p)):
        out = compose(p, out)
    return out


def _perm_order(p: Perm) -> int:
    k, cur = 1, p
    while cur != identity():
        cur = compose(p, cur)
        k += 1
    return k


def conj_perm(q: Perm, p: Perm) -> Perm:
    """q o p o q^{-1} as functions."""
    return compose(compose(q, p), inverse(q))


def right_action_perm(word: Sequence[str]) -> Perm:
    """Forward permutation of the product g1 g2 ... (c.(g1 g2) = (c.g1).g2)."""
    out = identity()
    for name in word:
        out = compose(FWD[name], out)      # apply earlier letter first
    return out


# ------------------------------------------------------------ first-principles P^1(F_5)

def i_residue(prime: str) -> int:
    """Image of i in Z[i]/prime = F_5."""
    return {"pi": 3, "pibar": 2}[prime]


def _red(z: G, prime: str) -> int:
    if not is_gaussian_integer(z):
        raise ValueError("not a Gaussian integer")
    return (int(z[0]) + i_residue(prime) * int(z[1])) % 5


def _inv5(x: int) -> int:
    return pow(x % 5, -1, 5)


def canon_point(c1: int, c2: int):
    c1, c2 = c1 % 5, c2 % 5
    if (c1, c2) == (0, 0):
        raise ValueError("zero row")
    if c2 != 0:
        return (c1 * _inv5(c2)) % 5
    return "inf"


def rep_row(label: int) -> tuple[int, int]:
    """Bottom row of the coset representative: label 0 -> I, label 1+k -> S T1^k."""
    if label == 0:
        return (0, 1)
    return (1, label - 1)


def label_of_point(t) -> int:
    if t == 0:
        return 0
    if t == "inf":
        return 1
    return 1 + _inv5(t)          # t = 1/k  =>  k = t^{-1}


def forward_perm_from_matrix(m: Mat, prime: str = "pi") -> Perm:
    """c -> c.m on P^1(F_5) = H_prime\\Gamma, from the Gaussian matrix itself."""
    a, b, c, d = (_red(x, prime) for x in m)
    out = []
    for label in range(N):
        r1, r2 = rep_row(label)
        out.append(label_of_point(canon_point(r1 * a + r2 * c, r1 * b + r2 * d)))
    return tuple(out)


def first_principles_generator_perms(prime: str = "pi") -> dict[str, Perm]:
    return {name: forward_perm_from_matrix(m, prime) for name, m in GENERATORS.items()}


def quarter_turn_perm(prime: str = "pi") -> Perm:
    """Label permutation sigma with Q gamma_c Q^{-1} in H_prime gamma_{sigma(c)}.

    Q = diag(i,1) conjugates [[a,b],[c1,c2]] to [[a, i b],[-i c1, c2]], so the
    bottom row (c1 : c2) goes to (-i c1 : c2).
    """
    iv = i_residue(prime)
    out = []
    for label in range(N):
        r1, r2 = rep_row(label)
        out.append(label_of_point(canon_point(-iv * r1, r2)))
    return tuple(out)


# ------------------------------------------------------------ exact linear algebra over Q

def to_fmat(rows: Sequence[Sequence[int | Fraction]]) -> FMat:
    return [[Fraction(x) for x in r] for r in rows]


def rho_matrix(perm: Perm) -> FMat:
    return perm_matrix_for_pullback(perm)


def mat_add(a: FMat, b: FMat) -> FMat:
    return [[x + y for x, y in zip(r, s)] for r, s in zip(a, b)]


def mat_sub(a: FMat, b: FMat) -> FMat:
    return [[x - y for x, y in zip(r, s)] for r, s in zip(a, b)]


def mat_scale(a: FMat, s: Fraction) -> FMat:
    return [[s * x for x in r] for r in a]


def mat_t(a: FMat) -> FMat:
    return [list(r) for r in zip(*a)]


def ident_f() -> FMat:
    return [[Fraction(int(i == j)) for j in range(N)] for i in range(N)]


def is_zero_f(a: FMat) -> bool:
    return all(x == 0 for r in a for x in r)


def rank_q(rows: FMat) -> int:
    m = [list(r) for r in rows]
    rank, ncols = 0, len(m[0]) if m else 0
    for col in range(ncols):
        piv = next((r for r in range(rank, len(m)) if m[r][col] != 0), None)
        if piv is None:
            continue
        m[rank], m[piv] = m[piv], m[rank]
        pv = m[rank][col]
        m[rank] = [x / pv for x in m[rank]]
        for r in range(len(m)):
            if r != rank and m[r][col] != 0:
                f = m[r][col]
                m[r] = [x - f * y for x, y in zip(m[r], m[rank])]
        rank += 1
    return rank


def nullspace_q(rows: FMat, ncols: int) -> list[list[Fraction]]:
    m = [list(r) for r in rows]
    pivots: list[int] = []
    rank = 0
    for col in range(ncols):
        piv = next((r for r in range(rank, len(m)) if m[r][col] != 0), None)
        if piv is None:
            continue
        m[rank], m[piv] = m[piv], m[rank]
        pv = m[rank][col]
        m[rank] = [x / pv for x in m[rank]]
        for r in range(len(m)):
            if r != rank and m[r][col] != 0:
                f = m[r][col]
                m[r] = [x - f * y for x, y in zip(m[r], m[rank])]
        pivots.append(col)
        rank += 1
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for f in free:
        v = [Fraction(0)] * ncols
        v[f] = Fraction(1)
        for r, pc in enumerate(pivots):
            v[pc] = -m[r][f]
        basis.append(v)
    return basis


def commutant_basis(perms: Sequence[Perm]) -> list[FMat]:
    """Exact rational basis of {M : M rho(p) = rho(p) M for all p}."""
    eqs: FMat = []
    for p in perms:
        rho = rho_matrix(p)
        for i in range(N):
            for j in range(N):
                row = [Fraction(0)] * (N * N)
                # (M rho - rho M)_{ij} = sum_k M_ik rho_kj - rho_ik M_kj
                for k in range(N):
                    row[i * N + k] += rho[k][j]
                    row[k * N + j] -= rho[i][k]
                eqs.append(row)
    return [[vec[i * N:(i + 1) * N] for i in range(N)] for vec in nullspace_q(eqs, N * N)]


# ------------------------------------------------------------ the certificate

def local_projector() -> FMat:
    """A_T = (1/5) sum_{j=0}^{4} rho(T1)^j  (the cusp-translation projector)."""
    rho_t = rho_matrix(FWD["T1"])
    acc = [[Fraction(0)] * N for _ in range(N)]
    cur = ident_f()
    for _ in range(5):
        acc = mat_add(acc, cur)
        cur = matmul(cur, rho_t)
    return mat_scale(acc, Fraction(1, 5))


def class_sizes(group: set[Perm]) -> list[int]:
    seen: set[Perm] = set()
    sizes = []
    elems = list(group)
    inv_cache = {x: inverse(x) for x in elems}
    for x in elems:
        if x in seen:
            continue
        cl = {compose(compose(h, x), inv_cache[h]) for h in elems}
        seen |= cl
        sizes.append(len(cl))
    return sorted(sizes)


def no_proper_normal_subgroup_from_classes(sizes: Sequence[int], order: int) -> bool:
    """A normal subgroup is a union of classes containing {1}; its order divides |G|."""
    rest = [s for s in sizes if s != 1]
    ones = sizes.count(1)
    for r in range(len(rest) + 1):
        for sub in itertools.combinations(range(len(rest)), r):
            tot = ones + sum(rest[i] for i in sub)
            if tot not in (1, order) and order % tot == 0:
                return False
    return True


def compute_certificate() -> dict:
    cert: dict = {"schema": "global-cusp-character-algebra/v1", "status": "PROVED (exact finite arithmetic)"}

    # 0. conventions re-derived from Gaussian matrices
    fp = first_principles_generator_perms("pi")
    cert["first_principles_perms_match_repository"] = {n: fp[n] == FWD[n] for n in FWD}
    cert["Ti_forward_from_matrix"] = list(forward_perm_from_matrix(GENERATORS["TiR"], "pi"))

    # 1. cyclotomic sums and the character basis
    cert["zeta_sums_exact"] = zeta_sums()
    t1 = FWD["T1"]
    eig = {}
    for k in range(1, 5):
        w = character_vector(k)
        eig[k] = vec_eq(apply_rho(t1, w), scale_vec(Zeta5.zeta_pow(k), w))
    cert["T1_eigenvector_equations"] = eig
    r = FWD["R"]
    cert["R_exchanges_k_and_5_minus_k"] = {
        k: vec_eq(apply_rho(r, character_vector(k)), character_vector(5 - k)) for k in range(1, 5)
    }

    # 2. rational projector
    a = local_projector()
    cert["A_T_idempotent"] = matmul(a, a) == a
    cert["A_T_symmetric"] = mat_t(a) == a
    cert["A_T_rank"] = rank_q(a)
    cert["I_minus_A_T_rank"] = rank_q(mat_sub(ident_f(), a))
    cert["A_T_equals_repository_singular_projector"] = a == projector_for_orbits([(0,), (1, 2, 3, 4, 5)])

    # 3. full group
    full = generated_group(tuple(FWD.values()))
    cert["group_order"] = len(full)
    sizes = class_sizes(full)
    cert["class_sizes"] = sizes
    cert["class_sizes_equal_A5"] = sizes == [1, 12, 12, 15, 20]
    cert["no_proper_normal_subgroup_by_class_union"] = no_proper_normal_subgroup_from_classes(sizes, len(full))
    cert["all_elements_even_on_6_points"] = all(_sign(p) == 1 for p in full)
    pairs = {(a_, b_) for a_ in range(N) for b_ in range(N) if a_ != b_}
    orbit_of_01 = {(p[0], p[1]) for p in full}
    cert["doubly_transitive"] = orbit_of_01 == pairs
    cert["stabilizer_of_a_point_order"] = len([p for p in full if p[0] == 0])

    # 4. character of the permutation representation and commutant
    chi = [sum(1 for c in range(N) if p[c] == c) for p in full]
    cert["inner_product_chi_chi"] = str(Fraction(sum(x * x for x in chi), len(full)))
    cert["inner_product_chi_trivial"] = str(Fraction(sum(chi), len(full)))
    comm = commutant_basis(list(FWD.values()))
    cert["commutant_dimension"] = len(comm)
    j_mat = [[Fraction(1)] * N for _ in range(N)]
    span_ok = all(
        rank_q([[x for r in m for x in r] for m in comm] + [[x for r in t for x in r]]) == len(comm)
        for t in (ident_f(), j_mat)
    )
    cert["commutant_equals_span_I_J"] = (len(comm) == 2) and span_ok

    # 5. noncommutation of the local projector with the global generators
    cert["A_T_commutes_with"] = {
        name: (matmul(a, rho_matrix(p)) == matmul(rho_matrix(p), a)) for name, p in FWD.items()
    }
    comm_s = mat_sub(matmul(a, rho_matrix(FWD["S"])), matmul(rho_matrix(FWD["S"]), a))
    cert["commutator_A_T_rho_S_nonzero"] = not is_zero_f(comm_s)
    cert["commutator_A_T_rho_S"] = [[str(x) for x in row] for row in comm_s]

    # 6. regular subspace V4 = Ran(I - A_T) is not invariant; its orbit spans V5
    v4 = [[(1 if c == m else 0) - a[c][m] for m in range(N)] for c in range(N)]  # columns span V4
    v4_cols = mat_t(v4)
    span_rows = [list(col) for col in v4_cols]
    for p in full:
        rho = rho_matrix(p)
        span_rows += [[sum(rho[i][k] * col[k] for k in range(N)) for i in range(N)] for col in v4_cols]
    cert["V4_dimension"] = rank_q(v4_cols)
    cert["dim_span_of_group_orbit_of_V4"] = rank_q(span_rows)
    cert["V4_group_invariant"] = rank_q(span_rows) == rank_q(v4_cols)

    # 7. obstruction: no nonzero constant equivariant E annihilates both singular directions
    #    E = x I + y J ; singular directions e_0 and u = e_1+...+e_5.
    sing = [[1, 0, 0, 0, 0, 0], [0, 1, 1, 1, 1, 1]]
    eqs = []
    for v in sing:
        for i in range(N):
            # (x I + y J) v  component i :  x v_i + y sum(v)
            eqs.append([Fraction(v[i]), Fraction(sum(v))])
    ns = nullspace_q(eqs, 2)
    cert["constant_equivariant_annihilators_of_singular_subspace"] = len(ns)

    # 8. the quarter-turn finite permutation (derived, then compared with the candidate)
    q = quarter_turn_perm("pi")
    cert["quarter_turn_perm_derived"] = list(q)
    cert["quarter_turn_perm_equals_candidate"] = q == (0, 1, 4, 2, 5, 3)
    cert["q_squared_equals_rho_R"] = compose(q, q) == FWD["R"]
    cert["q_T1_qinv_equals_T1_cubed_as_permutations"] = conj_perm(q, FWD["T1"]) == perm_pow(FWD["T1"], 3)
    rho_q = rho_matrix(q)
    rho_t = rho_matrix(FWD["T1"])
    rho_t2 = matmul(rho_t, rho_t)
    rho_t3 = matmul(rho_t2, rho_t)
    qi = rho_matrix(inverse(q))
    cert["matrix_level_rho_q_rho_T1_rho_qinv"] = {
        "equals_rho_T1_squared": matmul(matmul(rho_q, rho_t), qi) == rho_t2,
        "equals_rho_T1_cubed": matmul(matmul(rho_q, rho_t), qi) == rho_t3,
    }
    normalizes = {name: conj_perm(q, p) in full for name, p in FWD.items()}
    cert["q_conjugates_each_generator_into_the_order_60_group"] = normalizes
    cert["q_in_order_60_group"] = q in full
    big = generated_group(tuple(FWD.values()) + (q,))
    cert["group_generated_with_q_order"] = len(big)
    cert["q_is_odd_permutation"] = _sign(q) == -1
    # q acts on characters as k -> 3k (a 4-cycle on {1,2,3,4})
    cert["rho_q_w_k_equals_w_3k"] = {
        k: vec_eq(apply_rho(q, character_vector(k)), character_vector(3 * k % 5)) for k in range(1, 5)
    }
    orbit = [1]
    while (orbit[-1] * 3) % 5 != 1:
        orbit.append((orbit[-1] * 3) % 5)
    cert["q_orbit_on_nontrivial_characters"] = orbit

    cert["all_checks_pass"] = _all_pass(cert)
    return cert


def _sign(p: Perm) -> int:
    seen, sgn = set(), 1
    for i in range(len(p)):
        if i in seen:
            continue
        j, ln = i, 0
        while j not in seen:
            seen.add(j)
            j = p[j]
            ln += 1
        if ln % 2 == 0:
            sgn = -sgn
    return sgn


def _all_pass(c: dict) -> bool:
    return (
        all(c["first_principles_perms_match_repository"].values())
        and all(c["zeta_sums_exact"].values())
        and all(c["T1_eigenvector_equations"].values())
        and all(c["R_exchanges_k_and_5_minus_k"].values())
        and c["A_T_idempotent"] and c["A_T_symmetric"]
        and c["A_T_rank"] == 2 and c["I_minus_A_T_rank"] == 4
        and c["A_T_equals_repository_singular_projector"]
        and c["group_order"] == 60 and c["class_sizes_equal_A5"]
        and c["no_proper_normal_subgroup_by_class_union"] and c["all_elements_even_on_6_points"]
        and c["doubly_transitive"] and c["inner_product_chi_chi"] == "2"
        and c["inner_product_chi_trivial"] == "1"
        and c["commutant_dimension"] == 2 and c["commutant_equals_span_I_J"]
        and c["commutator_A_T_rho_S_nonzero"] and not c["V4_group_invariant"]
        and c["dim_span_of_group_orbit_of_V4"] == 5
        and c["constant_equivariant_annihilators_of_singular_subspace"] == 0
        and c["quarter_turn_perm_equals_candidate"] and c["q_squared_equals_rho_R"]
        and c["q_T1_qinv_equals_T1_cubed_as_permutations"]
        and c["matrix_level_rho_q_rho_T1_rho_qinv"]["equals_rho_T1_squared"]
        and not c["matrix_level_rho_q_rho_T1_rho_qinv"]["equals_rho_T1_cubed"]
        and all(c["q_conjugates_each_generator_into_the_order_60_group"].values())
        and not c["q_in_order_60_group"] and c["group_generated_with_q_order"] == 120
        and c["q_is_odd_permutation"]
        and all(c["rho_q_w_k_equals_w_3k"].values())
        and sorted(c["q_orbit_on_nontrivial_characters"]) == [1, 2, 3, 4]
    )


def main() -> int:
    cert = compute_certificate()
    out = Path("global_cusp_character_algebra_result.json")
    out.write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({k: cert[k] for k in ("status", "all_checks_pass", "group_order", "A_T_rank")}, indent=2))
    return 0 if cert["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
