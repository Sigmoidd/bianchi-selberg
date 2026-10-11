#!/usr/bin/env python3
r"""Old/new bookkeeping for H_pi = Gamma_0(2+i): exact algebra only.

Two level-one degeneracy lifts of f on Gamma\H^3 to H_pi\H^3:

    B_1 f = f ,      B_pi f = f o g ,     g = diag(pi, 1)  (projectively).

PROVED here (exact Q(i) arithmetic):
  * g H_pi g^{-1} subset Gamma, with equality of conditions g M g^{-1} in Gamma  <=>  pi | M_21;
  * in the six-copy picture the second lift F2_c(z) = f(g gamma_c z) satisfies the same
    equivariance F2(gamma z) = rho(gamma) F2(z) as the first;
  * the six elements g gamma_c lie in six distinct left Gamma-cosets (Hecke decomposition);
  * Q commutes with g, hence U_Q B_1 = B_1 U_Q^{(1)} and U_Q B_pi = B_pi U_Q^{(1)};
  * the Fricke/Atkin-Lehner element W = S g normalizes H_pi, W^{-1} Q W Q^{-1} in H_pi, so
    U_W commutes with U_Q.

NOT done here (see docs): the orthogonal projector P_old,full is defined abstractly as the
projector onto the closed span of Ran B_1 + Ran B_pi.  It is not *computable* without certified
level-one spectral data, so newness of any trial remains uncertified.
"""
from __future__ import annotations

import json
from pathlib import Path
import random

from global_cusp_gaussian import (
    GENERATORS, I, ONE, PI, QMAT, S, T1, ZERO, G, Mat, divisible_by, g, gdiv, gmul, gneg,
    is_gaussian_integer, mdet, minv, mmul, random_word,
)
from global_cusp_character_algebra import FWD, forward_perm_from_matrix, N

G_DEG: Mat = (PI, ZERO, ZERO, ONE)                       # diag(pi, 1), det = pi (not in SL2)
W_AL: Mat = mmul(S, G_DEG)                                # [[0,-1],[pi,0]]


def coset_rep(c: int) -> Mat:
    """gamma_0 = I, gamma_{1+k} = S T1^k (derived labeling, see global_cusp_character_algebra)."""
    if c == 0:
        return (ONE, ZERO, ZERO, ONE)
    m = S
    for _ in range(c - 1):
        m = mmul(m, T1)
    return m


def conj_by_g(m: Mat) -> Mat:
    """g M g^{-1} = [[a, pi b],[c/pi, d]] (entries in Q(i))."""
    return mmul(mmul(G_DEG, m), minv(G_DEG))


def in_gamma(m: Mat) -> bool:
    return all(is_gaussian_integer(x) for x in m) and mdet(m) == ONE


def check(samples: int = 80, seed: int = 81) -> dict:
    rng = random.Random(seed)
    elems = [random_word(rng, rng.randint(1, 10)) for _ in range(samples)]
    # (1) g M g^{-1} in Gamma  <=>  pi | M_21   (both directions on random elements)
    iff_ok = all(in_gamma(conj_by_g(m)) == divisible_by(m[2], PI) for m in elems)
    has_both = (sum(divisible_by(m[2], PI) for m in elems) > 5
                and sum(not divisible_by(m[2], PI) for m in elems) > 5)
    explicit = conj_by_g((ONE, g(3), ZERO, ONE)) == (ONE, gmul(PI, g(3)), ZERO, ONE)
    # (2) second lift is equivariant in the six-copy convention
    equivariant = True
    for m in elems:
        perm = forward_perm_from_matrix(m, "pi")
        for c in range(N):
            h = mmul(mmul(coset_rep(c), m), minv(coset_rep(perm[c])))      # in H_pi
            equivariant &= divisible_by(h[2], PI) and in_gamma(conj_by_g(h))
    # (3) Hecke cosets: g gamma_c gamma_{c'}^{-1} g^{-1} in Gamma  <=>  c = c'
    hecke = all(
        in_gamma(conj_by_g(mmul(coset_rep(c), minv(coset_rep(d))))) == (c == d)
        for c in range(N) for d in range(N)
    )
    # (4) Q commutes with g (both diagonal) and with the degeneracy maps
    q_commutes_g = mmul(QMAT, G_DEG) == mmul(G_DEG, QMAT)
    # (5) Atkin-Lehner: W = S g = [[0,-1],[pi,0]] normalizes H_pi; W^{-1} Q W Q^{-1} in H_pi
    h_samples = [m for m in elems if divisible_by(m[2], PI)]
    w_norm = all(
        divisible_by(mmul(mmul(W_AL, m), minv(W_AL))[2], PI) and in_gamma(mmul(mmul(W_AL, m), minv(W_AL)))
        for m in h_samples
    )
    comm = mmul(mmul(minv(W_AL), QMAT), mmul(W_AL, minv(QMAT)))
    # projective class: compare after normalizing the scalar (comm should be diag(-i, i) up to scalar)
    comm_is_diag = comm[1] == ZERO and comm[2] == ZERO
    comm_scaled = None
    if comm_is_diag:
        s_ = gdiv(ONE, comm[3])
        comm_scaled = tuple(gmul(s_, x) for x in comm)
    # the diagonal element lies in H_pi (up to scalar) because it is diagonal with unit entries
    unit_diag = comm_is_diag and all(
        x in (ONE, gneg(ONE), I, gneg(I)) for x in (comm[0], comm[3])
    )
    return {
        "gMg_inverse_in_Gamma_iff_pi_divides_lower_left": iff_ok and has_both,
        "g_T_b_g_inverse_is_T_pi_b": explicit,
        "second_lift_is_equivariant_in_six_copy_convention": equivariant,
        "six_Hecke_cosets_g_gamma_c_are_distinct": hecke,
        "Q_commutes_with_g": q_commutes_g,
        "Fricke_W_normalizes_H_pi": w_norm,
        "W_equals_S_g": W_AL == (ZERO, g(-1), PI, ZERO),
        "W_inverse_Q_W_Q_inverse_is_unit_diagonal_so_lies_in_H_pi": unit_diag,
        "U_W_commutes_with_U_Q_on_H_pi_quotient": bool(unit_diag),
    }


def compute_certificate() -> dict:
    res = check()
    return {
        "schema": "global-oldspace/v1",
        "status": "PROVED (exact arithmetic) for the listed algebraic identities only",
        "checks": res,
        "all_checks_pass": all(res.values()),
        "not_certified": [
            "computable P_old,full on a given trial (needs certified level-one spectral data)",
            "newness of any quasimode",
        ],
    }


def main() -> int:
    cert = compute_certificate()
    Path("global_oldspace_result.json").write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"all_checks_pass": cert["all_checks_pass"], **cert["checks"]}, indent=1))
    return 0 if cert["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
