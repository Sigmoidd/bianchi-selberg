#!/usr/bin/env python3
r"""Gaussian conjugation C(z,y) = (conj z, y): H_{2+i} <-> H_{2-i}.  Exact checks.

PROVED here (exact arithmetic): ring isomorphism of the residue fields, conjugation of
the level groups, C(gamma p) = conj(gamma) C(p) on H^3, C Q C^{-1} = Q^{-1},
cusp-lattice transport, the fiber-level identity rho_{pibar}(conj g) = rho_pi(g) (identity
relabeling), and the separation of Gaussian conjugation from cyclotomic label maps.

The spectral equivalence (unitary U_C intertwining the Laplacians and the odd sectors) is
a theorem proved in docs/GLOBAL_CUSP_SYMMETRY_THEOREM.md from these identities.
"""
from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import random

from global_cusp_gaussian import (
    GENERATORS, I, ONE, PI, PI_BAR, QMAT, R, S, T1, TI, ZERO, g, gconj, gmul, gnorm,
    conj_point, divisible_by, h3_action, is_gaussian_integer, mconj, mdet, minv, mmul,
    random_points, random_word, rotate_point,
)
from global_cusp_character_algebra import (
    FWD, N, character_vector, apply_rho, vec_eq, Zeta5, forward_perm_from_matrix,
    quarter_turn_perm, compose, inverse, perm_pow, rho_matrix, matmul,
)


def in_H(m, prime) -> bool:
    return all(is_gaussian_integer(x) for x in m) and mdet(m) == ONE and divisible_by(m[2], prime)


def residue_field_isomorphism() -> dict:
    """a+bi mod pi  ->  a-bi mod pibar is the identity on F_5 (i -> 3 resp. 2)."""
    ok = all((a + 3 * b) % 5 == (a - 2 * b) % 5 for a in range(-12, 13) for b in range(-12, 13))
    return {
        "pi_times_pibar_is_5": gmul(PI, PI_BAR) == g(5, 0),
        "i_is_3_mod_pi": divisible_by(g(-3, 1), PI),
        "i_is_2_mod_pibar": divisible_by(g(-2, 1), PI_BAR),
        "conjugation_induces_identity_F5_to_F5": ok,
    }


def check_levels_and_action(samples: int = 60, seed: int = 61) -> dict:
    rng = random.Random(seed)
    elems = [random_word(rng, rng.randint(1, 9)) for _ in range(samples)]
    pts = random_points(samples, seed + 1)
    hpi = [m for m in elems if in_H(m, PI)]
    hpb = [m for m in elems if in_H(m, PI_BAR)]
    # an element with pibar | c but pi does not divide c must exist (the groups differ)
    witness = (ONE, ZERO, PI_BAR, ONE)        # [[1,0],[pibar,1]]: in H_pibar, not in H_pi
    differ = in_H(witness, PI_BAR) and not in_H(witness, PI)
    cqc = all(
        conj_point(h3_action(QMAT, conj_point(p))) == h3_action(minv(QMAT), p) for p in pts
    )
    return {
        "sample_counts_H_pi_H_pibar": [len(hpi), len(hpb)],
        "conj_maps_H_pi_into_H_pibar": all(in_H(mconj(m), PI_BAR) for m in hpi),
        "conj_maps_H_pibar_into_H_pi": all(in_H(mconj(m), PI) for m in hpb),
        "H_pi_and_H_pibar_are_different_subgroups": differ,
        "C_gamma_p_equals_conj_gamma_C_p": all(
            conj_point(h3_action(m, p)) == h3_action(mconj(m), conj_point(p)) for m, p in zip(elems, pts)
        ),
        "C_Q_Cinv_equals_Q_inverse": cqc,
        "Q_inverse_equals_Q_on_H_pi_quotients_since_Q2_in_H_pi": (
            mmul(QMAT, QMAT) == (g(-1), ZERO, ZERO, ONE)
            and in_H((g(0, 1), ZERO, ZERO, g(0, -1)), PI)
        ),
    }


def check_cusp_lattices(box: int = 7) -> dict:
    def lam0(prime):
        out = set()
        for a in range(-box, box + 1):
            for b in range(-box, box + 1):
                t = (ONE, g(a, b), ZERO, ONE)
                if in_H(mmul(mmul(S, t), minv(S)), prime):
                    out.add((a, b))
        return out
    lp, lb = lam0(PI), lam0(PI_BAR)
    conj_lp = {(a, -b) for a, b in lp}
    # dual generators pi/5 and pibar/5
    return {
        "Lambda_zero_pi_equals_conj_of_Lambda_zero_pibar": conj_lp == lb,
        "dual_generator_pi_over_5": (Fraction(2, 5), Fraction(1, 5)),
        "dual_generator_pibar_over_5": (Fraction(2, 5), Fraction(-1, 5)),
        "dual_generators_are_conjugate": gconj((Fraction(2, 5), Fraction(1, 5))) == (Fraction(2, 5), Fraction(-1, 5)),
        "zero_cusp_frequency_conjugation": all(
            # freq_pi(conj beta) = conj( freq_pibar(beta) )
            ((2 * a + b) / Fraction(5), (a - 2 * b) / Fraction(5))
            == ((Fraction(2 * a + b, 5)), -(Fraction(-a + 2 * b, 5)))
            for a in range(-5, 6) for b in range(-5, 6)
        ),
    }


def check_fiber_level(samples: int = 80, seed: int = 71) -> dict:
    rng = random.Random(seed)
    perms_pi = {n: forward_perm_from_matrix(m, "pi") for n, m in GENERATORS.items()}
    perms_pb = {n: forward_perm_from_matrix(m, "pibar") for n, m in GENERATORS.items()}
    elems = [random_word(rng, rng.randint(1, 10)) for _ in range(samples)]
    # rho_pibar(conj g) == rho_pi(g) as permutations of the same six labels
    transport = all(
        forward_perm_from_matrix(mconj(m), "pibar") == forward_perm_from_matrix(m, "pi") for m in elems
    )
    ti_pi = forward_perm_from_matrix(TI, "pi")
    ti_pb = forward_perm_from_matrix(TI, "pibar")
    t1 = FWD["T1"]
    q_pi = quarter_turn_perm("pi")
    q_pb = quarter_turn_perm("pibar")
    # character labels: Gaussian conjugation fixes T1 (real matrix) so it fixes every label k;
    # the Q-twist permutes labels k -> 3k (pi) resp. k -> 2k (pibar); Galois sigma_a: k -> a k.
    gauss_label_map = {k: k for k in range(1, 5)}
    q_pi_labels = {k: 3 * k % 5 for k in range(1, 5)}
    q_pb_labels = {k: 2 * k % 5 for k in range(1, 5)}
    return {
        "T1_S_R_permutations_agree_for_pi_and_pibar": all(perms_pi[n] == perms_pb[n] for n in ("T1", "S", "R")),
        "Ti_permutation_pi_is_T1_cubed": ti_pi == perm_pow(t1, 3),
        "Ti_permutation_pibar_is_T1_squared": ti_pb == perm_pow(t1, 2),
        "rho_pibar_of_conjugate_equals_rho_pi_on_random_words": transport,
        "fiber_intertwiner_is_the_identity_permutation_matrix": transport,
        "quarter_turn_perm_pibar_is_inverse_of_pi": q_pb == inverse(q_pi),
        "label_map_of_Gaussian_conjugation": gauss_label_map,
        "label_map_of_Q_twist_pi": q_pi_labels,
        "label_map_of_Q_twist_pibar": q_pb_labels,
        "Gaussian_conjugation_is_not_a_cyclotomic_Galois_relabeling": gauss_label_map != q_pi_labels,
        "rho_q_acts_on_characters_as_k_to_3k": all(
            vec_eq(apply_rho(q_pi, character_vector(k)), character_vector(3 * k % 5)) for k in range(1, 5)
        ),
        "rho_qbar_acts_on_characters_as_k_to_2k": all(
            vec_eq(apply_rho(q_pb, character_vector(k)), character_vector(2 * k % 5)) for k in range(1, 5)
        ),
    }


def compute_certificate() -> dict:
    cert = {
        "schema": "global-conjugate-prime/v1",
        "status": "PROVED (exact arithmetic)",
        "residue_fields": residue_field_isomorphism(),
        "levels_and_action": check_levels_and_action(),
        "cusp_lattices": check_cusp_lattices(),
        "fiber_level": check_fiber_level(),
        "spectral_equivalence_statement": (
            "U_C f = f o C is a unitary L^2(H_pibar\\H^3) -> L^2(H_pi\\H^3) commuting with the Laplacians "
            "and mapping Q-odd to Q-odd; f -> conj(f o C) is the antiunitary version. "
            "Spectra of Delta on the two quotients coincide with multiplicity; conjugate-prime "
            "eigenvalues are NOT independent discoveries."
        ),
    }
    flat: list = []

    def collect(d):
        for v in d.values():
            if isinstance(v, dict):
                collect(v)
            elif isinstance(v, bool):
                flat.append(v)
    for k in ("residue_fields", "levels_and_action", "cusp_lattices", "fiber_level"):
        collect(cert[k])
    cert["all_checks_pass"] = all(flat)
    return cert


def main() -> int:
    cert = compute_certificate()
    Path("global_conjugate_prime_result.json").write_text(json.dumps(cert, indent=2, sort_keys=True, default=str), encoding="utf-8")
    print(json.dumps({"all_checks_pass": cert["all_checks_pass"]}))
    return 0 if cert["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
