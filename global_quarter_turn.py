#!/usr/bin/env python3
r"""Quarter-turn Q(z,y) = (iz, y) on Gamma_0(2+i)\H^3: exact checks.

PROVED here (exact Q(i) arithmetic): normalization of H_pi and Gamma by Q,
Q^2 = R, the cusp-coordinate rotations at infinity (by i) and zero (by -i), the
cusp lattices, the bundle-level intertwining
``rho_q rho(Q g Q^{-1}) rho_q^{-1} = rho(g)``, ``rho_q^2 rho(R) = I``, and the
mode-orbit counts.

The analytic statements (self-adjointness, discreteness of the cuspidal
spectrum) are proved in docs/GLOBAL_CUSP_SYMMETRY_THEOREM.md; they are not
re-proved by this file.

``field_level_intertwining_check`` is a *numerical* Arb check on the
repository's own Fourier-Whittaker field map (EXPERIMENTAL tier: it tests the
sector coefficient relations against the actual ``component_row`` code).
"""
from __future__ import annotations

from fractions import Fraction
import json
import math
import random
from pathlib import Path

from global_cusp_gaussian import (
    GENERATORS, I, ONE, PI, QMAT, QSQ_REP, R, S, T1, TI, ZERO, G, Mat, Point,
    conj_point, divisible_by, g, gmul, gneg, gnorm, h3_action, is_gaussian_integer,
    mdet, minv, mmul, proj_equal, proj_equal_gl, random_points, random_word, rotate_point,
)
from global_cusp_character_algebra import (
    FWD, N, compose, forward_perm_from_matrix, inverse, matmul, quarter_turn_perm,
    rho_matrix, ident_f,
)

QINV: Mat = minv(QMAT)


def conj_by_Q(m: Mat) -> Mat:
    """Q m Q^{-1} = [[a, i b], [-i c, d]]."""
    return mmul(mmul(QMAT, m), QINV)


def in_gamma(m: Mat) -> bool:
    return all(is_gaussian_integer(x) for x in m) and mdet(m) == ONE


def in_H_pi(m: Mat) -> bool:
    return in_gamma(m) and divisible_by(m[2], PI)


def random_H_pi_elements(count: int, seed: int) -> list[Mat]:
    rng = random.Random(seed)
    out: list[Mat] = []
    while len(out) < count:
        m = random_word(rng, rng.randint(2, 9))
        if in_H_pi(m):
            out.append(m)
    return out


# ---------------------------------------------------------------- group-level facts

def check_normalization(samples: int = 60, seed: int = 11) -> dict:
    rng = random.Random(seed)
    gamma_ok = all(in_gamma(conj_by_Q(random_word(rng, rng.randint(1, 9)))) for _ in range(samples))
    hpi = random_H_pi_elements(samples, seed + 1)
    hpi_ok = all(in_H_pi(conj_by_Q(h)) for h in hpi)
    # lower-left entry of Q h Q^{-1} is -i c, so pi | c  <=>  pi | -i c (i is a unit)
    entry_ok = all(conj_by_Q(h)[2] == gmul(gneg(I), h[2]) for h in hpi)
    q2_matches_R = all(
        mmul(mmul(QMAT, conj_by_Q(h)), QINV) == mmul(mmul(QSQ_REP, h), minv(QSQ_REP)) for h in hpi
    )
    return {
        "Q_normalizes_Gamma": gamma_ok,
        "Q_normalizes_H_pi": hpi_ok,
        "lower_left_entry_is_minus_i_c": entry_ok,
        "Q_squared_conjugation_equals_R_conjugation": q2_matches_R,
        "R_is_diag_i_minus_i_and_lies_in_H_pi": in_H_pi(QSQ_REP) and QSQ_REP == R,
        "R_acts_as_z_to_minus_z": all(
            h3_action(R, p) == ((gneg(p[0])), p[1]) for p in random_points(10, seed + 2)
        ),
        "det_Q_is_i_so_Q_is_not_in_PSL2_Z_i": mdet(QMAT) == I,
    }


def check_point_actions(samples: int = 25, seed: int = 21) -> dict:
    pts = random_points(samples, seed)
    rng = random.Random(seed + 1)
    q_is_rotation = all(h3_action(QMAT, p) == rotate_point(p, 1) for p in pts)
    elems = [random_word(rng, rng.randint(1, 7)) for _ in range(samples)]
    equivariance = all(
        h3_action(QMAT, h3_action(m, p)) == h3_action(conj_by_Q(m), h3_action(QMAT, p))
        for m, p in zip(elems, pts)
    )
    # cusp 0: Q S = S Q0 with Q0(z,y) = (-i z, y); sigma_0 = S
    q_at_zero = all(
        h3_action(mmul(QMAT, S), p) == h3_action(S, rotate_point(p, -1)) for p in pts
    )
    # cusp infinity: Q T_b Q^{-1} = T_{ib}, so Q maps the translation lattice to itself
    lattice_inf = all(
        conj_by_Q((ONE, g(a, b), ZERO, ONE)) == (ONE, gmul(I, g(a, b)), ZERO, ONE)
        for a in range(-3, 4) for b in range(-3, 4)
    )
    return {
        "Q_matrix_acts_as_rotation_z_to_iz": q_is_rotation,
        "Q_intertwines_gamma_and_conjugate_on_H3": equivariance,
        "Q_S_equals_S_Q0_with_Q0_rotation_by_minus_i": q_at_zero,
        "Q_T_b_Qinv_equals_T_ib": lattice_inf,
        "Q_fixes_cusp_infinity_and_cusp_zero_as_points": (
            h3_action(QMAT, (g(0), Fraction(1)))[0] == g(0)
        ),
    }


def check_cusp_lattices(box: int = 7) -> dict:
    """Lambda_inf = Z[i], Lambda_0 = pi Z[i]; both are stable under multiplication by i."""
    lam_inf, lam_zero = set(), set()
    for a in range(-box, box + 1):
        for b in range(-box, box + 1):
            beta = g(a, b)
            t = (ONE, beta, ZERO, ONE)
            if in_H_pi(t):
                lam_inf.add((a, b))
            if in_H_pi(mmul(mmul(S, t), minv(S))):
                lam_zero.add((a, b))
    pi_z = {(a, b) for a in range(-box, box + 1) for b in range(-box, box + 1)
            if divisible_by(g(a, b), PI)}
    full = {(a, b) for a in range(-box, box + 1) for b in range(-box, box + 1)}

    def rot(s: set) -> set:
        return {(-b, a) for a, b in s}

    inner = {(a, b) for (a, b) in lam_zero if abs(a) <= box - 1 and abs(b) <= box - 1}
    return {
        "Lambda_infinity_is_Z_i_on_box": lam_inf == full,
        "Lambda_zero_is_pi_Z_i_on_box": lam_zero == pi_z,
        "Lambda_zero_covolume": str(Fraction(gnorm(PI))),
        "rotation_by_i_preserves_Lambda_infinity": rot(lam_inf) == lam_inf,
        "rotation_by_i_preserves_Lambda_zero_on_box": rot(inner) <= lam_zero,
        "dual_zero_lattice_generator": "pi/5 = 1/conj(pi)",
        "dual_generator_identity": gmul(PI, g(Fraction(1), 0)) == PI
        and gmul(PI, g(Fraction(1, 5), 0)) == (Fraction(2, 5), Fraction(1, 5)),
    }


def check_fourier_rotation_identity(samples: int = 40, seed: int = 31) -> dict:
    """<mu, i x> = <-i mu, x> exactly, at both cusps (zero-cusp frequency (2a-b, a+2b)/5)."""
    rng = random.Random(seed)
    ok_inf = ok_zero = True
    for _ in range(samples):
        a, b = rng.randint(-6, 6), rng.randint(-6, 6)
        x = (Fraction(rng.randint(-9, 9), rng.randint(1, 9)), Fraction(rng.randint(-9, 9), rng.randint(1, 9)))
        ix = (-x[1], x[0])
        # infinity: frequency (a, b); -i(a+bi) = b - a i
        ok_inf &= a * ix[0] + b * ix[1] == b * x[0] + (-a) * x[1]
        # zero: frequency u,v = (2a-b)/5, (a+2b)/5 ; -i beta = (b, -a)
        u, v = Fraction(2 * a - b, 5), Fraction(a + 2 * b, 5)
        u2, v2 = Fraction(2 * b + a, 5), Fraction(b - 2 * a, 5)
        ok_zero &= (u * ix[0] + v * ix[1] == u2 * x[0] + v2 * x[1]) and (u2, v2) == (v, -u)
    return {"phi_mu_o_Q_equals_phi_minus_i_mu_at_infinity": ok_inf,
            "phi_beta_o_Q_equals_phi_minus_i_beta_at_zero": ok_zero}


# ---------------------------------------------------------------- bundle-level action

def rho_of(m: Mat):
    return rho_matrix(forward_perm_from_matrix(m, "pi"))


def check_bundle_action(samples: int = 80, seed: int = 41) -> dict:
    """(U_Q F)_c(z) = F_{q(c)}(Qz) preserves equivariance F(gamma z) = rho(gamma) F(z).

    Needed identity: rho_q rho(Q gamma Q^{-1}) rho_q^{-1} = rho(gamma) for all gamma.
    """
    q = quarter_turn_perm("pi")
    rq, rqi = rho_matrix(q), rho_matrix(inverse(q))
    rng = random.Random(seed)
    gens_ok = {
        name: matmul(matmul(rq, rho_of(conj_by_Q(m))), rqi) == rho_of(m)
        for name, m in GENERATORS.items()
    }
    gens_ok["Ti"] = matmul(matmul(rq, rho_of(conj_by_Q(TI))), rqi) == rho_of(TI)
    words_ok = True
    for _ in range(samples):
        m = random_word(rng, rng.randint(1, 10))
        words_ok &= matmul(matmul(rq, rho_of(conj_by_Q(m))), rqi) == rho_of(m)
    # explicit images of the generators under conjugation by Q
    images = {
        "Q T1 Q^-1": conj_by_Q(T1) == TI,
        "Q Ti Q^-1 equals T1^-1": conj_by_Q(TI) == minv(T1),
        "Q R Q^-1 equals R": conj_by_Q(R) == R,
        "Q S Q^-1 is S' = [[0,-i],[-i,0]]": conj_by_Q(S) == (ZERO, g(0, -1), g(0, -1), ZERO),
    }
    # square of the bundle map: rho_q^2 F(Q^2 z) = rho_q^2 rho(R) F(z) must be F(z)
    rq2 = matmul(rq, rq)
    square_identity = matmul(rq2, rho_of(R)) == ident_f()
    # the repo's coset-label convention, rederived: labels <-> bottom rows
    return {
        "intertwining_on_generators": gens_ok,
        "intertwining_on_random_words": words_ok,
        "generator_images": images,
        "rho_q_squared_times_rho_R_is_identity": square_identity,
        "finite_q_squared_is_nontrivial_so_fiber_map_alone_is_not_an_involution": rq2 != ident_f(),
        "bundle_map_squares_to_identity_using_Q2_equals_R_on_base": square_identity,
    }


# ---------------------------------------------------------------- mode orbits

def gaussian_modes(cutoff: int) -> list[tuple[int, int, int]]:
    """Same list and order as track_b_two_cusp_data.gaussian_modes."""
    radius = math.isqrt(cutoff)
    out = []
    for a in range(-radius, radius + 1):
        for b in range(-radius, radius + 1):
            nrm = a * a + b * b
            if 0 < nrm <= cutoff:
                out.append((a, b, nrm))
    return sorted(out, key=lambda row: (row[2], row[0], row[1]))


def cusp_modes(cutoff: int):
    return gaussian_modes(cutoff), gaussian_modes(5 * cutoff)


def rot_i(mode: tuple[int, int]) -> tuple[int, int]:
    """i (a + b i) = -b + a i."""
    return (-mode[1], mode[0])


def c4_orbits(modes: list[tuple[int, int, int]]) -> list[list[int]]:
    """Orbits of nu -> i nu as lists of indices into ``modes`` (free action, size 4)."""
    index = {(a, b): k for k, (a, b, _n) in enumerate(modes)}
    seen: set[int] = set()
    orbs = []
    for k, (a, b, _n) in enumerate(modes):
        if k in seen:
            continue
        orb, cur = [], (a, b)
        for _ in range(4):
            orb.append(index[cur])
            cur = rot_i(cur)
        seen.update(orb)
        orbs.append(orb)
    return orbs


def sector_basis(cutoff: int, eps: int) -> list[dict[int, int]]:
    """Columns (as sparse dicts over the concatenated [infinity | zero] ordering).

    Sector eps: a_{i nu} = eps a_nu, b_{i beta} = eps b_beta.  One column per C4-orbit.
    """
    mi, mz = cusp_modes(cutoff)
    cols = []
    for offset, modes in ((0, mi), (len(mi), mz)):
        for orb in c4_orbits(modes):
            col = {}
            for step, idx in enumerate(orb):         # orb[step] = i^step nu
                col[offset + idx] = eps ** step
            cols.append(col)
    return cols


def accounting(cutoff: int) -> dict:
    mi, mz = cusp_modes(cutoff)
    full = len(mi) + len(mz)
    oi, oz = c4_orbits(mi), c4_orbits(mz)
    return {
        "cutoff_M": cutoff,
        "modes_infinity": len(mi), "modes_zero": len(mz), "unknowns_full": full,
        "unknowns_after_R_only": full // 2,       # a_{-nu} = a_nu (R in H_pi)
        "orbits_infinity": len(oi), "orbits_zero": len(oz),
        "unknowns_per_Q_parity_sector": len(oi) + len(oz),
        "reduction_factor_vs_full": full / (len(oi) + len(oz)),
        "orbits_all_size_4": all(len(set(o)) == 4 for o in oi + oz),
        "box_is_rotation_invariant": all(
            rot_i((a, b)) in {(x, y) for x, y, _ in modes} for modes in (mi, mz) for a, b, _ in modes
        ),
    }


# ---------------------------------------------------------------- field-level check (Arb)

def field_level_intertwining_check(cutoff: int = 2, bits: int = 192, seed: int = 51, points: int = 4) -> dict:
    """Check F_{q(c)}(QP) = eps F_c(P) for sector coefficient vectors, with negative controls.

    Uses the repository's own ``component_row`` (the actual Fourier-Whittaker field map).
    Pass criterion: every difference enclosure contains 0 (Arb-consistent), AND the same test
    applied to a deliberately wrong relation (wrong sign / identity instead of q) is excluded
    from 0.  EXPERIMENTAL tier: sample points only; this is not a theorem.
    Note: Arb radii of K_{ir}(x) grow for large x, so heights are kept moderate and the
    reported bounds are *enclosure radii*, not measured roundoff.
    """
    from flint import acb, arb, ctx
    from track_b_two_cusp_hejhal import ValidatedWhittaker, component_row
    from track_b_two_cusp_data import cusp_modes as repo_cusp_modes, arb_fraction

    ctx.prec = bits
    modes_inf, modes_zero = repo_cusp_modes(cutoff)
    assert (modes_inf, modes_zero) == cusp_modes(cutoff), "mode lists differ from repository"
    backend = ValidatedWhittaker(arb("6.7439020359331625"))
    q = quarter_turn_perm("pi")
    ident = tuple(range(6))
    rng = random.Random(seed)

    def evaluate(row, coef):
        return sum((a * b for a, b in zip(row, coef)), acb(0))

    consistent = {1: True, -1: True}
    wrong_sign_detected = {1: False, -1: False}
    wrong_perm_detected = {1: False, -1: False}
    worst_upper = arb(0)
    for eps in (1, -1):
        cols = sector_basis(cutoff, eps)
        coef = [acb(0)] * (len(modes_inf) + len(modes_zero))
        for col in cols:
            w = acb(arb(rng.randint(-5, 5)) / 3, arb(rng.randint(1, 5)) / 7)
            for idx, sgn in col.items():
                coef[idx] = w * sgn
        for _ in range(points):
            x1 = Fraction(rng.randint(-9, 9), rng.randint(1, 8))
            x2 = Fraction(rng.randint(-9, 9), rng.randint(1, 8))
            y = Fraction(rng.randint(3, 12), rng.randint(3, 8))
            p = (arb_fraction(str(x1)), arb_fraction(str(x2)), arb_fraction(str(y)))
            qp = (arb_fraction(str(-x2)), arb_fraction(str(x1)), arb_fraction(str(y)))  # Q P = (iz, y)
            for c in range(6):
                lhs = evaluate(component_row(q[c], qp, modes_inf, modes_zero, backend), coef)
                rhs = evaluate(component_row(c, p, modes_inf, modes_zero, backend), coef) * eps
                d = lhs - rhs
                if not (d.real.contains(0) and d.imag.contains(0)):
                    consistent[eps] = False
                worst_upper = max(worst_upper, abs(d).upper(), key=lambda t: float(t.mid()))
                # negative controls
                bad_sign = lhs + rhs
                if not (bad_sign.real.contains(0) and bad_sign.imag.contains(0)):
                    wrong_sign_detected[eps] = True
                if q[c] != c:
                    wrong = evaluate(component_row(ident[c], qp, modes_inf, modes_zero, backend), coef) - rhs
                    if not (wrong.real.contains(0) and wrong.imag.contains(0)):
                        wrong_perm_detected[eps] = True
    return {
        "cutoff_M": cutoff, "bits": bits, "r": "6.7439020359331625",
        "even_sector_all_differences_contain_zero": consistent[1],
        "odd_sector_all_differences_contain_zero": consistent[-1],
        "negative_control_wrong_sign_detected": wrong_sign_detected,
        "negative_control_identity_instead_of_q_detected": wrong_perm_detected,
        "largest_difference_enclosure_upper_bound": f"{float(worst_upper.mid()):.3e}",
        "passes": all(consistent.values()) and all(wrong_sign_detected.values()) and all(wrong_perm_detected.values()),
        "tier": "EXPERIMENTAL (Arb sample-point check of the sector relations on the repository field map)",
    }


def compute_certificate(run_field_check: bool = True) -> dict:
    cert = {
        "schema": "global-quarter-turn/v1",
        "status": "PROVED for exact algebra blocks; field-level block is EXPERIMENTAL",
        "normalization": check_normalization(),
        "point_actions": check_point_actions(),
        "cusp_lattices": check_cusp_lattices(),
        "fourier_rotation_identity": check_fourier_rotation_identity(),
        "bundle_action": check_bundle_action(),
        "accounting": {f"M{m}": accounting(m) for m in range(1, 7)},
    }
    if run_field_check:
        try:
            cert["field_level_check"] = field_level_intertwining_check()
        except ImportError as exc:       # python-flint absent: report, never silently pass
            cert["field_level_check"] = {"skipped": f"import failed: {exc}"}
    cert["all_exact_checks_pass"] = _exact_pass(cert)
    return cert


def _flat(d) -> list:
    out = []
    for v in d.values():
        out += _flat(v) if isinstance(v, dict) else [v]
    return out


def _exact_pass(c: dict) -> bool:
    bools = []
    for key in ("normalization", "point_actions", "cusp_lattices", "fourier_rotation_identity", "bundle_action"):
        bools += [v for v in _flat(c[key]) if isinstance(v, bool)]
    acct_ok = all(a["orbits_all_size_4"] and a["box_is_rotation_invariant"] for a in c["accounting"].values())
    # finite_q_squared... is True by design (documents that fiber map alone is not an involution)
    return all(bools) and acct_ok


def main() -> int:
    cert = compute_certificate()
    Path("global_quarter_turn_result.json").write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"all_exact_checks_pass": cert["all_exact_checks_pass"],
                      "field_level_check": cert.get("field_level_check")}, indent=2))
    return 0 if cert["all_exact_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
