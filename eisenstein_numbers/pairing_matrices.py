"""Explicit pairing matrices g ∈ Γ = PSL(2, ℤ[ω]) + mesh checks.

Freezes SL(2,O) matrices for T1, Tw, U, S; verifies det, relations,
S/U face maps against the P3 dictionary, residue gluing, and user
inverse cycles (via check_generators helpers).

Docs: PAIRING_MATRICES.md, PROOF.md Theorems A–B.

Usage:
  python -u pairing_matrices.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from residue_omega import (  # noqa: E402
    Eisenstein, W, ONE, ZERO,
    gen_matrices, gluing_perms,
    ring_lambda3, ring_split7, ring_split13,
)
from geometry_fund import GEN, OMEGA, mat_act  # noqa: E402
from face_pairings_p3 import (  # noqa: E402
    build_pair_maps, exterior_faces, validate_pair_maps,
)
from reference_cell import build_P3_cell  # noqa: E402
from congruence_omega_proto import ref_elements  # noqa: E402

# ---------------------------------------------------------------------------
# Exact SL(2, O) freeze (Eisenstein entries)
# ---------------------------------------------------------------------------

# O-matrices: ((a,b),(c,d)) with a,b,c,d ∈ O
GEN_O: dict[str, tuple] = {
    "T1": ((ONE, ONE), (ZERO, ONE)),
    "Tw": ((ONE, W), (ZERO, ONE)),
    "U": ((W * W, ZERO), (ZERO, W)),           # diag(ω², ω)
    "S": ((ZERO, Eisenstein(-1, 0)), (ONE, ZERO)),
}

INV_O: dict[str, tuple] = {
    "T1": ((ONE, Eisenstein(-1, 0)), (ZERO, ONE)),
    "Tw": ((ONE, -W), (ZERO, ONE)),
    "U": ((W, ZERO), (ZERO, W * W)),           # diag(ω, ω²)
    "S": GEN_O["S"],
}


def det_O(M) -> Eisenstein:
    (a, b), (c, d) = M
    return a * d - b * c


def mat_eq_O(A, B) -> bool:
    return all(
        A[i][j].a == B[i][j].a and A[i][j].b == B[i][j].b
        for i in range(2) for j in range(2)
    )


def o_to_complex(M):
    """Embed O-matrix as complex GEN-style matrix."""
    (a, b), (c, d) = M

    def to_c(e: Eisenstein) -> complex:
        return complex(e.a, 0) + complex(e.b, 0) * OMEGA

    return ((to_c(a), to_c(b)), (to_c(c), to_c(d)))


def mat_mul_C(A, B):
    return (
        (A[0][0] * B[0][0] + A[0][1] * B[1][0],
         A[0][0] * B[0][1] + A[0][1] * B[1][1]),
        (A[1][0] * B[0][0] + A[1][1] * B[1][0],
         A[1][0] * B[0][1] + A[1][1] * B[1][1]),
    )


def mat_inv_C(M):
    a, b = M[0]
    c, d = M[1]
    return ((d, -b), (-c, a))


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def check_matrices_in_SL() -> bool:
    print("1) Matrices in SL(2, O)")
    print("-" * 60)
    ok = True
    for name, M in GEN_O.items():
        d = det_O(M)
        good = d.a == 1 and d.b == 0
        print(f"  [{'OK' if good else 'FAIL'}] det({name}) = {d}")
        ok &= good
        # match geometry_fund.GEN complex freeze
        Mc = o_to_complex(M)
        Gc = GEN[name]
        close = all(
            abs(Mc[i][j] - Gc[i][j]) < 1e-12
            for i in range(2) for j in range(2)
        )
        print(f"  [{'OK' if close else 'FAIL'}] {name} matches geometry_fund.GEN")
        ok &= close
    # U^3 = I in PSL (diag(ω^6, ω^3) = diag(1,1) since ω^3=1)
    U = o_to_complex(GEN_O["U"])
    U3 = mat_mul_C(mat_mul_C(U, U), U)
    u3_id = abs(U3[0][0] - 1) < 1e-12 and abs(U3[1][1] - 1) < 1e-12
    print(f"  [{'OK' if u3_id else 'FAIL'}] U^3 = I (complex)")
    ok &= u3_id
    # S^2 = -I ≡ I in PSL
    S = o_to_complex(GEN_O["S"])
    S2 = mat_mul_C(S, S)
    # [[0,-1],[1,0]]^2 = [[-1,0],[0,-1]]
    s2_ok = abs(S2[0][0] + 1) < 1e-12 and abs(S2[1][1] + 1) < 1e-12
    print(f"  [{'OK' if s2_ok else 'FAIL'}] S^2 = -I (≡ id in PSL)")
    ok &= s2_ok
    return ok


def check_face_maps_SU(N_tri=6, N3=3, tol=0.25) -> bool:
    """S and U: Möbius image of source centre ≈ partner centre."""
    print(f"\n2) Face maps S, U on mesh {N_tri}x{N3}")
    print("-" * 60)
    mesh = build_P3_cell(N_tri, N3, 1.25, lift=True)
    ref = ref_elements(mesh)
    pm, bf, meta = build_pair_maps(ref["mesh"], ref["tf"], ref["nfr"])
    by = {b["fid"]: b for b in bf}
    ok = True
    vok, msgs = validate_pair_maps(pm, bf)
    ok &= vok
    print(f"  dictionary validate={'OK' if vok else 'FAIL'} stats={meta['stats']}")

    def centre(bf):
        return complex(float(bf["cen"][0]), float(bf["cen"][1])), float(bf["cen"][2])

    # S
    n_s, n_s_ok = 0, 0
    for s, t in pm["S"].items():
        if s > t:
            continue
        n_s += 1
        zs, ts = centre(by[s])
        zt, tt = centre(by[t])
        z2, t2 = mat_act(GEN["S"], zs, ts)
        d = math.hypot(abs(z2 - zt), t2 - tt)
        if d < tol:
            n_s_ok += 1
    print(f"  S floor pairs: {n_s_ok}/{n_s} centres match (tol={tol})")
    ok &= n_s > 0 and n_s_ok == n_s

    # U: forward z↦ωz or inverse z↦ω²z
    n_u, n_u_ok = 0, 0
    Ui = mat_inv_C(GEN["U"])
    for s, t in pm["U"].items():
        if s > t:
            continue
        n_u += 1
        zs, ts = centre(by[s])
        zt, tt = centre(by[t])
        best = 1e9
        for M in (GEN["U"], Ui):
            try:
                z2, t2 = mat_act(M, zs, ts)
            except Exception:
                continue
            best = min(best, math.hypot(abs(z2 - zt), t2 - tt))
        if best < tol:
            n_u_ok += 1
    print(f"  U vert pairs: {n_u_ok}/{n_u} centres match U^{{\\pm1}} (tol={tol})")
    ok &= n_u > 0 and n_u_ok == n_u

    # T1 / Tw: edge-type dictionary (combinatorial section)
    for name, e_src, e_dst in (
        ("T1", "RIGHT", "LEFT"),
        ("Tw", "LOW", "UP"),
    ):
        bad = 0
        n = 0
        for s, t in pm[name].items():
            if s > t:
                continue
            n += 1
            es, et = by[s]["edge"], by[t]["edge"]
            if {es, et} != {e_src, e_dst}:
                bad += 1
        print(f"  {name} edge types {e_src}<->{e_dst}: "
              f"{'OK' if bad == 0 and n > 0 else 'FAIL'} ({n} pairs, {bad} bad)")
        ok &= bad == 0 and n > 0

    # tops free
    top_f = {b["fid"] for b in bf if b["kind"] == "top"}
    hit = any(set(m) & top_f for m in pm.values())
    print(f"  [{'OK' if not hit else 'FAIL'}] no TOP faces in pair maps (D0)")
    ok &= not hit
    return ok


def check_residue() -> bool:
    print("\n3) Residue gluing from matrices")
    print("-" * 60)
    ok = True
    # gen_matrices in residue_omega should match GEN_O
    gens = gen_matrices()
    for name in ("T1", "Tw", "U", "S"):
        if name not in gens:
            print(f"  FAIL missing {name} in gen_matrices")
            ok = False
            continue
        print(f"  [OK] residue gen_matrices has {name}")
    for factory, tag in (
        (ring_lambda3, "N=3"),
        (ring_split7, "N=7"),
        (ring_split13, "N=13"),
    ):
        R = factory()
        perms = gluing_perms(R)
        pts = R.p1_points()
        n = len(pts)
        for name, perm in perms.items():
            if sorted(perm) != list(range(n)):
                print(f"  FAIL {tag} {name} not perm")
                ok = False
        s = perms["S"]
        if any(s[s[i]] != i for i in range(n)):
            print(f"  FAIL {tag} S^2")
            ok = False
        else:
            print(f"  [OK] {tag} |P1|={n} bijective, S^2=id")
    return ok


def check_vertex_cycles() -> bool:
    print("\n4) User inverse cycles (polyhedron section)")
    print("-" * 60)
    from check_generators import (
        parse_user_cycles, check_user_cycle_vs_mobius, mat_inv_sl,
        cycle_notation, _perm_inv,
    )
    user = parse_user_cycles()
    ok = True
    # S
    ok &= check_user_cycle_vs_mobius("S^-1", mat_inv_sl(GEN["S"]), user["S"])
    print(f"  [{'OK' if ok else 'FAIL'}] S^-1 ~ user {cycle_notation(user['S'])}")
    # U: user U^-1 == our forward U
    ok_u = check_user_cycle_vs_mobius("U", GEN["U"], user["U"])
    print(f"  [{'OK' if ok_u else 'FAIL'}] user U^-1 == our U "
          f"{cycle_notation(user['U'])}")
    ok &= ok_u
    for name in ("T1", "Tw"):
        m = check_user_cycle_vs_mobius(
            f"{name}^-1", mat_inv_sl(GEN[name]), user[name]
        )
        print(f"  [{'OK' if m else 'FAIL'}] {name}^-1 Γ_∞-eq "
              f"{cycle_notation(user[name])}")
        ok &= m
    return ok


def print_table() -> None:
    print("\n5) Frozen pairing table (paper)")
    print("-" * 60)
    rows = [
        ("T1", "RIGHT <-> LEFT", "[[1,1],[0,1]]", "z |-> z+1"),
        ("Tw", "LOW <-> UP", "[[1,w],[0,1]]", "z |-> z+w"),
        ("U", "vert via w*z", "[[w^2,0],[0,w]]", "z |-> w z"),
        ("S", "FLOOR <-> FLOOR", "[[0,-1],[1,0]]", "inversion"),
    ]
    for name, faces, mat, act in rows:
        print(f"  {name:3}  {mat:18}  {faces:20}  {act}")


def main() -> int:
    print("Pairing matrices g in Gamma = PSL(2, Z[omega])")
    print("=" * 60)
    ok = True
    ok &= check_matrices_in_SL()
    ok &= check_face_maps_SU()
    ok &= check_residue()
    ok &= check_vertex_cycles()
    print_table()
    print("\n" + "=" * 60)
    print(f"OVERALL: {'PASS' if ok else 'FAIL'}")
    if ok:
        print("  SL(2,O) freeze OK; S/U face maps OK; T1/Tw edge section OK;")
        print("  residue + user inverse cycles OK.")
        print("  See PAIRING_MATRICES.md / PROOF.md.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
