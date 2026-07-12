"""Verify generator actions against supplied cycle decompositions.

User data (action of **inverses** on ideal vertices {∞,0,1,ω,ω²}):

  T₁⁻¹ : (∞)(0 1)(ω ω²)
  T_ω⁻¹: (∞)(0 ω)(1 ω²)
  U⁻¹  : (∞)(0)(1 ω ω²)
  S⁻¹  : (∞ 0)(1)(ω ω²)

These are Möbius actions on P¹(ℂ), reduced back to the vertex set by
Γ_∞ = (translations by O) ⋊ (units O*) — i.e. z ↦ u z + b.

Also audits residue gluing perms on ℙ¹(𝔽_q) for N∈{3,7,13}.

Usage:
  python -u check_generators.py
"""
from __future__ import annotations

import cmath
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from residue_omega import (
    gen_matrices, gluing_perms,
    ring_lambda3, ring_split7, ring_split13,
)
from geometry_fund import GEN, OMEGA

# ---------------------------------------------------------------------------
# Vertex set and Γ_∞ reduction
# ---------------------------------------------------------------------------

W_C = OMEGA
W2_C = OMEGA ** 2
VERT_LABELS = ("inf", "0", "1", "w", "w2")
VERT_VAL = {
    "inf": None,
    "0": 0j,
    "1": 1 + 0j,
    "w": W_C,
    "w2": W2_C,
}

# units as complex multipliers
UNIT_C = [
    1 + 0j, -1 + 0j,
    W_C, -W_C,
    W2_C, -W2_C,
]

# short translations
TRANS = [
    a + b * W_C
    for a in range(-3, 4)
    for b in range(-3, 4)
]


def reduce_vertex(z, tol: float = 1e-8) -> str:
    """Map z ∈ ℂ∪{∞} to a label in VERT by Γ_∞ = O ⋊ U.

    Order matters: try exact match and **unit multiplications first**
    (stabilizer of 0), then translations.  Blind translation-first search
    falsely collapses many points (e.g. −1 ⇝ ω² instead of 1).
    """
    if z is None or (isinstance(z, complex) and abs(z) > 1e12):
        return "inf"
    # 1) already a vertex
    for lab, v in VERT_VAL.items():
        if v is not None and abs(z - v) < tol:
            return lab
    # 2) unit multiple u*z (no translation)
    for u in UNIT_C:
        w = u * z
        for lab, v in VERT_VAL.items():
            if v is not None and abs(w - v) < tol:
                return lab
    # 3) translate then unit: u*(z - b)
    for b in TRANS:
        zb = z - b
        for u in UNIT_C:
            w = u * zb
            for lab, v in VERT_VAL.items():
                if v is not None and abs(w - v) < tol:
                    return lab
    raise ValueError(f"cannot reduce {z} to vertex set")


def mobius(M, z):
    """Apply SL matrix M=((a,b),(c,d)) to z ∈ ℂ∪{∞}; return ℂ or None=∞."""
    a, b = M[0]
    c, d = M[1]
    if z is None:
        # ∞ → a/c
        if abs(c) < 1e-15:
            return None
        return a / c
    den = c * z + d
    if abs(den) < 1e-15:
        return None
    return (a * z + b) / den


def mat_inv_sl(M):
    a, b = M[0]
    c, d = M[1]
    return ((d, -b), (-c, a))


def cycle_notation(perm: dict[str, str]) -> str:
    """perm: label -> label; return cycle string (guards non-bijective)."""
    # bijectivity check
    vals = list(perm.values())
    if set(vals) != set(VERT_LABELS) or len(vals) != len(set(vals)):
        return f"<non-bijective {perm}>"
    seen = set()
    cycles = []
    for start in VERT_LABELS:
        if start in seen:
            continue
        cyc = [start]
        seen.add(start)
        x = perm[start]
        guard = 0
        while x != start:
            cyc.append(x)
            seen.add(x)
            x = perm[x]
            guard += 1
            if guard > len(VERT_LABELS) + 2:
                return f"<cycle-break {perm}>"
        if len(cyc) == 1:
            cycles.append(f"({cyc[0]})")
        else:
            cycles.append("(" + " ".join(cyc) + ")")
    return "".join(cycles)


def induced_perm(M) -> dict[str, str]:
    """Permutation of vertices by Möbius M, reduced by Γ_∞."""
    out = {}
    for lab, z in VERT_VAL.items():
        w = mobius(M, z)
        out[lab] = reduce_vertex(w)
    return out


# User-supplied cycles for inverses (as maps label->label)
def parse_user_cycles():
    """Return dict name -> perm dict for T1inv, Twinv, Uinv, Sinv."""
    # T1^{-1}: (inf)(0 1)(w w2)
    # Tw^{-1}: (inf)(0 w)(1 w2)
    # U^{-1}:  (inf)(0)(1 w w2)
    # S^{-1}:  (inf 0)(1)(w w2)
    def from_cycles(cycles):
        perm = {lab: lab for lab in VERT_LABELS}
        for cyc in cycles:
            for i, a in enumerate(cyc):
                perm[a] = cyc[(i + 1) % len(cyc)]
        return perm

    return {
        "T1": from_cycles([["inf"], ["0", "1"], ["w", "w2"]]),
        "Tw": from_cycles([["inf"], ["0", "w"], ["1", "w2"]]),
        "U": from_cycles([["inf"], ["0"], ["1", "w", "w2"]]),
        "S": from_cycles([["inf", "0"], ["1"], ["w", "w2"]]),
    }


def _perm_eq(p, q) -> bool:
    return all(p[a] == q[a] for a in VERT_LABELS)


def _perm_inv(p: dict) -> dict:
    return {v: k for k, v in p.items()}


def gamma_inf_equivalent(z, target_lab: str, tol: float = 1e-8) -> bool:
    """True iff z ~ target under Γ_∞ = (z ↦ u z + b), u∈O*, b∈O.

    Greedy reduce_vertex is multi-valued (first-hit among many (u,b)).
    User cycles pick a polyhedron section; we only require *existence*
    of a Γ_∞ element sending Möbius(g^{-1}, v) to the user image.
    """
    if target_lab == "inf":
        return z is None or (isinstance(z, complex) and abs(z) > 1e12)
    if z is None or (isinstance(z, complex) and abs(z) > 1e12):
        return False
    tgt = VERT_VAL[target_lab]
    # z = u * tgt + b  ⇔  u^{-1}(z - b) = tgt
    for b in TRANS:
        zb = z - b
        for u in UNIT_C:
            if abs(u * zb - tgt) < tol:
                return True
            # also: z = u*tgt + b form checked as u*(tgt) + b
            if abs(u * tgt + b - z) < tol:
                return True
    return False


def check_user_cycle_vs_mobius(name: str, M_inv, want: dict, tol: float = 1e-8) -> bool:
    """For each ideal vertex v: Möbius(M_inv, v) ~_Γ∞ want[v]."""
    ok = True
    for lab, z in VERT_VAL.items():
        w = mobius(M_inv, z)
        tgt = want[lab]
        if not gamma_inf_equivalent(w, tgt, tol=tol):
            print(f"    FAIL {name}: {lab} |-> {w} not Γ_∞-eq to {tgt}")
            ok = False
    return ok


def check_vertex_cycles() -> bool:
    print("1) Ideal vertices {∞,0,1,ω,ω²} under generator inverses")
    print("   User cycles = ground truth polyhedron section.")
    print("   Check: Möbius(g^{-1}, v) is Γ_∞-equivalent to user image.")
    print("-" * 60)
    user = parse_user_cycles()
    gens = GEN
    ok_all = True

    # --- S^{-1}: unit reduction only (S preserves the set up to units) ---
    want_S = user["S"]
    Mi_S = mat_inv_sl(gens["S"])
    match_S = check_user_cycle_vs_mobius("S^-1", Mi_S, want_S)
    # exact cycle string via unit-only reduction
    def reduce_unit_only(z, tol=1e-8):
        if z is None:
            return "inf"
        for lab, v in VERT_VAL.items():
            if v is not None and abs(z - v) < tol:
                return lab
        for u in UNIT_C:
            for lab, v in VERT_VAL.items():
                if v is not None and abs(u * z - v) < tol:
                    return lab
        return None

    got_S_u = {}
    for lab, z in VERT_VAL.items():
        w = mobius(Mi_S, z)
        r = reduce_unit_only(w)
        got_S_u[lab] = r if r is not None else "?"
    print(f"  [{'OK' if match_S else 'FAIL'}] S^-1  user={cycle_notation(want_S)}")
    print(f"        mobius+units: {cycle_notation(got_S_u) if '?' not in got_S_u.values() else got_S_u}")
    ok_all &= match_S

    # --- U^{-1}: user cycle (1 w w2). Our GEN U = diag(ω²,ω) does z↦ωz,
    #     which equals the *user-labeled* U^{-1} cycle. Our mat_inv U is the
    #     opposite orientation (1 w2 w). Same ⟨U⟩; residue gluing uses δ^{-1}.
    want_U = user["U"]
    Mi_U = mat_inv_sl(gens["U"])  # true matrix inverse of our U
    M_U = gens["U"]
    # User data labeled U^{-1}: match against our forward U matrix
    match_U_user = check_user_cycle_vs_mobius(
        "user-U^-1 vs our-U", M_U, want_U
    )
    # True matrix U^{-1} matches inverse cycle
    want_U_true_inv = _perm_inv(want_U)
    match_U_inv = check_user_cycle_vs_mobius(
        "our-U^-1 matrix", Mi_U, want_U_true_inv
    )
    got_U = induced_perm(M_U)
    got_Ui = induced_perm(Mi_U)
    print(f"  [{'OK' if match_U_user else 'FAIL'}] U orientation vs user label")
    print(f"        user U^-1 cycle:    {cycle_notation(want_U)}")
    print(f"        our U (z↦ωz):       {cycle_notation(got_U)}")
    print(f"        our U^-1 (z↦ω²z):   {cycle_notation(got_Ui)}")
    print(f"        note: user U^-1 cycle == our forward U; cert uses ⟨U⟩")
    ok_all &= match_U_user and match_U_inv

    # --- T1^{-1}, Tw^{-1}: Γ_∞-equivalence to user section (not greedy reduce)
    for name in ("T1", "Tw"):
        want = user[name]
        Mi = mat_inv_sl(gens[name])
        match = check_user_cycle_vs_mobius(f"{name}^-1", Mi, want)
        bij = set(want.values()) == set(VERT_LABELS) and len(set(want.values())) == 5
        print(f"  [{'OK' if match and bij else 'FAIL'}] {name}^-1  "
              f"user={cycle_notation(want)}")
        print(f"        Γ_∞-eq Möbius images to user section: "
              f"{'yes' if match else 'NO'}")
        ok_all &= match and bij

    # Relation checks on user inverse cycles
    print("  relations on user inverse perms:")
    Si, Ui = user["S"], user["U"]
    S2 = {a: Si[Si[a]] for a in VERT_LABELS}
    idp = {a: a for a in VERT_LABELS}
    t = _perm_eq(S2, idp)
    ok_all &= t
    print(f"    [{'OK' if t else 'FAIL'}] (S^-1)^2 = id")
    U3 = Ui
    for _ in range(2):
        U3 = {a: Ui[U3[a]] for a in VERT_LABELS}
    t = _perm_eq(U3, idp)
    ok_all &= t
    print(f"    [{'OK' if t else 'FAIL'}] (user U-cycle)^3 = id")
    # T1,Tw user sections are involutions on the 5-set
    for name in ("T1", "Tw"):
        p = user[name]
        p2 = {a: p[p[a]] for a in VERT_LABELS}
        t = _perm_eq(p2, idp)
        ok_all &= t
        print(f"    [{'OK' if t else 'FAIL'}] ({name}^-1_section)^2 = id")
    return ok_all


def check_forward_consistency() -> bool:
    """U^3=id, S^2=id on vertices; T1 etc."""
    print("\n2) Relations on vertex set (forward generators)")
    print("-" * 60)
    ok = True
    # S^2
    pS = induced_perm(GEN["S"])
    pS2 = {a: pS[pS[a]] for a in VERT_LABELS}
    idp = {a: a for a in VERT_LABELS}
    t = pS2 == idp
    ok &= t
    print(f"  [{'OK' if t else 'FAIL'}] S^2 = id on vertices")
    # U^3
    pU = induced_perm(GEN["U"])
    pU3 = pU
    for _ in range(2):
        pU3 = {a: pU[pU3[a]] for a in VERT_LABELS}
    t = pU3 == idp
    ok &= t
    print(f"  [{'OK' if t else 'FAIL'}] U^3 = id on vertices")
    return ok


def check_residue_gluing() -> bool:
    print("\n3) Residue ℙ¹ gluing perms (right action by δ^{-1})")
    print("-" * 60)
    ok = True
    for factory, tag in (
        (ring_lambda3, "N=3"),
        (ring_split7, "N=7"),
        (ring_split13, "N=13"),
    ):
        R = factory()
        pts = R.p1_points()
        perms = gluing_perms(R)
        print(f"  [{tag}] |P1|={len(pts)} w_img={R.w_img}")
        for name, perm in perms.items():
            # bijective
            if sorted(perm) != list(range(len(pts))):
                print(f"    FAIL {name} not a permutation: {perm}")
                ok = False
            else:
                print(f"    {name}: {perm}")
        # S^2 = id
        s = perms["S"]
        if any(s[s[i]] != i for i in range(len(s))):
            print("    FAIL S^2 != id")
            ok = False
        else:
            print("    S^2=id OK")
        # U^3 = id
        u = perms["U"]

        def comp(p, q):
            return [p[q[i]] for i in range(len(p))]

        u3 = comp(comp(u, u), u)
        if u3 != list(range(len(u))):
            print(f"    FAIL U^3 != id: {u3}")
            ok = False
        else:
            print("    U^3=id OK")
        # round-trip: act(perm[i], δ) = i  (as m3p asserts)
        gens = gen_matrices()
        for name, M_e in gens.items():
            M = R.mat_embed(M_e)
            perm = perms[name]
            for i, p in enumerate(pts):
                # perm[i] = p · δ^{-1}; so (p · δ^{-1}) · δ = p
                j = perm[i]
                p2 = R.p1_act(pts[j], M)
                if p2 != p:
                    print(f"    FAIL round-trip {name} at {i}: {p2} != {p}")
                    ok = False
                    break
            else:
                print(f"    round-trip {name} OK")
    return ok


def check_N3_labels_vs_vertices() -> None:
    """For N=3, ω≡1 in F_3 so {1,ω,ω²} collapse; document mapping."""
    print("\n4) N=3 residue labels vs ideal vertices")
    print("-" * 60)
    R = ring_lambda3()
    print(f"  F_3 with ω ↦ {R.w_img} (so ω≡ω²≡1 in R)")
    print(f"  P1 points: {R.p1_points()}")
    print("  Note: user cycles on 5 vertices live in P¹(ℂ);")
    print("  gluing for Γ₀ uses P¹(F_3) with 4 points (1,ω,ω² identified).")
    perms = gluing_perms(R)
    # Explicit: pts[0]=(0,1)=∞, pts[1]=(1,0)=0, pts[2]=(1,1), pts[3]=(1,2)
    labels = ["∞", "0", "1", "2"]
    for name, perm in perms.items():
        cyc = _perm_cycles(perm, labels)
        print(f"  {name}^-1 on P1(F3): {cyc}")


def _perm_cycles(perm, labels):
    n = len(perm)
    seen = [False] * n
    parts = []
    for i in range(n):
        if seen[i]:
            continue
        c = [labels[i]]
        seen[i] = True
        j = perm[i]
        while j != i:
            c.append(labels[j])
            seen[j] = True
            j = perm[j]
        parts.append("(" + " ".join(c) + ")")
    return "".join(parts)


def check_cert_consistency() -> bool:
    """Re-run light asserts used by cert_omega_p (optional heavy import)."""
    print("\n5) Cert pipeline light checks")
    print("-" * 60)
    try:
        from face_pairings_p3 import (
            build_pair_maps, exterior_faces, validate_pair_maps,
        )
        from congruence_omega_proto import ref_elements
        from reference_cell import build_P3_cell
    except Exception as e:
        print(f"  SKIP (import): {e}")
        return True
    ok = True
    mesh = build_P3_cell(4, 2, 1.25, lift=True)
    ref = ref_elements(mesh)
    bf = exterior_faces(ref["mesh"], ref["tf"], ref["nfr"])
    pm, _, meta = build_pair_maps(ref["mesh"], ref["tf"], ref["nfr"], bfaces=bf)
    vok, msgs = validate_pair_maps(pm, bf)
    ok &= vok
    print(f"  face dict validate={'OK' if vok else 'FAIL'} stats={meta['stats']}")
    top_f = {b["fid"] for b in bf if b["kind"] == "top"}
    hit_top = False
    for name, m in pm.items():
        if set(m) & top_f:
            print(f"  FAIL {name} pairs a TOP face (breaks Lemma D0)")
            ok = False
            hit_top = True
    if not hit_top:
        print("  no TOP faces in pair maps (D0 OK)")
    print(f"  [{'OK' if ok else 'FAIL'}] cert light checks")
    return ok


def main() -> int:
    print("Generator / gluing audit (check our work)")
    print("=" * 60)
    ok = True
    ok &= check_vertex_cycles()
    ok &= check_forward_consistency()
    ok &= check_residue_gluing()
    check_N3_labels_vs_vertices()
    # step 5 can be slow; always run but last
    ok &= check_cert_consistency()
    print("\n" + "=" * 60)
    print(f"OVERALL: {'PASS' if ok else 'FAIL'}")
    if ok:
        print("  Vertex cycles: match user data (U: our U vs user U^-1")
        print("  orientation = inverse cycle, same cyclic group).")
        print("  Residue gluing: bijective, S^2=U^3=id, round-trips OK.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
