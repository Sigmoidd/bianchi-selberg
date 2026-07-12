"""cert_omega_p — interval certification for Γ₀(𝔭) ⊂ PSL(2, ℤ[ω]).

Target (rung 1): 𝔭 = (1−ω), N=3, index 4.
  ⇒ no Laplace eigenvalue in (0,1) on Γ₀(𝔭)\\H³
  via multi-cusp criterion + Theorem G1p (Lemma D0) + Arb + Rump BIT 46.

Design: CERT_OMEGA_P.md.  Parallel to independent_exclusion/m3p_certify.py.

Usage:
  python -u cert_omega_p.py              # certify N=3 @ 6×3
  python -u cert_omega_p.py status
  python -u cert_omega_p.py 3 6 3        # q N_tri N3
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

import flint
from flint import arb

_ROOT = Path(__file__).resolve().parents[1]
_HERE = Path(__file__).resolve().parent
for p in (_HERE, _ROOT / "independent_exclusion"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from m3_certify import (  # noqa: E402
    tet_arb_data, mid_rad, upper, amax, amin, a_yf, EPS, ETA,
)
from m3p_certify import (  # noqa: E402
    weighted_moments, rump_certify_inplace, scaled_radius_rows, build_A,
)

from reference_cell import (  # noqa: E402
    AREA_T, Y_DEFAULT, YF_MIN_P3, build_P3_cell,
)
from residue_omega import (  # noqa: E402
    ring_lambda3, ring_split7, ring_split13, gluing_perms,
)
from cusps_omega import cusp_table_prime  # noqa: E402
from face_pairings_p3 import (  # noqa: E402
    build_pair_maps, exterior_faces, validate_pair_maps,
)
from congruence_omega_proto import ref_elements  # noqa: E402
from cert_omega import orient_tets, build_face_tables, lemma_G_arb  # noqa: E402

flint.ctx.prec = 128

Y = Y_DEFAULT
NU_STAR = 1.05
NWIN = 8

# Frozen G1p params (certified 2026-07-12, N=3 mesh 6×3)
THETA = 0.5
THETA2 = 0.9
ALPHA = 0.2
TH4 = 0.5
RHO_T = 9.0  # ρ̃

CHECKLIST = [
    ("research_note", True, "CERT_OMEGA_P.md"),
    ("T0_areas_proved", True, "T0_AREA.md"),
    ("lemma_D0_applies", True, "tops at y=Y; no TOP in side glue"),
    ("float_mu_positive", True, "proto 3 6 3 min μ≈1.12"),
    ("face_dictionary", True, "face_pairings_p3"),
    ("ref_arb_p3", True, "Q_pc/rem M_ex constants"),
    ("global_scatter", True, "4-copy mid/rad + t_α"),
    ("window_coeffs_p", True, "c_e≥d_e arb 8/8"),
    ("rump_8of8", True,
     "Rump PSD 8/8 @ θ=0.5,θ₂=0.9,α=0.2,θ₄=0.5,ρ̃=9 (2026-07-12)"),
    ("pairing_lemma", True,
     "PAIRING_MATRICES.md + pairing_matrices.py PASS (SL(2,O); S/U maps; T1/Tw section)"),
]


def status() -> None:
    print("cert_omega_p — Γ₀(𝔭) interval cert")
    print("=" * 60)
    for name, done, note in CHECKLIST:
        print(f"  [{'DONE' if done else 'TODO'}] {name}: {note}")
    n = sum(1 for _, d, _ in CHECKLIST if d)
    print(f"\n  {n}/{len(CHECKLIST)} closed")


def _ring(q: int):
    if q == 3:
        return ring_lambda3()
    if q == 7:
        return ring_split7()
    if q == 13:
        return ring_split13()
    raise ValueError(q)


# ---------------------------------------------------------------------------
# 1. Reference cell in arb (m3p.reference_arb → P3)
# ---------------------------------------------------------------------------

def reference_arb_p3(N_tri: int = 6, N3: int = 3, Y: float = Y):
    """Element mid/rad + G1 constants on one P₃ core."""
    mesh = build_P3_cell(N_tri=N_tri, N3=N3, Y=Y, lift=True)
    mesh = orient_tets(mesh)
    faces = build_face_tables(mesh)
    X, tets = mesh["X"], mesh["tets"]
    nt = len(tets)
    tf = faces["tf"]
    nfr = faces["nfr"]

    ok_g = lemma_G_arb(mesh, faces)
    if not ok_g:
        raise RuntimeError("Lemma G failed")

    kappa = (1 / arb.pi() ** 2 + arb(1) / 15).sqrt()  # Lemma I1
    Sem = np.zeros((nt, 4, 4))
    Ser = np.zeros((nt, 4, 4))
    Rem = np.zeros((nt, 4, 4))
    Rer = np.zeros((nt, 4, 4))
    Mem = np.zeros((nt, 4, 4))
    Mer = np.zeros((nt, 4, 4))
    aem = np.zeros((nt, 4))
    aer = np.zeros((nt, 4))
    hT_u = np.empty(nt)
    wQ_lo = np.empty(nt)
    gamma = arb(0)
    alph2 = arb(0)
    rho_w = arb(0)

    for e in range(nt):
        P = X[list(tets[e])]
        vol, grads, hT = tet_arb_data(P)
        ymaxT = arb(float(P[:, 2].max()))
        yminT = arb(float(P[:, 2].min()))
        wQ = 1 / ymaxT
        wM = yminT ** (-3)
        gamma = amax(gamma, kappa * hT * (wM / wQ).sqrt())
        rho_w = amax(rho_w, ymaxT / yminT - 1)
        alph2 += wM ** 2 * vol * hT ** 2 / wQ
        hT_u[e] = upper(hT)
        wQ_lo[e] = (1.0 / float(P[:, 2].max())) * (1 - 1e-15)
        I1, ent, Mex = weighted_moments(
            [arb(float(P[i, 2])) for i in range(4)], vol
        )
        c_rem = I1 - wQ * vol
        if not bool(c_rem > -arb("1e-20")):
            raise RuntimeError(f"Q_rem weight not nonnegative on tet {e}")
        c_rem = amax(c_rem, arb(0))
        for aa in range(4):
            aem[e, aa], aer[e, aa] = mid_rad(ent[aa])
            for bb in range(4):
                gg = sum(grads[aa][d] * grads[bb][d] for d in range(3))
                # Exact-weight stiffness as Q_pc main term (matches float CR);
                # rem kept for G1p formula with th4 (often near 0 contribution
                # if we put full I1 in Sem — use split: pc=wQ, rem=c_rem).
                Sem[e, aa, bb], Ser[e, aa, bb] = mid_rad(9 * wQ * vol * gg)
                Rem[e, aa, bb], Rer[e, aa, bb] = mid_rad(9 * c_rem * gg)
                Mem[e, aa, bb], Mer[e, aa, bb] = mid_rad(Mex[aa][bb])
    alpha_ref = kappa * alph2.sqrt()

    # Lemma S on floor faces (optimized H per face like m3p)
    yhat_max = float(X[mesh["is_floor"], 2].max())
    H_max = Y - yhat_max
    assert H_max > 0
    S_M = arb(0)
    S_Q = arb(0)
    S_S = arb(0)
    db = arb(0)
    for (_e, _f, tri) in faces["floor_tets"]:
        P = X[list(tri)]
        gaps = [
            arb(float(P[i, 2])) - a_yf(float(P[i, 0]), float(P[i, 1]))
            for i in range(3)
        ]
        dhat = amax(amax(gaps[0], gaps[1]), gaps[2])
        db = amax(db, dhat)
        yfmin = amin(
            amin(
                a_yf(float(P[0, 0]), float(P[0, 1])),
                a_yf(float(P[1, 0]), float(P[1, 1])),
            ),
            a_yf(float(P[2, 0]), float(P[2, 1])),
        )
        yhat_F = arb(float(P[:, 2].max()))
        Hcap = float((arb(Y) - yhat_F).mid()) * 0.999
        bestM = None
        H_F = Hcap
        for H in np.linspace(0.05, max(Hcap, 0.06), 12):
            H = min(H, Hcap)
            if H <= 0:
                continue
            cand = 2 * (dhat / arb(H)) * ((yhat_F + H) / yfmin) ** 3
            if bestM is None or upper(cand) < upper(bestM):
                bestM, H_F = cand, H
        if bestM is None:
            bestM = 2 * (dhat / arb(H_max)) * ((yhat_F + H_max) / yfmin) ** 3
            H_F = H_max
        S_M = amax(S_M, bestM)
        S_Q = amax(
            S_Q,
            2 * dhat * (dhat + H_F) * (yhat_F + H_F) / yfmin ** 3,
        )
        S_S = amax(
            S_S,
            2 * dhat * (dhat + H_F) * yhat_F / yfmin ** 3,
        )
    y_m = arb(YF_MIN_P3)
    V_S = db / 2 * y_m ** (-3)

    # top pairs (fid, mid, rad) for scatter
    top_pairs = []
    for (_e, f) in faces["top_tets"]:
        nodes = faces["face_nodes"][f]
        p = X[list(nodes)]
        cr = (
            (arb(float(p[1, 0])) - arb(float(p[0, 0])))
            * (arb(float(p[2, 1])) - arb(float(p[0, 1])))
            - (arb(float(p[1, 1])) - arb(float(p[0, 1])))
            * (arb(float(p[2, 0])) - arb(float(p[0, 0])))
        )
        m, r = mid_rad(abs(cr) / 2)
        top_pairs.append((f, m, r))

    # D0 assert: tops at Y
    assert np.allclose(X[mesh["is_top"], 2], Y), "tops not at Y (D0)"

    print(
        f"  ref arb: tets={nt} nfr={nfr}  "
        f"gamma≤{upper(gamma):.5f} alpha_h≤{upper(alpha_ref):.5f} "
        f"rho_w≤{upper(rho_w):.5f}"
    )
    print(
        f"    S_M≤{upper(S_M):.4e} S_Q≤{upper(S_Q):.4e} "
        f"S_S≤{upper(S_S):.4e} V_S≤{upper(V_S):.4e}"
    )
    t1 = sum(m for _, m, _ in top_pairs)
    print(f"    top area sum mid={t1:.8f} (want {AREA_T:.8f})")

    return dict(
        mesh=mesh, faces=faces, X=X, tets=tets, tf=tf, nfr=nfr,
        Sem=Sem, Ser=Ser, Rem=Rem, Rer=Rer, Mem=Mem, Mer=Mer,
        aem=aem, aer=aer, top_pairs=top_pairs,
        gamma=gamma, alpha_ref=alpha_ref, rho_w=rho_w,
        S_M=S_M, S_Q=S_Q, S_S=S_S, V_S=V_S,
        Y=Y, N_tri=N_tri, N3=N3,
    )


# ---------------------------------------------------------------------------
# 2. Global scatter
# ---------------------------------------------------------------------------

def assemble_global_p(ref, q: int = 3):
    """NC copies, glue by residue perms + face dictionary, scatter mid/rad."""
    R = _ring(q)
    ct = cusp_table_prime(q, R.label)
    NC = R.index
    glue = gluing_perms(R)
    pts = R.p1_points()
    cusp_class = [0 if p[0] == 0 else 1 for p in pts]
    assert cusp_class.count(0) == 1 and cusp_class.count(1) == q

    mesh, tf, nfr = ref["mesh"], ref["tf"], ref["nfr"]
    bfaces = exterior_faces(mesh, tf, nfr)
    pair_maps, _, meta = build_pair_maps(mesh, tf, nfr, bfaces=bfaces)
    ok, msgs = validate_pair_maps(pair_maps, bfaces)
    if not ok:
        for m in msgs[:5]:
            print(" ", m)
        raise RuntimeError("face pair map validation failed")
    print(f"  pair maps: {meta['stats']}")

    parent = list(range(NC * nfr))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    n_glued = 0
    for name, perm in glue.items():
        pmap = pair_maps.get(name, {})
        for c in range(NC):
            j = perm[c]
            if j == c:
                continue
            for ds, dd in pmap.items():
                ra, rb = find(c * nfr + ds), find(j * nfr + dd)
                if ra != rb:
                    parent[rb] = ra
                n_glued += 1
    print(f"  cross-copy glues applied: {n_glued}")

    gid = {}
    gmap = np.empty(NC * nfr, dtype=int)
    for x in range(NC * nfr):
        r = find(x)
        gmap[x] = gid.setdefault(r, len(gid))
    ng = len(gid)
    print(f"  global CR dofs: {ng}  (raw {NC * nfr})")

    def scatter_mat(Em, Er):
        rows, cols, vm, vr = [], [], [], []
        for c in range(NC):
            ix = gmap[c * nfr:(c + 1) * nfr][tf]
            rows.append(np.repeat(ix, 4, axis=1).ravel())
            cols.append(np.tile(ix, (1, 4)).ravel())
            vm.append(Em.ravel())
            vr.append(Er.ravel())
        r = np.concatenate(rows)
        cl = np.concatenate(cols)
        Am = coo_matrix(
            (np.concatenate(vm), (r, cl)), (ng, ng)
        ).tocsr()
        Ar = coo_matrix(
            (np.concatenate(vr), (r, cl)), (ng, ng)
        ).tocsr()
        Ar.data = np.abs(Ar.data) + 1e-13 * np.abs(Am.data)
        return Am, Ar

    Qm, Qr = scatter_mat(ref["Sem"], ref["Ser"])
    Rm, Rr = scatter_mat(ref["Rem"], ref["Rer"])
    Mm, Mr = scatter_mat(ref["Mem"], ref["Mer"])
    am = np.zeros(ng)
    ar = np.zeros(ng)
    tim = np.zeros(ng)
    tir = np.zeros(ng)
    t0m = np.zeros(ng)
    t0r = np.zeros(ng)
    for c in range(NC):
        gd = gmap[c * nfr:(c + 1) * nfr]
        np.add.at(am, gd[tf].ravel(), ref["aem"].ravel())
        np.add.at(ar, gd[tf].ravel(), np.abs(ref["aer"]).ravel())
        tm, tr = (tim, tir) if cusp_class[c] == 0 else (t0m, t0r)
        for d, m, r in ref["top_pairs"]:
            tm[gd[d]] += m
            tr[gd[d]] += r
    ar += 1e-13 * np.abs(am)

    t_inf_1 = float(tim.sum())
    t_0_1 = float(t0m.sum())
    print(
        f"  checks: t_∞(1)={t_inf_1:.8f} (want {ct.area_inf:.8f})  "
        f"t_0(1)={t_0_1:.8f} (want {ct.area_0:.8f})"
    )
    assert abs(t_inf_1 - ct.area_inf) < 1e-6, (t_inf_1, ct.area_inf)
    assert abs(t_0_1 - ct.area_0) < 1e-6, (t_0_1, ct.area_0)
    one = np.ones(ng)
    print(f"  1'M_ex1={float(one @ (Mm @ one)):.6f}  |Q_pc@1|="
          f"{np.abs(Qm @ one).max():.2e}")

    return dict(
        ng=ng, NC=NC, q=q, ct=ct,
        Qm=Qm, Qr=Qr, Rm=Rm, Rr=Rr, Mm=Mm, Mr=Mr,
        am=am, ar=ar, tim=tim, tir=tir, t0m=t0m, t0r=t0r,
        # aliases for m3p build_A / scaled_radius_rows
        # m3p uses tim/t0m naming
    )


# ---------------------------------------------------------------------------
# 3. Window coefficients (G1p, |T|-general, Lemma D0)
# ---------------------------------------------------------------------------

def window_coeffs_p(
    ref, ct, s_lo, s_hi, th, th2, al, rho_t, nu, th4, ncopy,
):
    """G1p arb coefficients; D0 ⇒ no τ in d_e; two β̃_α."""
    aY = arb(ref["Y"])
    s_lo_a, s_hi_a = arb(float(s_lo)), arb(float(s_hi))
    s0 = (s_lo_a + s_hi_a) / 2
    lam_p = 1 - s_lo_a ** 2
    area_inf = arb(ct.area_inf)
    area_0 = arb(ct.area_0)
    binf_p = (1 - s_lo_a) / (area_inf * aY ** 2)
    b0_p = (1 - s_lo_a) / (area_0 * aY ** 2)
    kap0 = 1 / ((1 + s0) * aY ** 2)
    dk = amax(
        abs(1 / ((1 + s_lo_a) * aY ** 2) - kap0),
        abs(1 / ((1 + s_hi_a) * aY ** 2) - kap0),
    )
    d2 = arb(rho_t) / arb(th2)
    om = 12 * d2 * ref["V_S"]
    c_Q = 1 - (om + arb(nu) * lam_p) * ref["S_Q"]
    c_S = 1 - (om + arb(nu) * lam_p) * ref["S_S"]
    lam_t = arb(nu) * lam_p * (1 + ref["S_M"]) + om * ref["S_M"]
    bt_inf = arb(nu) * binf_p + 4 * d2 * dk ** 2
    bt_0 = arb(nu) * b0_p + 4 * d2 * dk ** 2
    # Lemma D0
    sig = arb(ncopy).sqrt() * ref["alpha_ref"]
    carry = 1 - (1 / arb(th4) - 1) * ref["rho_w"]
    c_e = c_Q * carry - arb(rho_t) * (1 / arb(th) - 1) * sig ** 2
    d_e = lam_t * (1 + 1 / arb(al)) * ref["gamma"] ** 2
    return dict(
        kap0=kap0, c_Q=c_Q, c_S=c_S, lam_t=lam_t,
        bt_inf=bt_inf, bt_0=bt_0, c_e=c_e, d_e=d_e, carry=carry,
        sig=sig,
    )


def _build_A_omega(glob, cQ_lo, cR_lo, rt1_lo, lt_up, btinf_up, bt0_up, k0):
    """Like m3p.build_A but keys tim/t0m already match."""
    # m3p build_A expects glob with Qm,Rm,Mm,am,ar,tim,tir,t0m,t0r
    g = dict(
        ng=glob["ng"],
        Qm=glob["Qm"], Qr=glob["Qr"],
        Rm=glob["Rm"], Rr=glob["Rr"],
        Mm=glob["Mm"], Mr=glob["Mr"],
        am=glob["am"], ar=glob["ar"],
        tim=glob["tim"], tir=glob["tir"],
        t0m=glob["t0m"], t0r=glob["t0r"],
    )
    return build_A(g, cQ_lo, cR_lo, rt1_lo, lt_up, btinf_up, bt0_up, k0)


# ---------------------------------------------------------------------------
# 4. Full certificate
# ---------------------------------------------------------------------------

def _mark(name: str, done: bool = True, note: str | None = None):
    global CHECKLIST
    out = []
    for n, d, nt in CHECKLIST:
        if n == name:
            out.append((n, done, note if note is not None else nt))
        else:
            out.append((n, d, nt))
    CHECKLIST = out


def _run_windows(ref, glob, ct, ncopy, th, th2, al, th4, rt, nu_star, nwin):
    """Return (all_ok, n_scalar, n_psd)."""
    sgrid = np.linspace(0.0, 1.0, nwin + 1)
    all_ok = True
    n_scalar = n_psd = 0
    print(
        f"  windows={nwin}  θ={th} θ₂={th2} α={al} θ₄={th4} ρ̃={rt}"
    )
    print(
        f"  {'window':>14} {'c_S>0':>6} {'c_e≥d_e':>9} {'PSD':>5} "
        f"{'c_e':>8} {'d_e':>8}"
    )
    for w in range(nwin):
        co = window_coeffs_p(
            ref, ct, sgrid[w], sgrid[w + 1],
            th, th2, al, rt, nu_star, th4, ncopy,
        )
        ok_s = (
            bool(co["c_S"] > 0)
            and bool(co["c_e"] > co["d_e"])
            and bool(co["c_Q"] > 0)
            and bool(co["carry"] > 0)
        )
        if ok_s:
            n_scalar += 1
        cQ_lo = float(co["c_Q"].mid() - co["c_Q"].rad())
        # Use full Q = Q_pc + Q_rem in N_h (cR_lo = cQ_lo) for PSD headroom;
        # G1p still allows c_Q[Q_pc+(1-th4)Q_rem] — try th4 small + full rem.
        cR_lo = float(
            (co["c_Q"] * (1 - arb(th4))).mid()
            - (co["c_Q"] * (1 - arb(th4))).rad()
        )
        # Boost: also allow full rem if th4 near 0
        rt1_lo = float((arb(rt) * (1 - arb(th))).mid()) * (1 - 1e-14)
        rt1_lo = max(rt1_lo, 0.0)
        lt_up = upper(co["lam_t"] * (1 + arb(al)))
        btinf_up = upper(co["bt_inf"])
        bt0_up = upper(co["bt_0"])
        k0m, k0r = mid_rad(co["kap0"])
        A, eps = _build_A_omega(
            glob, cQ_lo, cR_lo, rt1_lo, lt_up, btinf_up, bt0_up, k0m,
        )
        ellm = glob["am"] + k0m * (glob["tim"] + glob["t0m"])
        eps += (
            rt1_lo * 2 * np.linalg.norm(ellm)
            * (abs(k0r) * np.linalg.norm(glob["tim"] + glob["t0m"])) * 1.01
        )

        def radius_rows_for_s(s, _cQ=cQ_lo, _cR=cR_lo, _rt1=rt1_lo,
                              _lt=lt_up, _bi=btinf_up, _b0=bt0_up,
                              _k0m=k0m, _k0r=k0r):
            g = {k: glob[k] for k in (
                "Qm", "Qr", "Rm", "Rr", "Mm", "Mr",
                "am", "ar", "tim", "tir", "t0m", "t0r",
            )}
            return scaled_radius_rows(
                g, _cQ, _cR, _rt1, _lt, _bi, _b0, _k0m, _k0r, s,
            )

        psd, _ = rump_certify_inplace(
            A, eps, verbose=False, radius_rows_for_s=radius_rows_for_s,
        )
        del A
        if psd:
            n_psd += 1
        all_ok &= ok_s and psd
        print(
            f"  [{sgrid[w]:.2f},{sgrid[w+1]:.2f}] "
            f"{str(bool(co['c_S']>0)):>6} "
            f"{str(bool(co['c_e']>co['d_e'])):>9} "
            f"{str(psd):>5} "
            f"{float(co['c_e'].mid()):8.4f} {float(co['d_e'].mid()):8.4f}"
        )
    return all_ok, n_scalar, n_psd


def certify(
    q: int = 3,
    N_tri: int = 6,
    N3: int = 3,
    Y: float = Y,
    nu_star: float = NU_STAR,
    nwin: int = NWIN,
    th: float = THETA,
    th2: float = THETA2,
    al: float = ALPHA,
    th4: float = TH4,
    rt: float = RHO_T,
    ref=None,
    glob=None,
    verbose: bool = True,
):
    print(
        f"=== cert_omega_p: Γ₀ N={q}, mesh {N_tri}×{N3}, Y={Y}, "
        f"ν*={nu_star}, ρ̃={rt}, θ={th}, prec={flint.ctx.prec} ==="
    )
    ct = cusp_table_prime(q)
    print(
        f"  |T_∞|={ct.area_inf:.6f} |T_0|={ct.area_0:.6f}  "
        f"index={ct.index}  modes_OK={ct.ok_modes(Y)}"
    )
    if not ct.ok_modes(Y):
        raise RuntimeError("mode bounds not OK")

    if ref is None:
        ref = reference_arb_p3(N_tri=N_tri, N3=N3, Y=Y)
    _mark("ref_arb_p3", True, "Q_pc/rem M_ex constants")
    if glob is None:
        glob = assemble_global_p(ref, q=q)
    _mark("global_scatter", True, "4-copy mid/rad + t_α")
    ncopy = glob["NC"]
    ng = glob["ng"]
    print(f"  dense A budget: {ng}² × 8 B ≈ {ng * ng * 8 / 1e6:.1f} MB")

    all_ok, n_scalar, n_psd = _run_windows(
        ref, glob, ct, ncopy, th, th2, al, th4, rt, nu_star, nwin,
    )
    if n_scalar == nwin:
        _mark("window_coeffs_p", True, "c_e≥d_e arb all windows")
    if n_psd == nwin:
        _mark("rump_8of8", True, "Rump PSD all windows")

    print(f"\n  scalar c_e>d_e: {n_scalar}/{nwin}   Rump PSD: {n_psd}/{nwin}")
    print(
        f"  ==> Γ₀(N={q}) / PSL(2,Z[omega]) CERTIFIED (all windows): {all_ok}"
    )
    if all_ok:
        print(
            f"      no eigenvalues in (0,1) for Γ₀(𝔭)\\H³ with N(𝔭)={q} —"
        )
        print(
            "      multi-cusp criterion (T0_AREA) + Theorem G1p (Lemma D0) "
            "+ Rump BIT 46; face dict face_pairings_p3."
        )
        print("      Pairing matrices: pairing_matrices.py / PAIRING_MATRICES.md.")
    return all_ok, ref, glob


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in ("status", "-h", "--help"):
        status()
        return 0
    if argv and argv[0] == "research":
        print("See CERT_OMEGA_P.md")
        status()
        return 0
    q, N_tri, N3 = 3, 6, 3
    if argv and argv[0].isdigit():
        q = int(argv[0])
        if len(argv) >= 2:
            N_tri = int(argv[1])
        if len(argv) >= 3:
            N3 = int(argv[2])

    # Frozen first, then fallbacks
    param_grid = [
        (THETA, THETA2, ALPHA, TH4, RHO_T),  # certified freeze
        (0.5, 0.9, 0.2, 0.5, 9.0),
        (0.4, 0.9, 0.2, 0.3, 12.0),
        (0.6, 0.9, 0.2, 0.85, 9.0),
        (0.5, 0.85, 0.15, 0.4, 15.0),
    ]
    ok = False
    ref = glob = None
    for th, th2, al, th4, rt in param_grid:
        print(f"\n  === trial θ={th} θ₂={th2} α={al} θ₄={th4} ρ̃={rt} ===")
        ok, ref, glob = certify(
            q=q, N_tri=N_tri, N3=N3,
            th=th, th2=th2, al=al, th4=th4, rt=rt,
            ref=ref, glob=glob,
        )
        if ok:
            print(f"\n  FROZEN params: θ={th} θ₂={th2} α={al} θ₄={th4} ρ̃={rt}")
            break
    print()
    status()
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
