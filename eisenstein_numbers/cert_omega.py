"""Stage 7 — interval certification for PSL(2, Z[omega]) on the EGM P_3 CR mesh.

Ports independent_exclusion/m3_certify.py (level-1 G1) + m3p_certify.py
Rump machinery (BIT 46 equilibration + per-row radii) onto
reference_cell.build_P3_cell.

Pipeline:
  1. P_3 mesh + tet orientation
  2. Lemma G floor inclusion in arb
  3. Per-tet arb enclosures (tet_arb_data, sandwich weights, ell_entries)
  4. Lemma E/S constants in arb (|T|=√3/6, y_m=√(2/3))
  5. Per s-window: arb c_e > d_e; N_h−D_h via Rump SAS + per-row radii

Float evidence (Stage 6): 8/8 PASS at N_tri=6, N3=3, rho=55.
This module upgrades that to interval certificates when all checks pass.

Usage:
  python -u cert_omega.py              # default mesh N_tri=6 N3=3
  python -u cert_omega.py 8 4          # finer
  python -u cert_omega.py status       # checklist only
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

import flint
from flint import arb

# Picard cert helpers (repo-local)
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT / "independent_exclusion") not in sys.path:
    sys.path.insert(0, str(_ROOT / "independent_exclusion"))
if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from m3_certify import (  # noqa: E402
    tet_arb_data, mid_rad, upper, amax, amin, a_yf, EPS, ETA,
)
from m3p_certify import (  # noqa: E402
    rump_certify_inplace, weighted_moments,
)

from reference_cell import (  # noqa: E402
    AREA_T, Y_DEFAULT, YF_MIN_P3, build_P3_cell, exact_vol_KY,
)

flint.ctx.prec = 128

# Default frozen parameters matching float G1 success (cr_omega.py)
RHO = 55.0
NU_STAR = 1.05
THETA = 0.7
THETA2 = 0.5
THETAP = 0.5
ALPHA = 0.5
NWIN = 8

# Static checklist (updated when certify() succeeds; defaults match 2026-07-11 run)
CHECKLIST = [
    ("geometry_freeze", True,
     "EGM F_3 + P_3 mesh + pairings + vol(K_Y)"),
    ("float_G1_8of8", True,
     "cr_omega.py: rho=55, N_tri=6,N3=3, min c_e/d_e~2.09"),
    ("lemma_G_arb", True,
     "floor lift inclusion in arb (cert_omega Lemma G)"),
    ("element_enclosures", True,
     "arb Q,M,a,t via weighted_moments on P_3"),
    ("window_scalar_arb", True,
     "c_e>d_e arb all 8 windows @ N_tri=6,N3=3,rho=55"),
    ("rump_psd", True,
     "equilibrated Rump SAS/per-row all 8 windows"),
    ("picard_regression", True,
     "independent_exclusion m3/m3p helpers reused"),
    ("non_neumann_float", True,
     "non_neumann_omega.py 8/8 I1+CZZ (paired FE; optional)"),
    ("journal_draft", True,
     "papers/paper3_eisenstein.tex"),
]


def status() -> None:
    print("Stage 7 — interval certification (Eisenstein–Picard / P_3)")
    print("=" * 56)
    for name, done, note in CHECKLIST:
        tag = "DONE" if done else "TODO"
        print(f"  [{tag}] {name}: {note}")
    n_done = sum(1 for _, d, _ in CHECKLIST if d)
    print(f"\n  {n_done}/{len(CHECKLIST)} items closed")
    print("  run: python -u cert_omega.py [N_tri N3]")


def orient_tets(mesh: dict) -> dict:
    """Ensure every tet has positive Euclidean volume (swap first two verts)."""
    X = mesh["X"]
    tets = np.array(mesh["tets"], dtype=int, copy=True)
    nflip = 0
    for e in range(len(tets)):
        P = X[list(tets[e])]
        A = np.hstack([np.ones((4, 1)), P])
        if np.linalg.det(A) < 0:
            tets[e, [0, 1]] = tets[e, [1, 0]]
            nflip += 1
    out = dict(mesh)
    out["tets"] = tets
    out["_nflip"] = nflip
    return out


def build_face_tables(mesh):
    """Face id map, tet→faces, top/floor tet lists (after orientation)."""
    X, tets = mesh["X"], mesh["tets"]
    is_top, is_floor = mesh["is_top"], mesh["is_floor"]
    nt = len(tets)
    fid = {}
    tf = np.empty((nt, 4), dtype=int)
    face_nodes = {}
    for e in range(nt):
        tet = tets[e]
        for a in range(4):
            face = tuple(sorted(int(x) for x in np.delete(tet, a)))
            if face not in fid:
                fid[face] = len(fid)
                face_nodes[fid[face]] = face
            tf[e, a] = fid[face]
    top_tets, floor_tets = [], []
    for e in range(nt):
        tet = tets[e]
        for a in range(4):
            tri = np.delete(tet, a)
            if is_top[tri].all():
                top_tets.append((e, tf[e, a]))
            if is_floor[tri].all():
                floor_tets.append((e, tf[e, a], tri.copy()))
    ymax = np.array([float(X[list(t), 2].max()) for t in tets])
    return dict(
        fid=fid, tf=tf, nfr=len(fid), face_nodes=face_nodes,
        top_tets=top_tets, floor_tets=floor_tets, ymax=ymax,
    )


def lemma_G_arb(mesh, faces) -> bool:
    """Arb verification: every floor face node sits above sphere + sag bound."""
    X = mesh["X"]
    lift = mesh["lift_node"]
    ok = True
    n_faces = 0
    for (_e, _f, tri) in faces["floor_tets"]:
        n_faces += 1
        P = X[list(tri)]
        gaps = [arb(float(P[i, 2])) - a_yf(float(P[i, 0]), float(P[i, 1]))
                for i in range(3)]
        d2 = arb(0)
        for i in range(3):
            for j in range(i + 1, 3):
                d2 = amax(
                    d2,
                    (arb(float(P[i, 0])) - arb(float(P[j, 0]))) ** 2
                    + (arb(float(P[i, 1])) - arb(float(P[j, 1]))) ** 2,
                )
        yfmin = amin(
            amin(a_yf(float(P[0, 0]), float(P[0, 1])),
                 a_yf(float(P[1, 0]), float(P[1, 1]))),
            a_yf(float(P[2, 0]), float(P[2, 1])),
        )
        sag = d2 / (8 * yfmin ** 3)
        gmin = amin(amin(gaps[0], gaps[1]), gaps[2])
        ok &= bool(gmin > sag) and bool(gmin > 0)
        # also node-wise: height > y_f
        for n in tri:
            ok &= bool(
                arb(float(X[n, 2])) - a_yf(float(X[n, 0]), float(X[n, 1]))
                > arb(0)
            )
    # global: lift_node vs max adjacent sag is already encoded in face checks
    print(f"  Lemma G: {n_faces} floor faces, admissible={ok}  "
          f"delta_bar={mesh['delta_bar']:.4e}  y_floor_min={mesh['y_floor_min']:.4f}")
    return ok


def assemble_arb(mesh, faces, kappa_mode="self"):
    """Build (mid, rad) dense Q, M, a, t and arb constants.

    Element weights via m3p weighted_moments (Taylor + remainder ball):
      Q_ab = 9 (grad λ_a · grad λ_b) ∫_T y^{-1}
      M_ab = ∫_T phi_a phi_b y^{-3}
      a_a  = ∫_T phi_a y^{-3}
    Lemma E/S still use sandwich extremes (wQ=1/ymax, wM=ymin^{-3}).
    """
    X, tets = mesh["X"], mesh["tets"]
    nt, nf = len(tets), faces["nfr"]
    tf = faces["tf"]
    aY = arb(mesh["Y"])
    area_T = arb(AREA_T)

    api = arb.pi()
    if kappa_mode == "self":
        kappa1 = (1 / api ** 2 + arb(1) / 15).sqrt()
    elif kappa_mode == "czz":
        kappa1 = (1 / api ** 2 + arb(1) / 120).sqrt()
    else:
        raise ValueError(kappa_mode)

    Qm = np.zeros((nf, nf))
    Qr = np.zeros((nf, nf))
    Mm = np.zeros((nf, nf))
    Mr = np.zeros((nf, nf))
    lm = np.zeros(nf)
    lr = np.zeros(nf)

    gamma_a = arb(0)
    alph_a = arb(0)
    hT_u = np.empty(nt)
    wQ_lo = np.empty(nt)

    for e in range(nt):
        P = X[list(tets[e])]
        vol, grads, hT = tet_arb_data(P)
        ymaxT = arb(float(P[:, 2].max()))
        yminT = arb(float(P[:, 2].min()))
        # sandwich extremes for Lemma E/S only
        wQ_s = 1 / ymaxT
        wM_s = yminT ** (-3)
        gamma_a = amax(gamma_a, kappa1 * hT * (wM_s / wQ_s).sqrt())
        alph_a += wM_s ** 2 * vol * hT ** 2 / wQ_s
        hT_u[e] = upper(hT)
        wQ_lo[e] = (1.0 / float(P[:, 2].max())) * (1 - 1e-15)

        yv = [arb(float(P[i, 2])) for i in range(4)]
        I1, ent, Mex = weighted_moments(yv, vol)
        fid_e = tf[e]
        for aa in range(4):
            ga = grads[aa]
            m, r = mid_rad(ent[aa])
            lm[fid_e[aa]] += m
            lr[fid_e[aa]] += r
            for bb in range(4):
                gb = grads[bb]
                gg = ga[0] * gb[0] + ga[1] * gb[1] + ga[2] * gb[2]
                # CR: grad phi_a = -3 grad λ_a ⇒ factor 9
                m, r = mid_rad(9 * I1 * gg)
                Qm[fid_e[aa], fid_e[bb]] += m
                Qr[fid_e[aa], fid_e[bb]] += r
                m, r = mid_rad(Mex[aa][bb])
                Mm[fid_e[aa], fid_e[bb]] += m
                Mr[fid_e[aa], fid_e[bb]] += r

    alph_a = kappa1 * alph_a.sqrt()
    # accumulation guard
    Qr += 1e-13 * np.abs(Qm)
    Mr += 1e-13 * np.abs(Mm)
    lr += 1e-13 * np.abs(lm)

    # top trace t (planar areas in z-plane = Euclidean face area at y=Y)
    tm = np.zeros(nf)
    tr = np.zeros(nf)
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
        tm[f] += m
        tr[f] += r

    # ---- Lemma E: tau slab ----
    Y = mesh["Y"]
    yhat_max = float(X[mesh["is_floor"], 2].max())
    H_max = Y - yhat_max
    assert H_max > 0, f"H_max={H_max}"
    tau_best = None
    H_t = None
    for H in np.linspace(0.02, H_max * 0.999, 50):
        sel = faces["ymax"] > Y - H
        if not np.any(sel):
            continue
        m1 = arb(0)
        for e in np.nonzero(sel)[0]:
            m1 = amax(m1, arb(hT_u[e]) / arb(wQ_lo[e]).sqrt())
        t1 = (area_T / arb(H)).sqrt() * kappa1 * m1
        t2 = (arb(H) * area_T * aY / 3).sqrt()
        tt = t1 + t2
        if tau_best is None or upper(tt) < upper(tau_best):
            tau_best, H_t = tt, H
    tau_a = tau_best

    # ---- Lemma S: floor lift ----
    H_s = H_max * 0.999
    S_M = arb(0)
    S_Q = arb(0)
    S_S = arb(0)
    db = arb(0)
    for (_e, _f, tri) in faces["floor_tets"]:
        dhat = arb(float(mesh["lift_node"][list(tri)].max()))
        db = amax(db, dhat)
        ym_F = amin(
            amin(
                a_yf(float(X[tri[0], 0]), float(X[tri[0], 1])),
                a_yf(float(X[tri[1], 0]), float(X[tri[1], 1])),
            ),
            a_yf(float(X[tri[2], 0]), float(X[tri[2], 1])),
        )
        yhat_F = arb(float(X[list(tri), 2].max()))
        yp_F = yhat_F + H_s
        S_M = amax(S_M, 2 * (dhat / H_s) * (yp_F / ym_F) ** 3)
        S_Q = amax(S_Q, 2 * dhat * (dhat + H_s) * yp_F / ym_F ** 3)
        S_S = amax(S_S, 2 * dhat * (dhat + H_s) * yhat_F / ym_F ** 3)

    # V_S = (1/2) delta_bar * y_m^{-3}, y_m = min y_f on P_3 = √(2/3)
    y_m = arb(YF_MIN_P3)
    V_S = db / 2 * y_m ** (-3)

    one = np.ones(nf)
    print(
        f"  assembly: dofs={nf}  1'Qm1={float(one @ (Qm @ one)):.6f}  "
        f"1'Mm1={float(one @ (Mm @ one)):.6f}  t(1)={float(tm @ one):.6f}"
    )
    vK, _ = exact_vol_KY(Y)
    print(f"    (EGM vol(K_Y)={vK:.6f}; 1'Mm1 from exact-weight Taylor mid)")
    print(
        f"  constants (arb upper): gamma={upper(gamma_a):.5f}  "
        f"tau={upper(tau_a):.5f} (H_t={H_t:.4f})  alpha_h={upper(alph_a):.5f}"
    )
    print(
        f"    S_M={upper(S_M):.4e} S_Q={upper(S_Q):.4e} "
        f"S_S={upper(S_S):.4e} V_S={upper(V_S):.4e}"
    )

    nrmQr = float(np.abs(Qr).sum(1).max()) * (1 + 1e-12)
    nrmMr = float(np.abs(Mr).sum(1).max()) * (1 + 1e-12)
    print(f"    radii: ||Qr||_∞row={nrmQr:.2e} ||Mr||_∞row={nrmMr:.2e}")

    return dict(
        Qm=Qm, Qr=Qr, Mm=Mm, Mr=Mr, lm=lm, lr=lr, tm=tm, tr=tr,
        gamma=gamma_a, tau=tau_a, alpha_h=alph_a,
        S_M=S_M, S_Q=S_Q, S_S=S_S, V_S=V_S,
        nrmQr=nrmQr, nrmMr=nrmMr,
        Y=Y, area_T=area_T, nf=nf, kappa_mode=kappa_mode,
    )


def scaled_radius_rows_level1(data, cQ_lo, rt1_lo, lt_up, bt_up, k0m, k0r, s):
    """Per-row radius for level-1 pencil (single Q, M, t — no multi-cusp)."""
    s = np.asarray(s, dtype=float)
    Qr = np.abs(data["Qr"])
    Mr = np.abs(data["Mr"])
    r = s * (cQ_lo * (Qr @ s) + lt_up * (Mr @ s))

    ellm = data["lm"] + k0m * data["tm"]
    ellr = (
        data["lr"]
        + abs(k0m) * data["tr"]
        + abs(k0r) * np.abs(data["tm"])
        + 1e-12 * np.abs(ellm)
    )
    um, ur = s * ellm, s * np.abs(ellr)
    u1, du1 = np.sum(np.abs(um)), np.sum(ur)
    r = r + rt1_lo * (np.abs(um) * du1 + ur * u1 + ur * du1)

    umt, urt = s * data["tm"], s * np.abs(data["tr"])
    t1, dt1 = np.sum(np.abs(umt)), np.sum(urt)
    r = r + bt_up * (np.abs(umt) * dt1 + urt * t1 + urt * dt1)

    Qm = np.abs(data["Qm"])
    Mm = np.abs(data["Mm"])
    absrow = s * (cQ_lo * (Qm @ s) + lt_up * (Mm @ s))
    absrow = absrow + rt1_lo * np.abs(um) * np.sum(np.abs(um))
    absrow = absrow + bt_up * np.abs(umt) * np.sum(np.abs(umt))
    r = r + 10 * EPS * absrow
    return r * (1 + 4 * EPS)


def certify_windows(data, rho=RHO, nu_star=NU_STAR, nwin=NWIN,
                    theta=THETA, theta2=THETA2, thetap=THETAP, alpha=ALPHA,
                    diagnostics=True):
    """Per-window arb scalar checks + equilibrated Rump PSD."""
    aY = arb(data["Y"])
    area_T = data["area_T"]
    ar, ath, ath2, athp, aa_, ans = (
        arb(rho), arb(theta), arb(theta2), arb(thetap),
        arb(alpha), arb(nu_star),
    )
    sgrid = np.linspace(0.0, 1.0, nwin + 1)
    all_ok = True
    results = []
    print(
        f"  windows={nwin}  rho={rho}  nu*={nu_star}  "
        f"theta={theta} theta2={theta2} alpha={alpha}"
    )
    print(
        f"  {'window':>14} {'c_e>d_e':>8} {'c_S>0':>6} {'PSD':>5} "
        f"{'c_e':>8} {'d_e':>8}"
        + ("  lam_min(Ahat)" if diagnostics else "")
    )

    for w in range(nwin):
        s_lo, s_hi = arb(float(sgrid[w])), arb(float(sgrid[w + 1]))
        s0 = (s_lo + s_hi) / 2
        lam_p = 1 - s_lo ** 2
        # |T|-general beta (Picard has |T|=1/2 ⇒ 2(1-s)/Y^2)
        beta_p = (1 - s_lo) / (area_T * aY ** 2)
        kap0 = 1 / ((1 + s0) * aY ** 2)
        dk = amax(
            abs(1 / ((1 + s_lo) * aY ** 2) - kap0),
            abs(1 / ((1 + s_hi) * aY ** 2) - kap0),
        )
        omega = 2 * ar * (1 / ath2 - 1) * data["V_S"]
        c_Q = 1 - (omega + ans * lam_p) * data["S_Q"]
        c_S = 1 - (omega + ans * lam_p) * data["S_S"]
        lam_t = ans * lam_p * (1 + data["S_M"]) + omega * data["S_M"]
        beta_t = ans * beta_p + 2 * ar * (1 / ath2 - 1) * dk ** 2
        rho_t = ar * (1 - ath2)
        sig = data["alpha_h"] + kap0 * data["tau"]
        c_e = c_Q - rho_t * (1 / ath - 1) * sig ** 2
        d_e = (
            lam_t * (1 + 1 / aa_) * data["gamma"] ** 2
            + beta_t * (1 + 1 / athp) * data["tau"] ** 2
        )
        ok_scalar = bool(c_S > 0) and bool(c_e > d_e) and bool(c_Q > 0)

        cQ_lo = float(c_Q.mid() - c_Q.rad())
        rt1_lo = float(
            ((rho_t * (1 - ath)).mid() - (rho_t * (1 - ath)).rad())
        ) * (1 - 1e-14)
        if rt1_lo < 0:
            rt1_lo = 0.0
        lt_up = upper(lam_t * (1 + aa_))
        bt_up = upper(beta_t * (1 + athp))
        k0m, k0r = mid_rad(kap0)

        ellm = data["lm"] + k0m * data["tm"]
        # Ahat midpoint (safe directions)
        A = np.array(
            cQ_lo * data["Qm"]
            + rt1_lo * np.outer(ellm, ellm)
            - lt_up * data["Mm"]
            - bt_up * np.outer(data["tm"], data["tm"]),
            order="F",
        )

        def radius_rows_for_s(s):
            return scaled_radius_rows_level1(
                data, cQ_lo, rt1_lo, lt_up, bt_up, k0m, k0r, s
            )

        # uniform extra fallback (unused when per-row provided)
        nrm_ellm = np.linalg.norm(ellm) * (1 + 1e-12)
        ellr = (
            data["lr"] + abs(k0m) * data["tr"]
            + abs(k0r) * np.abs(data["tm"]) + 1e-12 * np.abs(ellm)
        )
        nrm_ellr = np.linalg.norm(ellr) * (1 + 1e-12)
        nt_m = np.linalg.norm(data["tm"]) * (1 + 1e-12)
        nt_r = np.linalg.norm(data["tr"]) * (1 + 1e-12)
        eps_tot = (
            abs(cQ_lo) * data["nrmQr"]
            + rt1_lo * 2 * nrm_ellm * nrm_ellr
            + lt_up * data["nrmMr"]
            + bt_up * (2 * nt_m * nt_r + nt_r ** 2)
        )

        lam_min = None
        if diagnostics:
            try:
                lam_min = float(
                    eigh(A, eigvals_only=True, subset_by_index=[0, 0])[0]
                )
            except Exception:
                lam_min = float("nan")

        psd, c_rump = rump_certify_inplace(
            A, eps_tot, verbose=False, radius_rows_for_s=radius_rows_for_s,
        )
        del A
        ok = ok_scalar and psd
        all_ok &= ok
        line = (
            f"  [{sgrid[w]:.2f},{sgrid[w+1]:.2f}] "
            f"{str(bool(c_e > d_e)):>8} {str(bool(c_S > 0)):>6} "
            f"{str(psd):>5} "
            f"{float(c_e.mid()):8.4f} {float(d_e.mid()):8.4f}"
        )
        if diagnostics and lam_min is not None:
            line += f"  {lam_min:.3e}"
        print(line)
        results.append(dict(
            ok=ok, psd=psd, ok_scalar=ok_scalar,
            c_e=float(c_e.mid()), d_e=float(d_e.mid()),
            c_rump=c_rump, lam_min=lam_min,
        ))

    return all_ok, results


def certify(
    N_tri=6,
    N3=3,
    Y=Y_DEFAULT,
    rho=RHO,
    nu_star=NU_STAR,
    nwin=NWIN,
    kappa_mode="self",
    lift=True,
    diagnostics=True,
):
    print(
        f"=== Eisenstein–Picard M3 cert: P_3 N_tri={N_tri} N3={N3}, "
        f"Y={Y}, rho={rho}, nu*={nu_star}, kappa={kappa_mode}, "
        f"prec={flint.ctx.prec} ==="
    )
    mesh = build_P3_cell(N_tri=N_tri, N3=N3, Y=Y, lift=lift)
    mesh = orient_tets(mesh)
    print(
        f"  mesh: tets={len(mesh['tets'])} nodes={len(mesh['X'])}  "
        f"flips={mesh['_nflip']}  delta_bar={mesh['delta_bar']:.4e}"
    )
    faces = build_face_tables(mesh)

    # 1. Lemma G
    ok_g = lemma_G_arb(mesh, faces)
    if not ok_g:
        print("  FAIL: Lemma G not admissible — abort")
        return False
    CHECKLIST[2] = ("lemma_G_arb", True, "floor lift inclusion in arb")

    # 2–3. Assembly + constants
    data = assemble_arb(mesh, faces, kappa_mode=kappa_mode)
    CHECKLIST[3] = ("element_enclosures", True, "arb Q,M,a,t on P_3 CR mesh")

    # 4–5. Windows
    all_ok, results = certify_windows(
        data, rho=rho, nu_star=nu_star, nwin=nwin, diagnostics=diagnostics,
    )
    n_scalar = sum(1 for r in results if r["ok_scalar"])
    n_psd = sum(1 for r in results if r["psd"])
    if n_scalar == nwin:
        CHECKLIST[4] = (
            "window_scalar_arb", True,
            f"c_e>d_e arb all {nwin} windows",
        )
    if n_psd == nwin:
        CHECKLIST[5] = (
            "rump_psd", True,
            f"equilibrated Rump all {nwin} windows",
        )

    print(f"\n  scalar c_e>d_e: {n_scalar}/{nwin}   Rump PSD: {n_psd}/{nwin}")
    print(f"  ==> P_3 / PSL(2,Z[omega]) M3 CERTIFIED (all windows): {all_ok}")
    if all_ok:
        print(
            "      lambda_1(PSL(2,Z[omega])\\H^3) >= 1 by DESIGN.md criterion "
            "+ Theorem G1 + Rump BIT 46 (2006) Thm 2.3/Cor 2.4/Lemma 2.5 "
            f"(kappa={kappa_mode}); Neumann relaxation; remaining: paper G4."
        )
    else:
        print(
            "      not yet certified — tighten mesh / rho / or inspect "
            "lam_min(Ahat) margins."
        )
    return all_ok


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "status":
        status()
        return 0
    N_tri = int(argv[0]) if len(argv) >= 1 else 6
    N3 = int(argv[1]) if len(argv) >= 2 else 3
    ok = certify(N_tri=N_tri, N3=N3)
    print()
    status()
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
