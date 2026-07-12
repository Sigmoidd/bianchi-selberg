"""Stages 5–6: CR assembly + Theorem G1 float probe for Eisenstein–Picard.

Ports independent_exclusion/cr_prototype.py to the omega reference cell:
  - CR face-mean dofs on tet mesh from reference_cell
  - proper Lemmas E/S constants (gamma, tau, alpha_h, S_Q/M/S, V_S)
  - per-window: c_e/d_e >= 1 and Cholesky PSD of N_h - D_h

Float only (Stage 7 = Arb+Rump). Geometry: GEOMETRY.md / geometry_fund.py.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.linalg import cholesky, LinAlgError, eigh
from scipy.sparse import coo_matrix

from reference_cell import AREA_T, Y_DEFAULT, build_reference_cell, y_floor
from constants import KAPPA_1, KAPPA_SC, beta_s, g1_coeffs

# deg-2 tet quadrature
QA, QB = (5 + 3 * math.sqrt(5)) / 20, (5 - math.sqrt(5)) / 20
TET_QP = np.array(
    [[QA if i == j else QB for j in range(3)] for i in range(3)]
    + [[QB, QB, QB]]
)
TET_QW = np.full(4, 0.25)


def _orient_vol(P):
    A = np.hstack([np.ones((4, 1)), P])
    det = np.linalg.det(A)
    if det < 0:
        P = P[[1, 0, 2, 3]]
        A = np.hstack([np.ones((4, 1)), P])
        det = np.linalg.det(A)
    return P, det / 6.0, np.linalg.inv(A)[1:4, :].T


def _tet_diameter(P):
    d2 = 0.0
    for a in range(4):
        for b in range(a + 1, 4):
            d2 = max(d2, float(np.sum((P[a] - P[b]) ** 2)))
    return math.sqrt(d2)


def assemble_cr(mesh):
    """CR global sparse Q, M and vectors a, t; plus per-tet geometry for G1."""
    X, tets = mesh["X"], mesh["tets"]
    nt = len(tets)
    fid = {}
    tf = np.empty((nt, 4), dtype=int)
    face_area = {}
    for e in range(nt):
        tet = tets[e]
        for a in range(4):
            face = tuple(sorted(np.delete(tet, a)))
            if face not in fid:
                fid[face] = len(fid)
                P3 = X[list(face)]
                face_area[fid[face]] = 0.5 * np.linalg.norm(
                    np.cross(P3[1] - P3[0], P3[2] - P3[0])
                )
            tf[e, a] = fid[face]
    nfr = len(fid)

    Se = np.zeros((nt, 4, 4))
    Me = np.zeros((nt, 4, 4))
    ae = np.zeros((nt, 4))
    Mloc = np.full((4, 4), -1 / 20.0) + np.eye(4) * (9 / 20.0)

    wQ = np.empty(nt)
    wM = np.empty(nt)
    vol = np.empty(nt)
    hT = np.empty(nt)
    ymin = np.empty(nt)
    ymax = np.empty(nt)
    grads = np.empty((nt, 4, 3))

    for e in range(nt):
        P = X[list(tets[e])]
        P, v, g = _orient_vol(P)
        vol[e] = v
        grads[e] = g
        hT[e] = _tet_diameter(P)
        ymin[e] = max(float(P[:, 2].min()), 1e-14)
        ymax[e] = float(P[:, 2].max())
        gphi = -3.0 * g
        wq = wm = 0.0
        av = np.zeros(4)
        for q, qw in zip(TET_QP, TET_QW):
            lam = np.array([1 - q.sum(), *q])
            y = max(float(lam @ P[:, 2]), 1e-14)
            wq += qw / y
            wm += qw / y ** 3
            av += qw / y ** 3 * (1 - 3 * lam)
        wQ[e], wM[e] = wq, wm
        Se[e] = wq * v * (gphi @ gphi.T)
        Me[e] = wm * v * Mloc
        ae[e] = v * av

    rows, cols, qv, mv = [], [], [], []
    avec = np.zeros(nfr)
    for e in range(nt):
        ix = tf[e]
        rows.append(np.repeat(ix, 4))
        cols.append(np.tile(ix, 4))
        qv.append(Se[e].ravel())
        mv.append(Me[e].ravel())
        np.add.at(avec, ix, ae[e])
    Q = coo_matrix(
        (np.concatenate(qv), (np.concatenate(rows), np.concatenate(cols))),
        (nfr, nfr),
    ).tocsr()
    M = coo_matrix(
        (np.concatenate(mv), (np.concatenate(rows), np.concatenate(cols))),
        (nfr, nfr),
    ).tocsr()

    # top trace t + floor/top tet lists for Lemma S / E
    tvec = np.zeros(nfr)
    top_tets = []
    floor_tets = []
    is_top, is_floor = mesh["is_top"], mesh["is_floor"]
    for e in range(nt):
        tet = tets[e]
        for a in range(4):
            tri = np.delete(tet, a)
            fa = tf[e, a]
            if is_top[tri].all():
                tvec[fa] += face_area[fa]
                top_tets.append((e, fa))
            if is_floor[tri].all():
                floor_tets.append((e, fa, tri.copy()))

    return dict(
        Q=Q, M=M, a=avec, t=tvec, nfr=nfr,
        AREA_T=mesh["AREA_T"], Y=mesh["Y"],
        wQ=wQ, wM=wM, vol=vol, hT=hT, ymin=ymin, ymax=ymax,
        grads=grads, tet_faces=tf, face_area=face_area,
        top_tets=top_tets, floor_tets=floor_tets,
        mesh=mesh,
    )


def compute_g1_constants(data, kappa=KAPPA_1):
    """Lemmas E + S constants (Picard cr_prototype.constants, |T|-adapted)."""
    mesh = data["mesh"]
    X = mesh["X"]
    Y = mesh["Y"]
    area_T = mesh["AREA_T"]
    wQ, wM, hT, vol = data["wQ"], data["wM"], data["hT"], data["vol"]
    ymax = data["ymax"]

    # gamma = kappa * max h_T * sqrt(wM/wQ)   (weighted CR interpolation)
    ratio = np.sqrt(np.maximum(wM, 1e-30) / np.maximum(wQ, 1e-30))
    gamma = float(kappa * np.max(hT * ratio))

    # alpha_h = kappa * sqrt( sum wM^2 vol h^2 / wQ )
    alpha_h = float(
        kappa * math.sqrt(
            np.sum(wM ** 2 * vol * hT ** 2 / np.maximum(wQ, 1e-30))
        )
    )

    # (E-t) slab form for tau, optimized over H_t
    yhat_max = float(X[mesh["is_floor"], 2].max())
    H_max = Y - yhat_max
    if H_max <= 0:
        raise RuntimeError(f"H_max={H_max}: floor reaches Y (need lift/Y fix)")
    tau = np.inf
    H_t_best = None
    for H in np.linspace(max(0.02, H_max * 0.05), H_max, 50):
        in_slab = ymax > Y - H - 1e-12
        if not np.any(in_slab):
            continue
        m1 = float(np.max(hT[in_slab] / np.sqrt(np.maximum(wQ[in_slab], 1e-30))))
        t1 = math.sqrt(area_T / H) * kappa * m1
        t2 = math.sqrt(H * area_T * Y / 3.0)
        if t1 + t2 < tau:
            tau, H_t_best = t1 + t2, H
    tau = float(tau)

    # Lemma S from floor faces (lift)
    db = float(mesh["delta_bar"])
    y_m = max(float(mesh.get("y_floor_min", 0.15)), 0.05)
    if db > 0 and data["floor_tets"]:
        H_s = H_max
        S_M = S_Q = S_S = 0.0
        lift_node = mesh["lift_node"]
        for (_e, _f, tri) in data["floor_tets"]:
            dhat = float(lift_node[tri].max())
            # min sphere height at the three corners (conservative)
            ym_F = min(y_floor(float(X[n, 0]), float(X[n, 1])) for n in tri)
            ym_F = max(ym_F, 1e-6)
            yhat_F = float(X[tri, 2].max())
            yp_F = yhat_F + H_s
            S_M = max(S_M, 2 * (dhat / H_s) * (yp_F / ym_F) ** 3)
            S_Q = max(S_Q, 2 * dhat * (dhat + H_s) * yp_F / ym_F ** 3)
            S_S = max(S_S, 2 * dhat * (dhat + H_s) * yhat_F / ym_F ** 3)
        V_S = 0.5 * db * y_m ** (-3)
    else:
        S_M = S_Q = S_S = V_S = 0.0

    return dict(
        gamma=gamma,
        tau=tau,
        alpha_h=alpha_h,
        H_t=H_t_best,
        S_M=S_M,
        S_Q=S_Q,
        S_S=S_S,
        V_S=V_S,
        y_m=y_m,
        H_max=H_max,
        kappa=kappa,
    )


def constrained_mu_dense(Q, M, a, t, area_T, Y, lam):
    """Float min of A on L=0 (dense; for modest CR dof counts)."""
    s = math.sqrt(1.0 - lam)
    beta = beta_s(s, area_T, Y)
    Qd = Q.toarray() if hasattr(Q, "toarray") else Q
    Md = M.toarray() if hasattr(M, "toarray") else M
    A = Qd - lam * Md - beta * np.outer(t, t)
    ell = a + t / ((1 + s) * Y ** 2)
    v = ell.copy()
    v[0] += np.sign(ell[0] if ell[0] != 0 else 1.0) * np.linalg.norm(ell)
    b = 2.0 / (v @ v)

    def hah(B):
        w = B @ v
        return (
            B - b * np.outer(v, w) - b * np.outer(w, v)
            + b ** 2 * (v @ w) * np.outer(v, v)
        )[1:, 1:]

    vals = eigh(hah(A), hah(Md), eigvals_only=True, subset_by_index=[0, 0])
    return float(vals[0])


def window_check(Q, M, a, t, co):
    """Thm G1 for one window: PSD of N-D and c_e/d_e >= 1."""
    ell = a + co["kap0"] * t
    th, thp, al = co["theta"], co["thetap"], co["alpha"]
    Nd = (
        co["c_Q"] * Q
        + co["rho_t"] * (1.0 - th) * np.outer(ell, ell)
    )
    Dd = (
        co["lam_t"] * (1.0 + al) * M
        + co["beta_t"] * (1.0 + thp) * np.outer(t, t)
    )
    ratio = co["c_e"] / max(co["d_e"], 1e-30)
    try:
        cholesky(Nd - Dd, lower=True)
        psd = True
        nu_h = None
    except LinAlgError:
        psd = False
        try:
            nu_h = float(
                eigh(Nd, Dd, eigvals_only=True, subset_by_index=[0, 0])[0]
            )
        except Exception:
            nu_h = float("nan")
    ok = (
        psd
        and co["c_S"] >= 0
        and co["c_e"] > 0
        and ratio >= 1.0
    )
    return dict(
        psd=psd, ratio=ratio, nu_h=nu_h, ok=ok,
        c_Q=co["c_Q"], c_e=co["c_e"], d_e=co["d_e"], c_S=co["c_S"],
    )


def window_float_checks(
    data, cst, nwin=8, rho=55.0, nu_star=1.05, do_mu=True,
):
    """Stage 6: per-window G1 float checks."""
    Y, area_T = data["Y"], data["AREA_T"]
    Q = data["Q"].toarray()
    M = data["M"].toarray()
    a, t = data["a"], data["t"]
    sgrid = np.linspace(0.0, 1.0, nwin + 1)
    results = []
    print(
        f"  windows={nwin}  rho={rho}  nu*={nu_star}\n"
        f"  gamma={cst['gamma']:.4f}  tau={cst['tau']:.4f}  "
        f"alpha_h={cst['alpha_h']:.4f}  H_t={cst['H_t']}\n"
        f"  S_M={cst['S_M']:.3e}  S_Q={cst['S_Q']:.3e}  "
        f"S_S={cst['S_S']:.3e}  V_S={cst['V_S']:.3e}"
    )
    print(f"  {'window':>14} {'c_Q':>7} {'c_e/d_e':>9} {'PSD':>5} {'mu':>9}")
    for w in range(nwin):
        co = g1_coeffs(
            sgrid[w], sgrid[w + 1], area_T, Y, cst,
            rho=rho, nu_star=nu_star,
        )
        r = window_check(Q, M, a, t, co)
        mu = float("nan")
        if do_mu and data["nfr"] <= 6000:
            lam = min(1.0 - sgrid[w] ** 2, 0.999)
            try:
                mu = constrained_mu_dense(Q, M, a, t, area_T, Y, lam)
            except Exception:
                mu = float("nan")
        tag = "PASS" if r["ok"] else "fail"
        print(
            f"  [{sgrid[w]:.2f},{sgrid[w+1]:.2f}] "
            f"{r['c_Q']:7.4f} {r['ratio']:9.3f} "
            f"{'yes' if r['psd'] else ' no':>5} {mu:9.4f}  {tag}"
        )
        results.append(dict(**r, mu=mu, co=co))
    return results


def run(
    N_tri=8, N3=4, Y=Y_DEFAULT, nwin=8,
    rho=55.0, nu_star=1.05, lift=True, do_mu=True,
    domain="P3", N1=None, N2=None,
):
    """N_tri = base edge subdivisions (P3); N3 = vertical layers.

    Legacy aliases: N1→N_tri, N2 ignored for P3.
    """
    if N1 is not None:
        N_tri = N1
    print("Eisenstein–Picard CR + G1 float probe (Stages 5–6)")
    print("=" * 64)
    print(
        f"reference cell domain={domain} N_tri={N_tri} N3={N3}, "
        f"Y={Y}, |T|={AREA_T:.6f}, lift={lift}"
    )
    mesh = build_reference_cell(N_tri, N_tri, N3, Y=Y, lift=lift, domain=domain)
    print(
        f"  tets={len(mesh['tets'])}  nodes={len(mesh['X'])}  "
        f"delta_bar={mesh['delta_bar']:.4e}  y_floor_min={mesh['y_floor_min']:.4f}"
    )
    data = assemble_cr(mesh)
    print(f"  CR dofs={data['nfr']}")
    one = np.ones(data["nfr"])
    print(
        f"  check 1'M1={one @ (data['M'] @ one):.6f}  "
        f"t(1)={data['t'] @ one:.6f} (want {AREA_T:.6f})  "
        f"|Q1|={np.abs(data['Q'] @ one).max():.2e}"
    )

    cst = compute_g1_constants(data, kappa=KAPPA_1)
    print(
        f"  G1 const: gamma={cst['gamma']:.4f} tau={cst['tau']:.4f} "
        f"alpha_h={cst['alpha_h']:.4f}"
    )

    if do_mu and data["nfr"] <= 4000:
        print("\n  M0-style mu (CR space, sample λ):")
        Qd, Md = data["Q"].toarray(), data["M"].toarray()
        for lam in (0.05, 0.5, 0.999):
            mu = constrained_mu_dense(
                Qd, Md, data["a"], data["t"],
                data["AREA_T"], data["Y"], lam,
            )
            print(f"    lam={lam:.3f}  mu={mu:.5f}")

    print("\n  Window G1 float checks:")
    res = window_float_checks(
        data, cst, nwin=nwin, rho=rho, nu_star=nu_star, do_mu=do_mu,
    )
    n_ok = sum(1 for r in res if r["ok"])
    n_mu = sum(1 for r in res if r.get("mu", 0) > 0)
    print(f"\n  G1 PASS (PSD + c_e/d_e>=1): {n_ok}/{nwin}")
    print(f"  windows with mu>0: {n_mu}/{nwin}")
    return data, cst, res


def scan_rho(N_tri=8, N3=4, Y=1.25, rhos=None, nwin=8, domain="P3"):
    """Find rho that maximises number of G1-passing windows (float)."""
    if rhos is None:
        rhos = [10, 20, 35, 55, 80, 120, 200]
    mesh = build_reference_cell(N_tri, N_tri, N3, Y=Y, lift=True, domain=domain)
    data = assemble_cr(mesh)
    cst = compute_g1_constants(data)
    Q, M = data["Q"].toarray(), data["M"].toarray()
    a, t = data["a"], data["t"]
    sgrid = np.linspace(0.0, 1.0, nwin + 1)
    print(
        f"rho scan P3 N_tri={N_tri} N3={N3}: gamma={cst['gamma']:.4f} "
        f"tau={cst['tau']:.4f} S_Q={cst['S_Q']:.3e}  dofs={data['nfr']}"
    )
    best = None
    for rho in rhos:
        n_ok = 0
        min_ratio = float("inf")
        all_psd = True
        for w in range(nwin):
            co = g1_coeffs(
                sgrid[w], sgrid[w + 1], data["AREA_T"], Y, cst, rho=rho,
            )
            r = window_check(Q, M, a, t, co)
            if r["ok"]:
                n_ok += 1
            min_ratio = min(min_ratio, r["ratio"])
            all_psd &= r["psd"]
        print(f"  rho={rho:6.1f}  pass={n_ok}/{nwin}  "
              f"min c_e/d_e={min_ratio:8.3f}  all_psd={all_psd}")
        if best is None or n_ok > best[0] or (
            n_ok == best[0] and min_ratio > best[1]
        ):
            best = (n_ok, min_ratio, rho)
    print(f"  best: rho={best[2]}  pass={best[0]}/{nwin}  "
          f"min_ratio={best[1]:.3f}")
    return best


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "scan":
        N = int(sys.argv[2]) if len(sys.argv) > 2 else 8
        N3 = int(sys.argv[3]) if len(sys.argv) > 3 else max(N // 2, 3)
        scan_rho(N, N3)
    else:
        print("=== rho scan (P3, N_tri=8, N3=4) ===")
        _n, _r, rho = scan_rho(8, 4)
        print("\n=== full run at best rho ===")
        run(N_tri=8, N3=4, Y=1.25, nwin=8, rho=rho, do_mu=True)
        if _n < 8:
            print("\n=== finer mesh probe N_tri=12 N3=6 ===")
            scan_rho(12, 6, rhos=[20, 35, 55, 80, 120, 200])
