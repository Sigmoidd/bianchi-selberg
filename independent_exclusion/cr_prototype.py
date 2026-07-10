"""Crouzeix-Raviart guaranteed-lower-bound pipeline for the exclusion
criterion — floating-point prototype of Theorem G1 in lower_bound_theory.md.

Pipeline per window W = [s-, s+] of the s-grid:
  constants (Lemmas G, E, S)  ->  coefficients c_Q, c_Sigma, lam~, beta~, rho~,
  c_e, d_e (Thm G1 (i)-(ii))  ->  checks:
     (a) N_h - D_h  positive semidefinite   (Cholesky)
     (b) c_e / d_e >= 1
Both pass on all windows covering (0,1)  =>  (P_s) holds with nu* = 1+eps0
=> lambda_1 >= 1 (given DESIGN.md criterion + Lemma P).

Float arithmetic: this measures the post-constants margin (M2); M3 redoes the
checks in interval arithmetic.
"""

import numpy as np
from scipy.linalg import eigh, cholesky, LinAlgError

KAPPA1 = np.sqrt(1.0 / np.pi ** 2 + 1.0 / 120.0)   # CP22 sec 2.4, 3D CR
C_KAPPA = KAPPA1 ** 2 + (2.0 / 3.0) * KAPPA1        # Lemma T + (I1)
Y_M = 1.0 / np.sqrt(2.0)                            # min y on K

# Kuhn split of the unit cube into 6 tets (corner bits), consistent diagonals
KUHN = [
    [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)],
    [(0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 1)],
    [(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)],
    [(0, 0, 0), (0, 1, 0), (0, 1, 1), (1, 1, 1)],
    [(0, 0, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1)],
    [(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)],
]

# degree-2 tetrahedron quadrature (4 points)
QA, QB = (5.0 + 3.0 * np.sqrt(5.0)) / 20.0, (5.0 - np.sqrt(5.0)) / 20.0
TET_QP = np.array([[QA if i == j else QB for j in range(3)] for i in range(3)]
                  + [[QB, QB, QB]])           # barycentric coords 1..3 (0th = 1-sum)
TET_QW = np.full(4, 0.25)


def yf(x1, x2):
    return np.sqrt(1.0 - x1 ** 2 - x2 ** 2)


def build_mesh(N1, N2, N3, Y, curved=True):
    """Tet mesh of K_h (or of a flat box if curved=False). Returns dict."""
    x1g = np.linspace(-0.5, 0.5, N1 + 1)
    x2g = np.linspace(0.0, 0.5, N2 + 1)

    # floor node heights with lift (Lemma G)
    floor_y = np.empty((N1 + 1, N2 + 1))
    lift = np.zeros((N1 + 1, N2 + 1))
    if curved:
        base = yf(x1g[:, None], x2g[None, :])
        dz2 = (x1g[1] - x1g[0]) ** 2 + (x2g[1] - x2g[0]) ** 2  # cell z-diam^2
        cell_lift = np.empty((N1, N2))
        for i in range(N1):
            for j in range(N2):
                yf_min = min(base[i, j], base[i + 1, j], base[i, j + 1],
                             base[i + 1, j + 1])
                # Lemma G: sag <= max|D2yf| * d^2 / 8 (barycentric variance)
                # (1+1e-9) guards float rounding so the arb re-check in
                # m3_certify.py has strict margin
                cell_lift[i, j] = 0.125 * dz2 / yf_min ** 3 * (1 + 1e-9)
        for i in range(N1 + 1):
            for j in range(N2 + 1):
                adj = cell_lift[max(i - 1, 0):min(i + 1, N1),
                                max(j - 1, 0):min(j + 1, N2)]
                lift[i, j] = adj.max()
        floor_y = base + lift
        delta_bar = lift.max()          # Lemma G: delta <= max node lift
    else:
        floor_y[:] = 0.5
        delta_bar = 0.0

    # nodes
    nn = (N1 + 1) * (N2 + 1) * (N3 + 1)
    X = np.empty((nn, 3))
    is_floor = np.zeros(nn, bool)
    is_top = np.zeros(nn, bool)
    lift_node = np.zeros(nn)

    def nid(i, j, k):
        return (i * (N2 + 1) + j) * (N3 + 1) + k

    for i in range(N1 + 1):
        for j in range(N2 + 1):
            for k in range(N3 + 1):
                tt = k / N3
                n = nid(i, j, k)
                X[n] = (x1g[i], x2g[j], floor_y[i, j] * (1 - tt) + Y * tt)
                is_floor[n] = (k == 0)
                is_top[n] = (k == N3)
                if k == 0:
                    lift_node[n] = lift[i, j]

    # verify K_h subset K on floor faces (self-check of Lemma G)
    if curved:
        for i in range(N1):
            for j in range(N2):
                for (a, b, c) in (((i, j), (i + 1, j), (i + 1, j + 1)),
                                  ((i, j), (i, j + 1), (i + 1, j + 1))):
                    P = np.array([[x1g[a[0]], x2g[a[1]], floor_y[a]],
                                  [x1g[b[0]], x2g[b[1]], floor_y[b]],
                                  [x1g[c[0]], x2g[c[1]], floor_y[c]]])
                    for u in np.linspace(0, 1, 4):
                        for v in np.linspace(0, 1 - u, 4):
                            p = (1 - u - v) * P[0] + u * P[1] + v * P[2]
                            assert p[2] >= yf(p[0], p[1]) - 1e-13, "K_h not in K"

    # tets
    tets = []
    for i in range(N1):
        for j in range(N2):
            for k in range(N3):
                for tt in KUHN:
                    tets.append([nid(i + a, j + b, k + c) for (a, b, c) in tt])
    tets = np.array(tets)

    return dict(X=X, tets=tets, is_floor=is_floor, is_top=is_top,
                lift_node=lift_node, delta_bar=delta_bar, Y=Y, N=(N1, N2, N3))


def geometry(mesh):
    """Per-tet geometric data + face table."""
    X, tets = mesh["X"], mesh["tets"]
    nt = len(tets)
    vol = np.empty(nt)
    grads = np.empty((nt, 4, 3))     # gradients of barycentric lambdas
    hT = np.empty(nt)
    ymin = np.empty(nt)
    ymax = np.empty(nt)

    face_id = {}
    tet_faces = np.empty((nt, 4), int)     # face opposite local vertex a
    face_area = []
    face_nodes = []
    for e, tet in enumerate(tets):
        P = X[tet]
        A = np.hstack([np.ones((4, 1)), P])
        det = np.linalg.det(A)
        if det < 0:                       # fix orientation
            tet[[0, 1]] = tet[[1, 0]]
            mesh["tets"][e] = tet
            P = X[tet]
            A = np.hstack([np.ones((4, 1)), P])
            det = np.linalg.det(A)
        assert det > 1e-15, "degenerate tet"
        vol[e] = det / 6.0
        G = np.linalg.inv(A)
        grads[e] = G[1:4, :].T
        d = P[:, None, :] - P[None, :, :]
        hT[e] = np.sqrt((d ** 2).sum(-1)).max()
        ymin[e], ymax[e] = P[:, 2].min(), P[:, 2].max()
        for a in range(4):
            tri = tuple(sorted(np.delete(tet, a)))
            if tri not in face_id:
                face_id[tri] = len(face_id)
                p = X[list(tri)]
                face_area.append(0.5 * np.linalg.norm(
                    np.cross(p[1] - p[0], p[2] - p[0])))
                face_nodes.append(tri)
            tet_faces[e, a] = face_id[tri]

    wQ = 1.0 / ymax                # stiffness weight lower bound
    wM = 1.0 / ymin ** 3           # mass weight upper bound
    return dict(vol=vol, grads=grads, hT=hT, ymin=ymin, ymax=ymax,
                wQ=wQ, wM=wM, tet_faces=tet_faces, nfaces=len(face_id),
                face_area=np.array(face_area), face_nodes=face_nodes)


def assemble(mesh, geo, unit_weights=False):
    """Dense CR matrices Q_pc, M_pc; vectors t, a; classify faces."""
    X, tets = mesh["X"], mesh["tets"]
    nf = geo["nfaces"]
    Q = np.zeros((nf, nf))
    M = np.zeros((nf, nf))
    avec = np.zeros(nf)
    tvec = np.zeros(nf)
    Mloc = (np.full((4, 4), -1.0 / 20.0) + np.eye(4) * (9.0 / 20.0))

    top_tets = []       # (tet index, global face) with a top face
    floor_tets = []
    for e, tet in enumerate(tets):
        wq = 1.0 if unit_weights else geo["wQ"][e]
        wm = 1.0 if unit_weights else geo["wM"][e]
        gphi = -3.0 * geo["grads"][e]                 # grad of CR basis
        Se = wq * geo["vol"][e] * (gphi @ gphi.T)
        Me = wm * geo["vol"][e] * Mloc
        fid = geo["tet_faces"][e]
        Q[np.ix_(fid, fid)] += Se
        M[np.ix_(fid, fid)] += Me

        # a-vector: exact weight y^-3 by quadrature (degree 2)
        P = X[tet]
        for (qp, qw) in zip(TET_QP, TET_QW):
            lam = np.array([1.0 - qp.sum(), *qp])
            y = lam @ P[:, 2]
            phi = 1.0 - 3.0 * lam
            w = 1.0 if unit_weights else y ** -3
            avec[fid] += qw * geo["vol"][e] * w * phi

        # top / floor faces of this tet
        for a in range(4):
            tri = np.delete(tet, a)
            if mesh["is_top"][tri].all():
                tvec[fid[a]] += geo["face_area"][fid[a]]
                top_tets.append((e, fid[a]))
            if mesh["is_floor"][tri].all():
                floor_tets.append((e, fid[a], tri.copy()))

    return Q, M, avec, tvec, top_tets, floor_tets


def constants(mesh, geo, top_tets, floor_tets):
    """Lemma E and Lemma S constants (slab / column-mean forms)."""
    wQ, wM, hT, vol = geo["wQ"], geo["wM"], geo["hT"], geo["vol"]
    Y = mesh["Y"]
    X = mesh["X"]
    R_area = 0.5
    gamma = KAPPA1 * np.max(hT * np.sqrt(wM / wQ))
    alpha_h = KAPPA1 * np.sqrt(np.sum(wM ** 2 * vol * hT ** 2 / wQ))

    # (E-t) slab form, optimized over admissible H_t
    yhat_max = X[mesh["is_floor"], 2].max() if mesh["delta_bar"] > 0 else \
        X[mesh["is_floor"], 2].max()
    H_max = Y - yhat_max
    assert H_max > 0
    tau = np.inf
    H_t_best = None
    for H in np.linspace(0.02, H_max, 40):
        in_slab = geo["ymax"] > Y - H - 1e-12
        m1 = np.max(hT[in_slab] / np.sqrt(wQ[in_slab]))
        t1 = np.sqrt(R_area / H) * KAPPA1 * m1
        t2 = np.sqrt(H * R_area * Y / 3.0)
        if t1 + t2 < tau:
            tau, H_t_best = t1 + t2, H
    tau = float(tau)

    # Lemma S per floor face, H_s = maximal admissible column height
    db = mesh["delta_bar"]
    if db > 0:
        H_s = H_max
        S_M = S_Q = S_S = 0.0
        for (e, f, tri) in floor_tets:
            dhat = mesh["lift_node"][tri].max()
            ym_F = np.min(yf(X[tri, 0], X[tri, 1]))     # min y_f at vertices (concave)
            yhat_F = X[tri, 2].max()
            yp_F = yhat_F + H_s
            S_M = max(S_M, 2 * (dhat / H_s) * (yp_F / ym_F) ** 3)
            S_Q = max(S_Q, 2 * dhat * (dhat + H_s) * yp_F / ym_F ** 3)
            S_S = max(S_S, 2 * dhat * (dhat + H_s) * yhat_F / ym_F ** 3)
        V_S = 0.5 * db * Y_M ** -3
    else:
        S_M = S_Q = S_S = V_S = 0.0
    return dict(gamma=gamma, tau=tau, alpha_h=alpha_h, H_t=H_t_best,
                S_M=S_M, S_Q=S_Q, S_S=S_S, V_S=V_S)


def window_check(Q, M, avec, tvec, cst, Y, s_lo, s_hi, rho, nu_star,
                 theta=0.7, theta2=0.5, thetap=0.5, alpha=0.5,
                 diagnostics=False):
    """Theorem G1 for one window. Returns dict with pass/fail and numbers."""
    s0 = 0.5 * (s_lo + s_hi)
    lam_p = 1.0 - s_lo ** 2
    beta_p = 2.0 * (1.0 - s_lo) / Y ** 2
    kap0 = 1.0 / ((1.0 + s0) * Y ** 2)
    dkap = max(abs(1.0 / ((1.0 + s_lo) * Y ** 2) - kap0),
               abs(1.0 / ((1.0 + s_hi) * Y ** 2) - kap0))

    omega = 2 * rho * (1 / theta2 - 1) * cst["V_S"]
    c_Q = 1.0 - (omega + nu_star * lam_p) * cst["S_Q"]
    c_S = 1.0 - (omega + nu_star * lam_p) * cst["S_S"]
    lam_t = nu_star * lam_p * (1 + cst["S_M"]) + omega * cst["S_M"]
    beta_t = nu_star * beta_p + 2 * rho * (1 / theta2 - 1) * dkap ** 2
    rho_t = rho * (1 - theta2)

    sigma_h = cst["alpha_h"] + kap0 * cst["tau"]
    c_e = c_Q - rho_t * (1 / theta - 1) * sigma_h ** 2
    d_e = (lam_t * (1 + 1 / alpha) * cst["gamma"] ** 2
           + beta_t * (1 + 1 / thetap) * cst["tau"] ** 2)

    ell = avec + kap0 * tvec
    N = c_Q * Q + rho_t * (1 - theta) * np.outer(ell, ell)
    D = lam_t * (1 + alpha) * M + beta_t * (1 + thetap) * np.outer(tvec, tvec)

    res = dict(s_lo=s_lo, s_hi=s_hi, c_Q=c_Q, c_S=c_S, c_e=c_e, d_e=d_e,
               ratio_ed=c_e / d_e)
    try:
        cholesky(N - D, lower=True)
        res["psd"] = True
    except LinAlgError:
        res["psd"] = False
    res["pass"] = res["psd"] and c_S >= 0 and c_e > 0 and c_e / d_e >= 1.0
    if diagnostics:
        res["nu_h"] = eigh(N, D, eigvals_only=True, subset_by_index=[0, 0])[0]
    return res


def run_K(N1, N2, N3, Y=1.25, rho=55.0, nu_star=1.05, nwin=8,
          diagnostics=False):
    print(f"--- K_h mesh {N1}x{N2}x{N3} (tets={6*N1*N2*N3}), Y={Y}, "
          f"rho={rho}, nu*={nu_star} ---")
    mesh = build_mesh(N1, N2, N3, Y)
    geo = geometry(mesh)
    Q, M, avec, tvec, top_tets, floor_tets = assemble(mesh, geo)
    nf = geo["nfaces"]
    one = np.ones(nf)

    vol_true = sum(geo["vol"])      # Euclidean volume of K_h
    # Euclidean vol(K) by fine 2D quadrature of Y - y_f over R
    xs = np.linspace(-0.5, 0.5, 2001)
    ys = np.linspace(0.0, 0.5, 1001)
    volK = np.trapezoid(np.trapezoid(Y - yf(xs[:, None], ys[None, :]), ys), xs)
    db = mesh["delta_bar"]
    ok_vol = vol_true <= volK + 1e-9 and volK <= vol_true + db * 0.5 + 1e-9
    print(f"  faces (CR dofs): {nf};  delta_bar = {db:.2e}")
    print(f"  check: vol(K_h)={vol_true:.6f} <= vol(K)={volK:.6f}"
          f" <= vol(K_h)+delta*|R|={vol_true + db*0.5:.6f} : {ok_vol}")
    print(f"  check: |Q@1| = {np.abs(Q @ one).max():.2e}   t(1) = {tvec@one:.6f}")
    cst = constants(mesh, geo, top_tets, floor_tets)
    print(f"  constants: gamma={cst['gamma']:.4f} tau={cst['tau']:.4f} "
          f"(H_t={cst['H_t']:.3f}) alpha_h={cst['alpha_h']:.4f} "
          f"S_M={cst['S_M']:.2e} S_Q={cst['S_Q']:.2e} "
          f"S_S={cst['S_S']:.2e} V_S={cst['V_S']:.2e}")

    sgrid = np.linspace(0.0, 1.0, nwin + 1)
    all_pass = True
    print(f"  {'window':>16} {'c_Q':>7} {'c_e/d_e':>9} {'PSD':>5}"
          + ("  nu_h" if diagnostics else ""))
    for w in range(nwin):
        r = window_check(Q, M, avec, tvec, cst, Y, sgrid[w], sgrid[w + 1],
                         rho, nu_star, diagnostics=diagnostics)
        all_pass &= r["pass"]
        line = (f"  [{r['s_lo']:.3f},{r['s_hi']:.3f}] {r['c_Q']:7.4f} "
                f"{r['ratio_ed']:9.2f} {str(r['psd']):>5}")
        if diagnostics:
            line += f"  {r['nu_h']:.4f}"
        print(line)
    print(f"  ==> ALL WINDOWS PASS: {all_pass}  "
          f"(certified-shape statement: (P_s) with nu*={nu_star} on (0,1))")
    return all_pass


def validate_box(N1=8, N2=4, N3=4):
    """Plain Neumann Laplacian on a box: CR GLB <= exact (unit weights)."""
    print("--- validation: Neumann box [-.5,.5]x[0,.5]x[0.5,1.25], "
          "unit weights ---")
    mesh = build_mesh(N1, N2, N3, Y=1.25, curved=False)
    geo = geometry(mesh)
    Q, M, _, _, _, _ = assemble(mesh, geo, unit_weights=True)
    lam_cr = eigh(Q, M, eigvals_only=True, subset_by_index=[0, 4])
    hmax = geo["hT"].max()
    a, b, c = 1.0, 0.5, 0.75
    exact = sorted(np.pi ** 2 * (k1 ** 2 / a ** 2 + k2 ** 2 / b ** 2
                                 + k3 ** 2 / c ** 2)
                   for k1 in range(4) for k2 in range(4) for k3 in range(4))[:5]
    print(f"  {'CR':>10} {'GLB':>10} {'exact':>10}  GLB<=exact?")
    ok = True
    for lcr, ex in zip(lam_cr, exact):
        glb = lcr / (1 + KAPPA1 ** 2 * hmax ** 2 * lcr)
        good = glb <= ex + 1e-9
        ok &= good
        print(f"  {lcr:10.4f} {glb:10.4f} {ex:10.4f}  {good}")
    print(f"  box validation pass: {ok}")
    return ok


if __name__ == "__main__":
    ok_box = validate_box()
    print()
    # coarse with eigenvalue diagnostics, then production Cholesky-only
    run_K(8, 4, 4, diagnostics=True)
    print()
    run_K(12, 6, 6, diagnostics=False)
