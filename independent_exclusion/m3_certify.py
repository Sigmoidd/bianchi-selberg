"""M3: interval-arithmetic certification of Theorem G1's finite checks.

Everything Theorem G1 needs, re-verified rigorously:
  1. geometry admissibility in arb (tet volumes > 0, floor inclusion:
     node heights >= y_f + sag bound, Lemma G);
  2. per-tet quantities and all Lemma E / Lemma S constants as arb balls;
  3. the ell_0 vector enclosed by a degree-p Taylor expansion of y^-3
     around the tet mean height with a rigorous remainder ball;
  4. scalar window checks (c_Sigma > 0, c_e > d_e) by rigorous arb
     comparison;
  5. the matrix check N_h - D_h >= 0 via: interval terms reduced to
     lambda_min(A_hat) > eps_total with Lemma R (factored rank-one
     enclosures) and matrix-radius norms, then certified by a literal
     implementation of Rump, "Verification of positive definiteness",
     BIT Numer. Math. 46 (2006) 433-452: Theorem 2.3 (Cholesky
     completion => lambda_min > -||Delta(A)||_2), Corollary 2.4 and
     Lemma 2.5 (the phi-trick diagonal shift, no directed rounding).
     Rump op. cit. p.1-2: "any library routine can be used" -- we call
     LAPACK dpotrf via scipy, as INTLAB's isspd does.

Interpolation constant kappa: default is the self-contained Lemma I1
(lower_bound_theory.md), kappa_sc = sqrt(1/pi^2 + 1/15), valid for
arbitrary tetrahedra with only Payne-Weinberger (1960/Bebendorf 2003)
as input.  kappa_mode="czz" switches to the sharper published
sqrt(1/pi^2 + 1/120) (CP22 sec 2.4 / CZZ20).

The mesh coordinates are IEEE doubles = exact rationals; all arb inputs
are exact.  Float->arb->float extraction guards are documented inline.
"""

import numpy as np
import itertools
from math import factorial
from scipy.linalg import cholesky, LinAlgError, eigh

import flint
from flint import arb, arb_mat

from cr_prototype import build_mesh, geometry, assemble

flint.ctx.prec = 128
EPS = 2.0 ** -53        # Rump (1.1): relative rounding error unit, IEEE double
ETA = 2.0 ** -1074      # Rump (1.1): underflow unit

TAYLOR_P = 5            # ell_0 Taylor order


def rump_psd_certificate(A, extra):
    """Certify lambda_min(A) > extra for symmetric float A, after Rump,
    BIT 46 (2006), Thm 2.3 + Cor 2.4 + Lemma 2.5.

    Thm 2.3 constants: alpha_ij = gamma_{s(i,j)+2} <= gamma_{n+1} (dense
    pattern upper bound, Remark 1); d_j^2 = a_jj/(1-alpha_jj);
    M = 3(2n + max a_jj); Delta_ij = alpha_ij d_i d_j + M*eta;
    ||Delta||_2 <= gamma_{n+1} sum_j d_j^2 + n*M*eta  (entrywise
    domination + Perron-Frobenius, Remark 3).

    Cor 2.4 chain with c >= ||Delta(A)||_2 + extra: Cholesky completion on
    A~ (Lemma 2.5 diagonal, a~_ii <= a_ii - c) gives
    lambda_min(A) >= c + lambda_min(A~) > c - ||Delta(A~)||_2 >= extra.

    Returns (certified: bool, c: float).
    """
    n = A.shape[0]
    g = (n + 1) * EPS / (1.0 - (n + 1) * EPS)          # gamma_{n+1}
    assert g < 1.0
    diag = np.diag(A)
    d2 = np.maximum(diag, 0.0) / (1.0 - g)
    M = 3.0 * (2 * n + max(diag.max(), 0.0))
    # float-evaluation padding of the scalar bound itself: (1+n*EPS) on the
    # sum, (1+8*EPS) on products/adds
    c_delta = (g * d2.sum() * (1.0 + n * EPS) + n * M * ETA) * (1.0 + 8 * EPS)
    c = (c_delta + extra) * (1.0 + 8 * EPS)
    # Lemma 2.5: a~_ii = fl(d' - phi|d'|), d' = fl(a_ii - c)  =>  <= a_ii - c
    phi = EPS * (1.0 + 2.0 * EPS)
    dprime = diag - c
    At = A.copy()
    np.fill_diagonal(At, dprime - phi * np.abs(dprime))
    try:
        cholesky(At, lower=False)
        return True, c
    except LinAlgError:
        return False, c


def amax(a, b):
    return (a + b + abs(a - b)) / 2


def amin(a, b):
    return (a + b - abs(a - b)) / 2


def upper(a):
    """Rigorous float upper bound of arb a (outward guard on extraction)."""
    m, r = float(a.mid()), float(a.rad())
    return (m + r) * (1 + 1e-14) + 1e-300


def mid_rad(a):
    m = float(a.mid())
    r = float(a.rad()) * (1 + 1e-14) + abs(m) * 2.3e-16
    return m, r


def a_yf(x1, x2):
    return (arb(1) - arb(x1) ** 2 - arb(x2) ** 2).sqrt()


def tet_arb_data(P):
    """P: 4x3 float (exact). Returns vol, grads(4x3), hT, all arb."""
    A = arb_mat([[arb(1), arb(P[i, 0]), arb(P[i, 1]), arb(P[i, 2])]
                 for i in range(4)])
    det = A.det()
    vol = det / 6
    assert bool(vol > 0), "tet volume not certifiably positive"
    G = A.inv()
    grads = [[G[1, i], G[2, i], G[3, i]] for i in range(4)]
    h2 = arb(0)
    for i in range(4):
        for j in range(i + 1, 4):
            d = sum((arb(P[i, k]) - arb(P[j, k])) ** 2 for k in range(3))
            h2 = amax(h2, d)
    return vol, grads, h2.sqrt()


def ell_entries(yv, vol, p=TAYLOR_P):
    """arb enclosures of int_T phi_a y^-3 dx (a=0..3) by Taylor around ybar.

    int phi_a (y-ybar)^k = 6V k! sum_{|alpha|=k} prod d^alpha *
                           [1/(k+3)! - 3(alpha_a+1)/(k+4)!]
    g = y^-3, g^(k)/k! = (-1)^k (k+2)!/(2 k!) ybar^(-3-k).
    Remainder: |R_p| <= (p+3)!/(2 (p+1)!) ymin^(-4-p) * Dmax^(p+1),
    |phi_a| <= 2.
    """
    ybar = sum(yv) / 4
    d = [y - ybar for y in yv]
    out = [arb(0)] * 4
    for k in range(p + 1):
        gk = arb((-1) ** k) * factorial(k + 2) / (2 * factorial(k)) \
            * ybar ** (-3 - k)
        base = arb(0)
        corr = [arb(0)] * 4
        for alpha in itertools.product(range(k + 1), repeat=4):
            if sum(alpha) != k:
                continue
            mono = arb(factorial(k))
            for i in range(4):
                mono *= d[i] ** alpha[i] / factorial(alpha[i])
            # note: k!/prod(alpha!) * prod d^alpha, then the two brackets
            base += mono * (arb(factorial(k)) / factorial(k + 3))
            for a in range(4):
                corr[a] += mono * (arb(3 * factorial(k)) * (alpha[a] + 1)
                                   / factorial(k + 4))
        for a in range(4):
            out[a] += gk * 6 * vol * (base - corr[a])
    ymin = yv[0]
    for y in yv[1:]:
        ymin = amin(ymin, y)
    Dmax = amax(amax(yv[0], amax(yv[1], amax(yv[2], yv[3]))) - ybar,
                ybar - ymin)
    rem = (arb(factorial(p + 3)) / (2 * factorial(p + 1))
           * ymin ** (-4 - p) * Dmax ** (p + 1))
    ball = arb(0).union(2 * vol * rem).union(-2 * vol * rem)
    return [out[a] + ball for a in range(4)]


def certify(N1, N2, N3, Y=1.25, rho=55.0, nu_star=1.05, nwin=8,
            theta=0.7, theta2=0.5, thetap=0.5, alpha=0.5,
            kappa_mode="self", diagnostics=True):
    print(f"=== M3 certification: mesh {N1}x{N2}x{N3}, Y={Y}, rho={rho}, "
          f"nu*={nu_star}, Taylor p={TAYLOR_P}, prec={flint.ctx.prec}, "
          f"kappa={kappa_mode} ===")
    mesh = build_mesh(N1, N2, N3, Y)
    geo = geometry(mesh)                       # combinatorics reused; floats
    _, _, _, _, top_tets, floor_tets = assemble(mesh, geo)  # face classif.
    X, tets = mesh["X"], mesh["tets"]
    nf = geo["nfaces"]
    nt = len(tets)
    aY = arb(Y)

    # ---- 1. geometry admissibility (Lemma G) in arb -------------------
    N1_, N2_, N3_ = mesh["N"]
    x1g = np.linspace(-0.5, 0.5, N1_ + 1)
    x2g = np.linspace(0.0, 0.5, N2_ + 1)
    lift = mesh["lift_node"]

    def nid(i, j, k):
        return (i * (N2_ + 1) + j) * (N3_ + 1) + k

    ok_geom = True
    for i in range(N1_):
        for j in range(N2_):
            dz2 = arb(x1g[i + 1] - x1g[i]) ** 2 + arb(x2g[j + 1] - x2g[j]) ** 2
            yfm = amin(amin(a_yf(x1g[i], x2g[j]), a_yf(x1g[i + 1], x2g[j])),
                       amin(a_yf(x1g[i], x2g[j + 1]),
                            a_yf(x1g[i + 1], x2g[j + 1])))
            sag = dz2 / (8 * yfm ** 3)
            for (ii, jj) in ((i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)):
                ok_geom &= bool(arb(lift[nid(ii, jj, 0)]) > sag)
    # floor nodes exactly at yf+lift up to float rounding: verify node >= yf
    for i in range(N1_ + 1):
        for j in range(N2_ + 1):
            n = nid(i, j, 0)
            ok_geom &= bool(arb(X[n, 2]) - a_yf(X[n, 0], X[n, 1])
                            > arb(0))
    print(f"  geometry admissible (Lemma G, arb-verified): {ok_geom}")
    assert ok_geom

    # ---- 2. per-tet arb data, assembly with (mid, rad) ----------------
    Qm = np.zeros((nf, nf)); Qr = np.zeros((nf, nf))
    Mm = np.zeros((nf, nf)); Mr = np.zeros((nf, nf))
    lm = np.zeros(nf); lr = np.zeros(nf)
    gamma_a = arb(0)
    alph_a = arb(0)
    hT_u = np.empty(nt)          # float upper bounds of h_T (for tau slab max)
    wQ_lo = np.empty(nt)         # float lower bounds of w_T^Q
    api = arb.pi()
    if kappa_mode == "self":       # Lemma I1 (self-contained, arbitrary tets)
        kappa1 = (1 / api ** 2 + arb(1) / 15).sqrt()
    elif kappa_mode == "czz":      # CP22 sec 2.4 / CZZ20 (sharper, published)
        kappa1 = (1 / api ** 2 + arb(1) / 120).sqrt()
    else:
        raise ValueError(kappa_mode)

    MLOC = [[(arb(9) if aa == bb else arb(-1)) / 20 for bb in range(4)]
            for aa in range(4)]
    for e in range(nt):
        P = X[tets[e]]
        vol, grads, hT = tet_arb_data(P)
        ymaxT = arb(P[:, 2].max())            # exact double
        yminT = arb(P[:, 2].min())
        wQ = 1 / ymaxT
        wM = yminT ** -3
        gamma_a = amax(gamma_a, kappa1 * hT * (wM / wQ).sqrt())
        alph_a += wM ** 2 * vol * hT ** 2 / wQ
        hT_u[e] = upper(hT)
        wQ_lo[e] = 1.0 / P[:, 2].max()        # exact reciprocal bound dir?
        wQ_lo[e] = wQ_lo[e] * (1 - 1e-15)     # guard downward

        fid = geo["tet_faces"][e]
        for aa in range(4):
            ga = grads[aa]
            for bb in range(4):
                gb = grads[bb]
                s = 9 * wQ * vol * (ga[0] * gb[0] + ga[1] * gb[1]
                                    + ga[2] * gb[2])
                m, r = mid_rad(s)
                Qm[fid[aa], fid[bb]] += m; Qr[fid[aa], fid[bb]] += r
                m, r = mid_rad(wM * vol * MLOC[aa][bb])
                Mm[fid[aa], fid[bb]] += m; Mr[fid[aa], fid[bb]] += r
        ent = ell_entries([arb(P[i, 2]) for i in range(4)], vol)
        for aa in range(4):
            m, r = mid_rad(ent[aa])
            lm[fid[aa]] += m; lr[fid[aa]] += r
    alph_a = kappa1 * alph_a.sqrt()
    # assembly float-accumulation guard (<= ~60 contributions per entry)
    Qr += 1e-13 * np.abs(Qm); Mr += 1e-13 * np.abs(Mm)
    lr += 1e-13 * np.abs(lm)

    # t vector (top faces): arb areas
    tm = np.zeros(nf); tr = np.zeros(nf)
    for (e, f) in top_tets:
        nodes = list(geo["face_nodes"][f])
        p = X[nodes]
        cr = ((arb(p[1, 0]) - arb(p[0, 0])) * (arb(p[2, 1]) - arb(p[0, 1]))
              - (arb(p[1, 1]) - arb(p[0, 1])) * (arb(p[2, 0]) - arb(p[0, 0])))
        m, r = mid_rad(abs(cr) / 2)
        tm[f] += m; tr[f] += r

    # ---- 3. Lemma E / S constants in arb ------------------------------
    yhat_max = X[mesh["is_floor"], 2].max()      # exact double
    H_max = Y - yhat_max
    assert H_max > 0
    R_area = arb(1) / 2
    tau_best = None
    for H in np.linspace(0.02, H_max * 0.999, 40):
        sel = geo["ymax"] > Y - H                # float on exact data
        m1 = arb(0)
        for e in np.nonzero(sel)[0]:
            m1 = amax(m1, arb(hT_u[e]) / arb(wQ_lo[e]).sqrt())
        t1 = (R_area / arb(H)).sqrt() * kappa1 * m1
        t2 = (arb(H) * R_area * aY / 3).sqrt()
        tt = t1 + t2
        if tau_best is None or upper(tt) < upper(tau_best):
            tau_best, H_t = tt, H
    tau_a = tau_best

    H_s = H_max * 0.999
    S_M = arb(0); S_Q = arb(0); S_S = arb(0)
    db = arb(0)
    for (e, f, tri) in floor_tets:
        dhat = arb(float(mesh["lift_node"][list(tri)].max()))
        db = amax(db, dhat)
        ym_F = amin(amin(a_yf(X[tri[0], 0], X[tri[0], 1]),
                         a_yf(X[tri[1], 0], X[tri[1], 1])),
                    a_yf(X[tri[2], 0], X[tri[2], 1]))
        yhat_F = arb(X[list(tri), 2].max())
        yp_F = yhat_F + H_s
        S_M = amax(S_M, 2 * (dhat / H_s) * (yp_F / ym_F) ** 3)
        S_Q = amax(S_Q, 2 * dhat * (dhat + H_s) * yp_F / ym_F ** 3)
        S_S = amax(S_S, 2 * dhat * (dhat + H_s) * yhat_F / ym_F ** 3)
    V_S = db / 2 * arb(2).sqrt() ** 3            # (1/2) delta_bar y_m^-3

    print(f"  constants (arb, upper): gamma={upper(gamma_a):.5f} "
          f"tau={upper(tau_a):.5f} (H_t={H_t:.3f}) alpha_h={upper(alph_a):.5f}")
    print(f"    S_M={upper(S_M):.4e} S_Q={upper(S_Q):.4e} "
          f"S_S={upper(S_S):.4e} V_S={upper(V_S):.4e}")

    # norms for the matrix reduction (float uppers, guarded)
    nrmQr = np.abs(Qr).sum(1).max() * (1 + 1e-12)
    nrmMr = np.abs(Mr).sum(1).max() * (1 + 1e-12)
    nl_m = np.linalg.norm(lm) * (1 + 1e-12)
    nl_r = np.linalg.norm(lr) * (1 + 1e-12)
    nt_m = np.linalg.norm(tm) * (1 + 1e-12)
    nt_r = np.linalg.norm(tr) * (1 + 1e-12)
    print(f"    radii: ||Qr||={nrmQr:.2e} ||Mr||={nrmMr:.2e} "
          f"||l_rad||={nl_r:.2e} (||l||={nl_m:.4f}) ||t_rad||={nt_r:.2e}")

    # ---- 4./5. windows -------------------------------------------------
    ar, at, aa_, ath, ath2, athp, ans = (arb(rho), arb(theta), arb(alpha),
                                         arb(theta), arb(theta2), arb(thetap),
                                         arb(nu_star))
    sgrid = np.linspace(0.0, 1.0, nwin + 1)
    all_ok = True
    print(f"  {'window':>16} {'c_e>d_e':>8} {'PSD':>5} {'Rump shift c':>13}"
          + ("  diag lam_min" if diagnostics else ""))
    for w in range(nwin):
        s_lo, s_hi = arb(sgrid[w]), arb(sgrid[w + 1])
        s0 = (s_lo + s_hi) / 2
        lam_p = 1 - s_lo ** 2
        beta_p = 2 * (1 - s_lo) / aY ** 2
        kap0 = 1 / ((1 + s0) * aY ** 2)
        dk = amax(abs(1 / ((1 + s_lo) * aY ** 2) - kap0),
                  abs(1 / ((1 + s_hi) * aY ** 2) - kap0))
        omega = 2 * ar * (1 / ath2 - 1) * V_S
        c_Q = 1 - (omega + ans * lam_p) * S_Q
        c_S = 1 - (omega + ans * lam_p) * S_S
        lam_t = ans * lam_p * (1 + S_M) + omega * S_M
        beta_t = ans * beta_p + 2 * ar * (1 / ath2 - 1) * dk ** 2
        rho_t = ar * (1 - ath2)
        sig = alph_a + kap0 * tau_a
        c_e = c_Q - rho_t * (1 / ath - 1) * sig ** 2
        d_e = (lam_t * (1 + 1 / aa_) * gamma_a ** 2
               + beta_t * (1 + 1 / athp) * tau_a ** 2)
        ok_scalar = bool(c_S > 0) and bool(c_e > d_e) and bool(c_Q > 0)

        # coefficient bounds in the safe directions
        cQ_lo = float((c_Q.mid() - c_Q.rad()))
        rt_lo = float(((rho_t * (1 - ath)).mid() - (rho_t * (1 - ath)).rad()))
        lt_up = upper(lam_t * (1 + aa_))
        bt_up = upper(beta_t * (1 + athp))
        k0_m, k0_r = mid_rad(kap0)

        ellm = lm + k0_m * tm
        ellr = lr + k0_m * tr + (abs(k0_r) * (np.abs(tm) + tr))
        nem = np.linalg.norm(ellm) * (1 + 1e-12)
        ner = np.linalg.norm(ellr) * (1 + 1e-12)

        Ahat = (cQ_lo * Qm + rt_lo * np.outer(ellm, ellm)
                - lt_up * Mm - bt_up * np.outer(tm, tm))
        eps_tot = (cQ_lo * nrmQr + rt_lo * 2 * nem * ner
                   + lt_up * nrmMr + bt_up * (2 * nt_m * nt_r + nt_r ** 2))
        # float build error of Ahat (<= ~8 flops/entry)
        absA = (abs(cQ_lo) * np.abs(Qm) + abs(rt_lo) * np.outer(
            np.abs(ellm), np.abs(ellm)) + lt_up * np.abs(Mm)
            + bt_up * np.outer(np.abs(tm), np.abs(tm)))
        eps_tot += 8 * EPS * np.abs(absA).sum(1).max()

        # Rump BIT 46 (2006) Thm 2.3 + Cor 2.4 + Lemma 2.5
        psd, c_rump = rump_psd_certificate(Ahat, eps_tot)
        ok = ok_scalar and psd
        all_ok &= ok
        line = (f"  [{sgrid[w]:.3f},{sgrid[w+1]:.3f}] "
                f"{str(bool(c_e > d_e)):>8} {str(psd):>5} "
                f"{(c_rump if psd else 0.0):13.2e}")
        if diagnostics:
            evmin = eigh(Ahat, eigvals_only=True,
                         subset_by_index=[0, 0])[0]
            line += f"  {evmin:.3e}"
        print(line)
    print(f"  ==> M3 CERTIFIED (all windows): {all_ok}")
    if all_ok:
        print("      lambda_1(PSL(2,Z[i])\\H^3) >= 1, by DESIGN.md criterion"
              " + Theorem G1")
        print("      + Rump BIT 46 (2006) Thm 2.3/Cor 2.4/Lemma 2.5;"
              " kappa via " + ("self-contained Lemma I1"
                               if kappa_mode == "self" else "CZZ20/CP22")
              + ". Remaining: spectral-theory citations (G4).")
    return all_ok


if __name__ == "__main__":
    ok_self = certify(12, 6, 6, kappa_mode="self")
    print()
    ok_czz = certify(12, 6, 6, kappa_mode="czz", diagnostics=False)
