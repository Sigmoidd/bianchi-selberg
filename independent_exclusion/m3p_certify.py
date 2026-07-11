"""M-p2: certified run of Theorem G1p for Gamma_0(p) (CONGRUENCE.md sec 7).

Reference-Cell Principle in action: all interval-enclosed data (element
matrices with sandwiched weights, ell_0 Taylor enclosures, Lemma E/G/S
constants) live on the ONE 24-split reference cell; the N(p)+1 copies
enter through exact coset gluing over R = Z[i]/p only.  The matrix
inequality is certified by the literal Rump (BIT 46, 2006) test, built
memory-carefully (single dense array, blocked rank-one updates, in-place
Cholesky).

Pipeline:
  1. reference cell in arb: geometry admissibility (Lemma G on quarter
     triangles), element enclosures, constants gamma/alpha_h/tau,
     sliver constants with per-face column heights;
  2. global (mid, rad) sparse matrices + vectors by exact scatter;
  3. float parameter search over (theta, theta2, theta', alpha, rho~)
     using the rigorous scalar formulas + Rump-in-float PSD trials;
  4. frozen parameters re-evaluated in arb; scalar checks by rigorous
     comparison; Rump certificate per window.

Level selection: call set_level('(3)') etc. before certify, or pass
level= to certify / CLI argv.
"""

import numpy as np
from scipy.linalg import cholesky, LinAlgError
from scipy.sparse import coo_matrix

import flint
from flint import arb

from congruence_prototype import (build_reference, ref_geometry,
                                  build_gluing, PAIRINGS, tri_key, Y,
                                  set_level)
import congruence_prototype as _cproto  # live NP/LEVEL after set_level
from m3_certify import (tet_arb_data, mid_rad, upper,
                        amax, amin, a_yf, EPS, ETA)

flint.ctx.prec = 128
NU_STAR = 1.02
NWIN = 8
TAYP = 5

import itertools
from math import factorial


def weighted_moments(yv, vol, p=TAYP):
    """Exact-weight enclosures on a tet with vertex heights yv (arb):
    I1 = int y^-1, ell[a] = int phi_a y^-3, Mex[a,b] = int phi_a phi_b y^-3.
    Degree-p Taylor about the mean height + rigorous remainder balls."""
    ybar = sum(yv) / 4
    d = [y - ybar for y in yv]
    I1 = arb(0)
    ell = [arb(0)] * 4
    Mex = [[arb(0)] * 4 for _ in range(4)]
    for k in range(p + 1):
        S0 = arb(0)
        S1 = [arb(0)] * 4
        S2 = [[arb(0)] * 4 for _ in range(4)]
        for beta in itertools.product(range(k + 1), repeat=4):
            if sum(beta) != k:
                continue
            Db = arb(1)
            for i in range(4):
                Db *= d[i] ** beta[i]
            S0 += Db
            for a in range(4):
                S1[a] += Db * (beta[a] + 1)
                for b in range(a, 4):
                    f = ((beta[a] + 1) * (beta[a] + 2) if a == b
                         else (beta[a] + 1) * (beta[b] + 1))
                    S2[a][b] += Db * f
        kf = factorial(k)
        g1 = arb((-1) ** k) * ybar ** (-1 - k)
        g3 = arb((-1) ** k) * factorial(k + 2) / (2 * kf) * ybar ** (-3 - k)
        I1 += g1 * 6 * vol * kf * S0 / factorial(k + 3)
        for a in range(4):
            ell[a] += g3 * 6 * vol * kf * (S0 / factorial(k + 3)
                                           - 3 * S1[a] / factorial(k + 4))
            for b in range(a, 4):
                Mex[a][b] += g3 * 6 * vol * kf * (
                    S0 / factorial(k + 3)
                    - 3 * (S1[a] + S1[b]) / factorial(k + 4)
                    + 9 * S2[a][b] / factorial(k + 5))
    ymin = amin(amin(yv[0], yv[1]), amin(yv[2], yv[3]))
    ymax = amax(amax(yv[0], yv[1]), amax(yv[2], yv[3]))
    Dm = amax(ymax - ybar, ybar - ymin)
    rem1 = vol * ymin ** (-2 - p) * Dm ** (p + 1)
    rem3 = (arb(factorial(p + 3)) / (2 * factorial(p + 1))
            * ymin ** (-4 - p) * Dm ** (p + 1))
    b1 = arb(0).union(rem1).union(-rem1)
    b2 = arb(0).union(2 * vol * rem3).union(-2 * vol * rem3)
    b4 = arb(0).union(4 * vol * rem3).union(-4 * vol * rem3)
    I1 += b1
    for a in range(4):
        ell[a] += b2
        for b in range(a, 4):
            Mex[a][b] += b4
            if b > a:
                Mex[b][a] = Mex[a][b]
    # consistency: row sums of Mex must overlap ell (sum_b phi_b = 1)
    for a in range(4):
        assert sum(Mex[a]).overlaps(ell[a]), "Mex/ell consistency"
    return I1, ell, Mex


# ----------------------------------------------------------------------
# 1. reference cell in arb
# ----------------------------------------------------------------------

def reference_arb(N1, N2, N3):
    X, tets, btri = build_reference(N1, N2, N3)
    vol_f, grads_f = ref_geometry(X, tets)     # fixes orientations in tets
    nt = len(tets)

    fid = {}
    tf = np.empty((nt, 4), int)
    for e in range(nt):
        for a in range(4):
            key = tuple(sorted(np.delete(tets[e], a)))
            tf[e, a] = fid.setdefault(key, len(fid))
    nfr = len(fid)

    kappa = (1 / arb.pi() ** 2 + arb(1) / 15).sqrt()   # Lemma I1
    MLOC = [[(arb(9) if aa == bb else arb(-1)) / 20 for bb in range(4)]
            for aa in range(4)]

    Sem = np.zeros((nt, 4, 4)); Ser = np.zeros((nt, 4, 4))     # Q_pc
    Rem = np.zeros((nt, 4, 4)); Rer = np.zeros((nt, 4, 4))     # Q_rem
    Mem = np.zeros((nt, 4, 4)); Mer = np.zeros((nt, 4, 4))     # M_ex
    aem = np.zeros((nt, 4)); aer = np.zeros((nt, 4))           # ell_0 part a
    hT_u = np.empty(nt)
    wQ_lo = np.empty(nt)
    gamma = arb(0)
    alph2 = arb(0)
    rho_w = arb(0)
    for e in range(nt):
        P = X[list(tets[e])]
        vol, grads, hT = tet_arb_data(P)
        ymaxT, yminT = arb(P[:, 2].max()), arb(P[:, 2].min())
        wQ, wM = 1 / ymaxT, yminT ** -3
        gamma = amax(gamma, kappa * hT * (wM / wQ).sqrt())
        rho_w = amax(rho_w, ymaxT / yminT - 1)
        alph2 += wM ** 2 * vol * hT ** 2 / wQ
        hT_u[e] = upper(hT)
        wQ_lo[e] = (1.0 / P[:, 2].max()) * (1 - 1e-15)
        I1, ent, Mex = weighted_moments([arb(P[i, 2]) for i in range(4)],
                                        vol)
        c_rem = I1 - wQ * vol          # >= 0: int (y^-1 - wQ)
        assert bool(c_rem > 0), "Q_rem cell weight not positive"
        for aa in range(4):
            aem[e, aa], aer[e, aa] = mid_rad(ent[aa])
            for bb in range(4):
                gg = sum(grads[aa][d] * grads[bb][d] for d in range(3))
                Sem[e, aa, bb], Ser[e, aa, bb] = mid_rad(9 * wQ * vol * gg)
                Rem[e, aa, bb], Rer[e, aa, bb] = mid_rad(9 * c_rem * gg)
                Mem[e, aa, bb], Mer[e, aa, bb] = mid_rad(Mex[aa][bb])
    alpha_ref = kappa * alph2.sqrt()

    # tau (slab form) on the reference cell
    yhat_max = X[np.array(sorted({n for t in btri["floor"] for n in t})), 2].max()
    H_max = Y - yhat_max
    assert H_max > 0
    ymax_f = np.array([X[list(t), 2].max() for t in tets])
    tau = None
    for H in np.linspace(0.02, H_max * 0.999, 40):
        sel = ymax_f > Y - H
        m1 = arb(0)
        for e in np.nonzero(sel)[0]:
            m1 = amax(m1, arb(hT_u[e]) / arb(wQ_lo[e]).sqrt())
        tt = ((arb(1) / 2) / arb(H)).sqrt() * kappa * m1 \
            + (arb(H) * (arb(1) / 2) * arb(Y) / 3).sqrt()
        if tau is None or upper(tt) < upper(tau):
            tau = tt

    # Lemma G admissibility + Lemma S constants on floor quarter-triangles
    ok = True
    S_M = arb(0); S_Q = arb(0); S_S = arb(0); db = arb(0)
    for t in btri["floor"]:
        P = X[list(t)]
        gaps = [arb(P[i, 2]) - a_yf(P[i, 0], P[i, 1]) for i in range(3)]
        d2 = arb(0)
        for i in range(3):
            for j in range(i + 1, 3):
                d2 = amax(d2, (arb(P[i, 0]) - arb(P[j, 0])) ** 2
                          + (arb(P[i, 1]) - arb(P[j, 1])) ** 2)
        yfmin = amin(amin(a_yf(P[0, 0], P[0, 1]), a_yf(P[1, 0], P[1, 1])),
                     a_yf(P[2, 0], P[2, 1]))
        sag = d2 / (8 * yfmin ** 3)
        gmin = amin(amin(gaps[0], gaps[1]), gaps[2])
        ok &= bool(gmin > sag) and bool(gmin > 0)
        dhat = amax(amax(gaps[0], gaps[1]), gaps[2])
        db = amax(db, dhat)
        yhat_F = arb(P[:, 2].max())
        Hcap = float((arb(Y) - yhat_F).mid()) * 0.999
        bestM = None
        for H in np.linspace(0.05, max(Hcap, 0.06), 12):
            H = min(H, Hcap)
            cand = 2 * (dhat / arb(H)) * ((yhat_F + H) / yfmin) ** 3
            if bestM is None or upper(cand) < upper(bestM):
                bestM, H_F = cand, H
        S_M = amax(S_M, bestM)
        S_Q = amax(S_Q, 2 * dhat * (dhat + H_F) * (yhat_F + H_F) / yfmin ** 3)
        S_S = amax(S_S, 2 * dhat * (dhat + H_F) * yhat_F / yfmin ** 3)
    V_S = db / 2 * arb(2).sqrt() ** 3
    assert ok, "Lemma G admissibility failed"

    top_pairs = []
    for t in btri["top"]:
        P = X[list(t)]
        cr = ((arb(P[1, 0]) - arb(P[0, 0])) * (arb(P[2, 1]) - arb(P[0, 1]))
              - (arb(P[1, 1]) - arb(P[0, 1])) * (arb(P[2, 0]) - arb(P[0, 0])))
        m, r = mid_rad(abs(cr) / 2)
        top_pairs.append((fid[tuple(sorted(t))], m, r))

    return dict(X=X, tets=tets, btri=btri, fid=fid, tf=tf, nfr=nfr,
                Sem=Sem, Ser=Ser, Rem=Rem, Rer=Rer,
                Mem=Mem, Mer=Mer, aem=aem, aer=aer,
                gamma=gamma, alpha_ref=alpha_ref, tau=tau, rho_w=rho_w,
                S_M=S_M, S_Q=S_Q, S_S=S_S, V_S=V_S, top_pairs=top_pairs)


# ----------------------------------------------------------------------
# 2. global scatter (exact combinatorics)
# ----------------------------------------------------------------------

def assemble_global(ref):
    X, tets, btri = ref["X"], ref["tets"], ref["btri"]
    fid, tf, nfr = ref["fid"], ref["tf"], ref["nfr"]
    pair_maps = {}
    for name, src, dst, mp in PAIRINGS:
        dst_lookup = {tri_key(X, t): fid[tuple(sorted(t))]
                      for t in btri[dst]}
        pair_maps[name] = {fid[tuple(sorted(t))]:
                           dst_lookup[tri_key(X, t, mapping=mp)]
                           for t in btri[src]}
    pts, glue, cusp_class = build_gluing()
    NC = getattr(_cproto, "INDEX", _cproto.NP + 1)
    parent = list(range(NC * nfr))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for name, src, dst, mp in PAIRINGS:
        perm = glue[name]
        for c in range(NC):
            j = perm[c]
            if j == c:
                continue
            for ds, dd in pair_maps[name].items():
                parent[find(c * nfr + ds)] = find(j * nfr + dd)
    gid = {}
    gmap = np.empty(NC * nfr, int)
    for x in range(NC * nfr):
        r = find(x)
        gmap[x] = gid.setdefault(r, len(gid))
    ng = len(gid)

    def scatter_mat(Em, Er):
        rows, cols, vm, vr = [], [], [], []
        for c in range(NC):
            ix = gmap[c * nfr:(c + 1) * nfr][tf]
            rows.append(np.repeat(ix, 4, axis=1).ravel())
            cols.append(np.tile(ix, (1, 4)).ravel())
            vm.append(Em.ravel()); vr.append(Er.ravel())
        r = np.concatenate(rows); cl = np.concatenate(cols)
        Am = coo_matrix((np.concatenate(vm), (r, cl)), (ng, ng)).tocsr()
        Ar = coo_matrix((np.concatenate(vr), (r, cl)), (ng, ng)).tocsr()
        Ar.data = np.abs(Ar.data) + 1e-13 * np.abs(Am.data)  # accum guard
        return Am, Ar

    Qm, Qr = scatter_mat(ref["Sem"], ref["Ser"])
    Rm, Rr = scatter_mat(ref["Rem"], ref["Rer"])
    Mm, Mr = scatter_mat(ref["Mem"], ref["Mer"])
    am = np.zeros(ng); ar = np.zeros(ng)
    tim = np.zeros(ng); tir = np.zeros(ng)
    t0m = np.zeros(ng); t0r = np.zeros(ng)
    for c in range(NC):
        gd = gmap[c * nfr:(c + 1) * nfr]
        np.add.at(am, gd[tf].ravel(), ref["aem"].ravel())
        np.add.at(ar, gd[tf].ravel(), np.abs(ref["aer"]).ravel())
        tm, tr = (tim, tir) if cusp_class[c] == 0 else (t0m, t0r)
        for d, m, r in ref["top_pairs"]:
            tm[gd[d]] += m
            tr[gd[d]] += r
    ar += 1e-13 * np.abs(am)
    ncopy = getattr(_cproto, "INDEX", _cproto.NP + 1)
    print(f"  global CR dofs: {ng}  level={_cproto.LEVEL} "
          f"N(n)={_cproto.NP} copies={ncopy}")
    print(f"  checks: t_inf(1)={tim.sum():.6f} (0.5) "
          f"t_0(1)={t0m.sum():.6f} ({_cproto.NP / 2} if 2-cusp prime) "
          f"1'M_ex 1={np.ones(ng) @ (Mm @ np.ones(ng)):.6f} "
          f"({ncopy}*vol_w(K_ref))")
    return dict(ng=ng, Qm=Qm, Qr=Qr, Rm=Rm, Rr=Rr, Mm=Mm, Mr=Mr,
                am=am, ar=ar, tim=tim, tir=tir, t0m=t0m, t0r=t0r)


# ----------------------------------------------------------------------
# 3./4. window coefficients, Rump certificate
# ----------------------------------------------------------------------

def window_coeffs(ref, s_lo, s_hi, th, th2, al, rho_t, nu, th4):
    """All Theorem G1p coefficients as arb. rho_t = rho~ = rho(1-theta2).
    Trace functionals are exact under I_CR (Lemma D0): no tau, no theta'."""
    aY = arb(Y)
    s_lo_a, s_hi_a = arb(s_lo), arb(s_hi)
    s0 = (s_lo_a + s_hi_a) / 2
    lam_p = 1 - s_lo_a ** 2
    binf_p = 2 * (1 - s_lo_a) / aY ** 2
    b0_p = binf_p / _cproto.NP
    kap0 = 1 / ((1 + s0) * aY ** 2)
    dk = amax(abs(1 / ((1 + s_lo_a) * aY ** 2) - kap0),
              abs(1 / ((1 + s_hi_a) * aY ** 2) - kap0))
    d2 = rho_t / arb(th2)
    om = 12 * d2 * ref["V_S"]
    c_Q = 1 - (om + nu * lam_p) * ref["S_Q"]
    c_S = 1 - (om + nu * lam_p) * ref["S_S"]
    lam_t = nu * lam_p * (1 + ref["S_M"]) + om * ref["S_M"]
    bt_inf = nu * binf_p + 4 * d2 * dk ** 2
    bt_0 = nu * b0_p + 4 * d2 * dk ** 2
    # Lemma D0: sigma_h = sqrt(#copies) * alpha_h^ref = sqrt(index)*...
    ncopy = getattr(_cproto, "INDEX", _cproto.NP + 1)
    sig = arb(ncopy).sqrt() * ref["alpha_ref"]
    carry = 1 - (1 / arb(th4) - 1) * ref["rho_w"]
    c_e = c_Q * carry - rho_t * (1 / arb(th) - 1) * sig ** 2
    d_e = lam_t * (1 + 1 / arb(al)) * ref["gamma"] ** 2
    return dict(kap0=kap0, c_Q=c_Q, c_S=c_S, lam_t=lam_t,
                bt_inf=bt_inf, bt_0=bt_0, c_e=c_e, d_e=d_e, carry=carry)


def build_A(glob, cQ_lo, cR_lo, rt1_lo, lt_up, btinf_up, bt0_up, k0):
    """Dense A_hat (safe midpoint directions) + eps_tot, memory-careful:
    single Fortran-ordered array (in-place LAPACK potrf, no copy),
    filled blockwise from sparse, rank-ones added in blocks."""
    ng = glob["ng"]
    Scomb = (cQ_lo * glob["Qm"] + cR_lo * glob["Rm"]
             - lt_up * glob["Mm"]).tocsr()
    A = np.zeros((ng, ng), order="F")
    for i0 in range(0, ng, 2048):
        i1 = min(i0 + 2048, ng)
        A[i0:i1, :] = Scomb[i0:i1].toarray()
    del Scomb
    ellm = glob["am"] + k0 * glob["tim"] + k0 * glob["t0m"]
    ellr = (glob["ar"] + abs(k0) * (glob["tir"] + glob["t0r"])
            + 1e-12 * np.abs(ellm))
    for i0 in range(0, ng, 2048):
        i1 = min(i0 + 2048, ng)
        A[i0:i1] += rt1_lo * ellm[i0:i1, None] * ellm[None, :]
        A[i0:i1] -= btinf_up * glob["tim"][i0:i1, None] * glob["tim"][None, :]
        A[i0:i1] -= bt0_up * glob["t0m"][i0:i1, None] * glob["t0m"][None, :]
    # eps_tot: matrix radii + Lemma R rank-one shifts + build error
    nrm = lambda v: np.linalg.norm(v) * (1 + 1e-12)
    rows_abs = np.abs(glob["Qr"]).sum(1).max() * cQ_lo \
        + np.abs(glob["Rr"]).sum(1).max() * cR_lo \
        + np.abs(glob["Mr"]).sum(1).max() * lt_up
    eps = float(rows_abs)
    eps += rt1_lo * 2 * nrm(ellm) * nrm(ellr)
    for tm, tr, bu in ((glob["tim"], glob["tir"], btinf_up),
                       (glob["t0m"], glob["t0r"], bt0_up)):
        eps += bu * (2 * nrm(tm) * nrm(tr) + nrm(tr) ** 2)
    # float build error (each entry <= ~10 flops); row-sum bound w/o
    # densifying: |A| parts bounded analytically
    absrow = (np.abs(cQ_lo * glob["Qm"]) + np.abs(cR_lo * glob["Rm"])
              + np.abs(lt_up * glob["Mm"])).sum(1).max()
    absrow += rt1_lo * np.abs(ellm).max() * np.abs(ellm).sum()
    absrow += btinf_up * np.abs(glob["tim"]).max() * np.abs(glob["tim"]).sum()
    absrow += bt0_up * np.abs(glob["t0m"]).max() * np.abs(glob["t0m"]).sum()
    eps += 10 * EPS * float(absrow)
    return A, eps


def power_of_two_diag_scale(diag, tiny=None):
    """s_i = 2^round(-1/2 log2(a_ii)). Powers of two ⇒ exact in IEEE-754."""
    if tiny is None:
        tiny = np.finfo(float).tiny
    d = np.maximum(np.asarray(diag, dtype=float), tiny)
    e = np.rint(-0.5 * np.log2(d)).astype(np.int32)
    return np.ldexp(np.ones(d.shape[0], dtype=float), e)


def _abs_csr(M):
    out = M.tocsr(copy=True)
    out.data = np.abs(out.data)
    return out


def scaled_radius_rows(glob, cQ_lo, cR_lo, rt1_lo, lt_up, btinf_up, bt0_up,
                       k0m, k0r, s):
    """Per-row radius sums r_i ≥ Σ_j s_i |ΔA|_ij s_j (scaled radius model).

    Sufficient for Rump: subtract r_i from diagonal i (not max_i r_i from all).
    Matrix part exact for |R|@s ⊙ s.  Rank-one |Δ(uuᵀ)| row-sum:
      r_i ≥ |u_i|‖du‖₁ + |du_i|‖u‖₁ + |du_i|‖du‖₁  (and same for −bt t tᵀ).
    """
    s = np.asarray(s, dtype=float)
    n = s.shape[0]
    Qr, Rr, Mr = (_abs_csr(glob["Qr"]), _abs_csr(glob["Rr"]),
                  _abs_csr(glob["Mr"]))
    # s_i * (|R| s)_i
    r = s * (cQ_lo * (Qr @ s) + cR_lo * (Rr @ s) + lt_up * (Mr @ s))

    ellm = glob["am"] + k0m * (glob["tim"] + glob["t0m"])
    ellr = (glob["ar"] + abs(k0m) * (glob["tir"] + glob["t0r"])
            + abs(k0r) * np.abs(glob["tim"] + glob["t0m"])
            + 1e-12 * np.abs(ellm))
    um, ur = s * ellm, s * np.abs(ellr)
    # outer-product radius row sums, times rt1
    u1, du1 = np.sum(np.abs(um)), np.sum(ur)
    r = r + rt1_lo * (np.abs(um) * du1 + ur * u1 + ur * du1)

    def add_tt(tm, tr, bu):
        nonlocal r
        umt, urt = s * tm, s * np.abs(tr)
        t1, dt1 = np.sum(np.abs(umt)), np.sum(urt)
        r = r + bu * (np.abs(umt) * dt1 + urt * t1 + urt * dt1)

    add_tt(glob["tim"], glob["tir"], btinf_up)
    add_tt(glob["t0m"], glob["t0r"], bt0_up)

    # float build padding: 10 EPS * scaled |mid| row-sums
    Qm, Rm, Mm = (_abs_csr(glob["Qm"]), _abs_csr(glob["Rm"]),
                  _abs_csr(glob["Mm"]))
    absrow = s * (cQ_lo * (Qm @ s) + cR_lo * (Rm @ s) + lt_up * (Mm @ s))
    absrow = absrow + rt1_lo * np.abs(um) * np.sum(np.abs(um))
    absrow = absrow + btinf_up * np.abs(s * glob["tim"]) * np.sum(
        np.abs(s * glob["tim"]))
    absrow = absrow + bt0_up * np.abs(s * glob["t0m"]) * np.sum(
        np.abs(s * glob["t0m"]))
    r = r + 10 * EPS * absrow
    return r * (1 + 4 * EPS)


def rump_certify_inplace(A, extra, verbose=False, radius_rows_for_s=None):
    """Rump BIT 46 (2006) Thm 2.3 + Cor 2.4 + Lemma 2.5; destroys A.

    Exact power-of-two diagonal congruence A ↦ S A S, then per-row radius
    diagonal reduction when radius_rows_for_s is provided:
      d'_i = a_ii - c_Δ - r_i
    with r = radius_rows_for_s(s) ≥ row sums of |S ΔA S|.
    Fallback: uniform extra scaled by max(s)² (legacy).
    """
    n = A.shape[0]
    diag0 = np.diag(A).copy()
    if not np.all(diag0 > 0):
        if verbose:
            print(f"      rump: nonpositive diagonal "
                  f"(min={diag0.min():.3e})", flush=True)
        return False, 0.0
    s = power_of_two_diag_scale(diag0)
    s_max = float(s.max())
    # In-place SAS
    for i0 in range(0, n, 2048):
        i1 = min(i0 + 2048, n)
        A[i0:i1, :] *= s[i0:i1, None]
        A[i0:i1, :] *= s[None, :]
    g = (n + 1) * EPS / (1.0 - (n + 1) * EPS)
    diag = np.diag(A).copy()
    d2 = np.maximum(diag, 0.0) / (1.0 - g)
    M = 3.0 * (2 * n + max(float(diag.max()), 0.0))
    # pure rounding part of Rump (shared); radius handled per-row below
    c_round = (g * d2.sum() * (1 + n * EPS) + n * M * ETA) * (1 + 8 * EPS)
    if radius_rows_for_s is not None:
        r = np.asarray(radius_rows_for_s(s), dtype=float)
        assert r.shape == (n,)
        c_vec = (c_round + r) * (1 + 8 * EPS)
        c_report = float(c_vec.max())
    else:
        extra_S = float(extra) * (s_max * s_max) * (1 + 4 * EPS)
        c_vec = np.full(n, (c_round + extra_S) * (1 + 8 * EPS))
        c_report = float(c_vec[0])
    phi = EPS * (1 + 2 * EPS)
    dprime = diag - c_vec
    if np.any(dprime <= 0):
        if verbose:
            print(f"      scale s∈[{s.min():.3e},{s_max:.3e}] "
                  f"c_max={c_report:.3e} r_max="
                  f"{float(np.max(c_vec - c_round)):.3e} "
                  f"diag' min={dprime.min():.3e} (nonpos after shift)",
                  flush=True)
        return False, c_report
    np.fill_diagonal(A, dprime - phi * np.abs(dprime))
    if verbose:
        rmax = float(np.max(c_vec - c_round))
        print(f"      scale s∈[{s.min():.3e},{s_max:.3e}] "
              f"c_max={c_report:.3e} r_max={rmax:.3e} "
              f"diag∈[{diag.min():.3e},{diag.max():.3e}] "
              f"diag'∈[{dprime.min():.3e},{dprime.max():.3e}]",
              flush=True)
    try:
        cholesky(A, lower=False, overwrite_a=True, check_finite=False)
        return True, c_report
    except LinAlgError:
        return False, c_report


def certify(N1, N2, N3, params=None, level=None):
    if level is not None:
        set_level(level)
    print(f"=== M-p2: Gamma_0{_cproto.LEVEL} N(p)={_cproto.NP}, "
          f"reference {N1}x{N2}x{N3} (24-split), Y={Y}, nu*={NU_STAR}, "
          f"kappa=Lemma I1 ===")
    ref = reference_arb(N1, N2, N3)
    print(f"  reference constants (arb upper): gamma={upper(ref['gamma']):.4f}"
          f" tau={upper(ref['tau']):.4f} alpha_ref={upper(ref['alpha_ref']):.4f}")
    print(f"    S_M={upper(ref['S_M']):.3f} S_Q={upper(ref['S_Q']):.2e} "
          f"S_S={upper(ref['S_S']):.2e} V_S={upper(ref['V_S']):.2e}")
    glob = assemble_global(ref)
    sgrid = np.linspace(0.0, 1.0, NWIN + 1)

    # ---- float parameter search: scalar conditions + LOBPCG pencil ----
    def pencil_mineig(th, th2, al, th4, rt, w):
        """Float smallest eig of (N_h, D_h) for window w via sparse LOBPCG."""
        from scipy.sparse.linalg import splu, lobpcg, LinearOperator
        co = window_coeffs(ref, sgrid[w], sgrid[w + 1], th, th2, al,
                           arb(rt), arb(NU_STAR), th4)
        cQ = float(co["c_Q"].mid())
        k0 = float(co["kap0"].mid())
        ell = glob["am"] + k0 * (glob["tim"] + glob["t0m"])
        rt1 = rt * (1 - th)
        lt = float(co["lam_t"].mid()) * (1 + al)
        bi = float(co["bt_inf"].mid())
        b0 = float(co["bt_0"].mid())
        Nsp = (cQ * (glob["Qm"] + (1 - th4) * glob["Rm"])).tocsr()
        Dsp = (lt * glob["Mm"]).tocsr()
        ng = glob["ng"]

        def nmv(x):
            x = np.asarray(x).ravel()
            return Nsp @ x + rt1 * (ell @ x) * ell

        def dmv(x):
            x = np.asarray(x).ravel()
            return (Dsp @ x + bi * (glob["tim"] @ x) * glob["tim"]
                    + b0 * (glob["t0m"] @ x) * glob["t0m"])

        Nop = LinearOperator((ng, ng), matvec=nmv, dtype=float)
        Dop = LinearOperator((ng, ng), matvec=dmv, dtype=float)
        prec = splu((Nsp + Dsp).tocsc())
        Mprec = LinearOperator((ng, ng), matvec=prec.solve, dtype=float)
        rng = np.random.default_rng(1)
        X = rng.standard_normal((ng, 4))
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            vals, _ = lobpcg(Nop, X, B=Dop, M=Mprec, largest=False,
                             tol=1e-5, maxiter=200)
        return np.sort(vals)[0]

    if params is None:
        volw = float(np.ones(glob["ng"]) @ (glob["Mm"] @ np.ones(glob["ng"])))
        t1, t2 = glob["tim"].sum(), glob["t0m"].sum()
        a1 = glob["am"].sum()
        cands = []
        # Wider grid than the Np=5 freeze: rho~ demand shrinks with N(p),
        # but sigma_h^2 = (N(p)+1) alpha^2 grows, so c_e pressure rises.
        # Scalar slack anti-correlates with pencil via rho~ — do not only
        # trial the highest-slack candidates (they prefer large rho~).
        for th in (0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.85):
            for th2 in (0.6, 0.7, 0.8, 0.9):
                for al in (0.15, 0.2, 0.25, 0.35, 0.5):
                    for th4 in (0.8, 0.85, 0.9):
                        for rt in (1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0,
                                   6.0, 7.5, 9.0):
                            slack = np.inf
                            for w in range(NWIN):
                                co = window_coeffs(
                                    ref, sgrid[w], sgrid[w + 1],
                                    th, th2, al, arb(rt),
                                    arb(NU_STAR), th4)
                                ce = float(co["c_e"].mid()
                                           - co["c_e"].rad())
                                de = upper(co["d_e"])
                                if float(co["c_S"].mid()) <= 0 or ce <= 0:
                                    slack = -np.inf
                                    break
                                slack = min(slack, ce / de - 1)
                                # constant-direction heuristic
                                k0 = float(co["kap0"].mid())
                                L1 = a1 + k0 * (t1 + t2)
                                D1 = (float(co["lam_t"].mid())
                                      * (1 + al) * volw
                                      + float(co["bt_inf"].mid()) * t1 ** 2
                                      + float(co["bt_0"].mid()) * t2 ** 2)
                                slack = min(slack, rt * (1 - th)
                                            * L1 ** 2 / D1 - 1)
                            if slack > 0:
                                cands.append((slack, th, th2, al, th4, rt))
        cands.sort(reverse=True, key=lambda c: c[0])
        print(f"  scalar-feasible candidates: {len(cands)}", flush=True)
        assert cands, "no scalar-feasible parameters at this mesh"
        # Trial order: low rho~ first (pencil-friendly), then high slack
        trial = sorted(cands, key=lambda c: (c[5], -c[0]))
        best = None
        n_try = min(40, len(trial))
        for c in trial[:n_try]:
            ev = min(pencil_mineig(*c[1:], w) for w in (0, NWIN - 1))
            print(f"    cand {c[1:]}: scalar slack {c[0]:.3f}, "
                  f"pencil min-eig {ev:.4f}", flush=True)
            score = min(c[0], ev - 1)
            if best is None or score > best[0]:
                best = (score, *c[1:])
            if score > 0.15:
                break
        print(f"  chosen: score={best[0]:.3f} params={best[1:]}", flush=True)
        assert best[0] > 0, "no feasible parameters at this mesh"
        params = best[1:]
    th, th2, al, th4, rt = params

    # ---- rigorous per-window certification ----
    all_ok = True
    print(f"  {'window':>16} {'c_S>0':>6} {'c_e>=d_e':>9} {'PSD':>5} "
          f"{'c_e':>7} {'d_e':>7}", flush=True)
    for w in range(NWIN):
        co = window_coeffs(ref, sgrid[w], sgrid[w + 1], th, th2, al,
                           arb(rt), arb(NU_STAR), th4)
        ok_s = bool(co["c_S"] > 0) and bool(co["c_e"] > co["d_e"]) \
            and bool(co["c_Q"] > 0) and bool(co["carry"] > 0)
        cQ_lo = float(co["c_Q"].mid() - co["c_Q"].rad())
        cR_lo = float((co["c_Q"] * (1 - arb(th4))).mid()
                      - (co["c_Q"] * (1 - arb(th4))).rad())
        rt1_lo = float((arb(rt) * (1 - arb(th))).mid()) * (1 - 1e-14)
        lt_up = upper(co["lam_t"] * (1 + arb(al)))
        btinf_up = upper(co["bt_inf"])
        bt0_up = upper(co["bt_0"])
        k0m, k0r = mid_rad(co["kap0"])
        A, eps = build_A(glob, cQ_lo, cR_lo, rt1_lo, lt_up, btinf_up,
                         bt0_up, k0m)
        # kap0 radius folded into ell radius inside build_A via 1e-12 guard;
        # add explicitly (unscaled fallback path only):
        eps += rt1_lo * 2 * np.linalg.norm(
            glob["am"] + k0m * (glob["tim"] + glob["t0m"])) * \
            (k0r * np.linalg.norm(glob["tim"] + glob["t0m"])) * 1.01

        def radius_rows_for_s(s):
            return scaled_radius_rows(
                glob, cQ_lo, cR_lo, rt1_lo, lt_up, btinf_up, bt0_up,
                k0m, k0r, s)

        psd, c = rump_certify_inplace(
            A, eps, verbose=True, radius_rows_for_s=radius_rows_for_s)
        del A
        ok = ok_s and psd
        all_ok &= ok
        print(f"  [{sgrid[w]:.3f},{sgrid[w+1]:.3f}] {str(bool(co['c_S']>0)):>6}"
              f" {str(bool(co['c_e']>co['d_e'])):>9} {str(psd):>5} "
              f"{float(co['c_e'].mid()):7.3f} {float(co['d_e'].mid()):7.3f}",
              flush=True)
    print(f"  ==> Gamma_0{_cproto.LEVEL} CERTIFIED (all windows): {all_ok}")
    if all_ok:
        print(f"      no eigenvalues in (0,1) for "
              f"Gamma_0{_cproto.LEVEL}\\H^3 —")
        print("      by the multi-cusp criterion (CONGRUENCE.md sec 2) +"
              " Theorem G1p + Rump BIT 46 (2006).")
    return all_ok, params


if __name__ == "__main__":
    import sys
    # usage: python m3p_certify.py coarse [(2+i)|(3)|(3+2i)]
    args = sys.argv[1:]
    mesh = args[0] if args else "coarse"
    lev = args[1] if len(args) > 1 else "(2+i)"
    if mesh == "coarse":
        certify(6, 3, 3, level=lev)
    else:
        certify(8, 4, 3)
