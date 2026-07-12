"""M0 float prototype: independent-exclusion criterion for
PSL(2, Z[omega]) level 1 (Eisenstein–Picard).

Same architecture as independent_exclusion/fem_prototype.py (Picard),
with geometry adapted to the hexagonal cusp lattice:

  Lambda = Z[omega],  omega = -1/2 + i sqrt(3)/2
  parallelogram slice P: z = u + v*omega,  u in [0,1], v in [0,1/3]
  |T| = area(P) = sqrt(3)/6   (= covol(Lambda)/3, matching [Ginf:G'inf]=3)
  core K: y_f(z) <= y <= Y,  y_f = sqrt(1-|z|^2)

Criterion (DESIGN.md, |T|-parametrized):
  A_s = Q_K - lam M_K - beta t^2,  beta = (1-s)/(|T| Y^2)
  L_s = a + t / ((1+s) Y^2)
  need min A_s on {L_s=0} > 0 for all lam in (0,1).

Conforming FEM => upper bound on true mu (evidence / margin, not proof).
Neumann relaxation: no side/floor identifications (same as Picard M0).

Mode-bound check: shortest dual of hexagonal lattice ~ 2/sqrt(3);
  4 pi^2 Y^2 / |mu_min|^2 >> 1 for Y>=1.25.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.linalg import eigh

Y_DEFAULT = 1.25
# Full lattice parallelogram has area sqrt(3)/2.  With [Gamma_inf:Gamma'_inf]=3
# the orbifold cusp section has area |T| = sqrt(3)/6 (translations by Lambda
# mod a order-3 unit action).  We mesh one fundamental parallelogram slice
# u in [0,1], v in [0, VMAX] with VMAX=1/3 so area = (sqrt3/2)*(1/3)=sqrt3/6.
# This R_comp is volume-/cusp-equivalent to the EGM truncated core K_Y; see
# geometry_fund.py and GEOMETRY.md.
VMAX = 1.0 / 3.0
AREA_T = (math.sqrt(3) / 2) * VMAX   # = sqrt(3)/6


def _vol_checks(Y):
    """EGM closed-form vol(K_Y) and quadrature vol(K_comp) for this mesh."""
    try:
        from geometry_fund import (
            vol_KY_exact, vol_F_exact, cusp_tail_volume, vol_K_comp_quad,
        )
        vK_egm, rad = vol_KY_exact(Y)
        vF, _ = vol_F_exact()
        vK_comp = vol_K_comp_quad(Y, n=500)
        return dict(
            vK_egm=vK_egm, rad=rad, vF=vF,
            tail=cusp_tail_volume(Y), vK_comp=vK_comp,
        )
    except Exception:
        vF = 0.16915693440160895
        tail = AREA_T / (2.0 * Y * Y)
        return dict(
            vK_egm=vF - tail, rad=1e-12, vF=vF, tail=tail, vK_comp=None,
        )


def z_from_uv(u, v):
    """z = u + v*omega = (u - v/2) + i v sqrt(3)/2."""
    return u - 0.5 * v, (math.sqrt(3) / 2) * v


def y_floor(x1, x2):
    r2 = x1 ** 2 + x2 ** 2
    return np.sqrt(np.maximum(1.0 - r2, 0.0))


def build_matrices(N1, N2, N3, Y):
    """Q1 FEM on parameter box (u,v,tau): u,tau in [0,1], v in [0,VMAX]."""
    ug = np.linspace(0.0, 1.0, N1 + 1)
    vg = np.linspace(0.0, VMAX, N2 + 1)
    tg = np.linspace(0.0, 1.0, N3 + 1)
    hu, hv, ht = 1.0 / N1, VMAX / N2, 1.0 / N3
    nn = (N1 + 1) * (N2 + 1) * (N3 + 1)

    def nid(i, j, k):
        return (i * (N2 + 1) + j) * (N3 + 1) + k

    corners = np.array([[a, b, c] for a in (0, 1) for b in (0, 1)
                        for c in (0, 1)])
    gp1 = np.array([0.5 - 0.5 / np.sqrt(3), 0.5 + 0.5 / np.sqrt(3)])
    gpts = np.array([[p, q, r] for p in gp1 for q in gp1 for r in gp1])
    wgt = np.full(8, 0.125)

    def shp(xi):
        G = xi.shape[0]
        N = np.ones((G, 8))
        dN = np.ones((G, 8, 3))
        for d in range(3):
            f = np.where(corners[None, :, d] == 1, xi[:, None, d],
                         1.0 - xi[:, None, d])
            df = np.where(corners[None, :, d] == 1, 1.0, -1.0)
            N *= f
            for dd in range(3):
                dN[:, :, dd] *= f if dd != d else df
        return N, dN

    Nsh, dNsh = shp(gpts)
    S = np.zeros((nn, nn))
    M = np.zeros((nn, nn))
    avec = np.zeros(nn)
    tvec = np.zeros(nn)

    els = [(i, j, k) for i in range(N1) for j in range(N2) for k in range(N3)]
    conn = np.array([[nid(i + c[0], j + c[1], k + c[2]) for c in corners]
                     for (i, j, k) in els])
    i_arr = np.array([e[0] for e in els])
    j_arr = np.array([e[1] for e in els])
    k_arr = np.array([e[2] for e in els])

    # gauss points in (u,v,tau)
    u = ug[i_arr][:, None] + gpts[None, :, 0] * hu
    v = vg[j_arr][:, None] + gpts[None, :, 1] * hv
    tt = tg[k_arr][:, None] + gpts[None, :, 2] * ht
    x1 = u - 0.5 * v
    x2 = (math.sqrt(3) / 2) * v
    yf = y_floor(x1, x2)
    y = yf * (1.0 - tt) + Y * tt

    # Jacobian: (x1,x2,y) w.r.t (xi_u, xi_v, xi_t) where
    # u = u0 + xi_u*hu, etc.
    # dx1/du=1, dx1/dv=-1/2, dx2/du=0, dx2/dv=sqrt3/2
    # dy/du = (1-tt)*dyf/dx1*dx1/du + ...; chain through yf(x1,x2)
    r2 = x1 ** 2 + x2 ** 2
    # dyf/dx1 = -x1/yf, dyf/dx2 = -x2/yf when yf>0
    safe = np.maximum(yf, 1e-14)
    dyf_dx1 = np.where(yf > 1e-12, -x1 / safe, 0.0)
    dyf_dx2 = np.where(yf > 1e-12, -x2 / safe, 0.0)
    dyf_du = dyf_dx1 * 1.0 + dyf_dx2 * 0.0
    dyf_dv = dyf_dx1 * (-0.5) + dyf_dx2 * (math.sqrt(3) / 2)

    # y = yf*(1-tt)+Y*tt
    # partial y / partial xi_u = dy/du * hu
    # partial y / partial xi_v = dy/dv * hv
    # partial y / partial xi_t = (Y-yf)*ht
    dy_du = (1.0 - tt) * dyf_du
    dy_dv = (1.0 - tt) * dyf_dv
    # Euclidean (x1,x2) Jacobian times y column:
    # |det d(x1,x2)/d(u,v)| = sqrt3/2
    # Full 3D det for (x1,x2,y)/(xi_u,xi_v,xi_t):
    # det = (sqrt3/2)*hu*hv * (Y-yf)*ht
    detJ = (math.sqrt(3) / 2) * hu * hv * (Y - yf) * ht
    if detJ.min() <= 0:
        raise RuntimeError(
            f"degenerate elements: min detJ={detJ.min()} "
            f"(need Y>1 and |z|<=1 on mesh)")

    # Gradients in physical (x1,x2,y): invert Jacobian of map
    # J columns: dX/dxi_u, dX/dxi_v, dX/dxi_t
    # d(x1,x2)/d(u,v) * (hu,hv) for first two columns in x1,x2
    # For CR we need dN/dx. Use chain rule via inverse of
    # [[dx1/du, dx1/dv, dx1/dtt],
    #  [dx2/du, dx2/dv, dx2/dtt],
    #  [dy/du,  dy/dv,  dy/dtt ]] with dtt scaled by ht, etc.
    #
    # Parameter xi: (u,v,tt) = corner + xi * h
    # d(x1)/dxi_u = dx1/du * hu = hu
    # d(x1)/dxi_v = dx1/dv * hv = -0.5*hv
    # d(x1)/dxi_t = 0
    # d(x2)/dxi_u = 0
    # d(x2)/dxi_v = (sqrt3/2)*hv
    # d(x2)/dxi_t = 0
    # d(y)/dxi_u = dy_du * hu
    # d(y)/dxi_v = dy_dv * hv
    # d(y)/dxi_t = (Y-yf)*ht

    # Invert 3x3 for each gauss point is expensive in pure python loops;
    # closed form: x1,x2 independent of tt.
    # Let A = d(x1,x2)/d(u,v) = [[1,-1/2],[0,sqrt3/2]], detA=sqrt3/2
    # dN/dx1, dN/dx2 from (u,v) only for horizontal, then correct for y.
    # Standard FEM on product: 
    #   grad_{x1,x2} N = A^{-T} grad_{u,v} N
    #   dN/dy = dN/dtt / ((Y-yf)*ht) with tt chain
    #
    # A^{-1} = (2/sqrt3) [[sqrt3/2, 1/2],[0, 1]] = [[1, 1/sqrt3],[0, 2/sqrt3]]
    # A^{-T} = [[1, 0],[1/sqrt3, 2/sqrt3]]

    dNu = dNsh[None, :, :, 0] / hu   # dN/du
    dNv = dNsh[None, :, :, 1] / hv
    dNt = dNsh[None, :, :, 2] / ht   # dN/dtt

    # dN/dx1, dN/dx2 ignoring y-tilt of iso-u,v surfaces:
    # more carefully, surfaces of const (u,v) vary y with u,v.
    # Full inverse of J:
    # J = [[hu, -0.5 hv, 0],
    #      [0, (sqrt3/2)hv, 0],
    #      [dy_du*hu, dy_dv*hv, (Y-yf)*ht]]
    # det J = hu * (sqrt3/2 hv) * (Y-yf)ht = detJ (matches)
    # J^{-1} via block form:
    # top-left 2x2 inv of A_h = [[hu,-0.5 hv],[0, s3/2 hv]]
    # A_h^{-1} = [[1/hu, 1/(sqrt3 hu)],[0, 2/(sqrt3 hv)]]
    inv00 = 1.0 / hu
    inv01 = 1.0 / (math.sqrt(3) * hu)
    inv11 = 2.0 / (math.sqrt(3) * hv)
    # dN/dx1 = inv00 * dN/dxi_u + 0 * dN/dxi_v + ...
    # Actually dN/d(x1,x2,y) = J^{-T} dN/dxi
    # For lower triangular block structure:
    # Let p = dN/dxi (3-vector)
    # Solve J^T g = p for g = grad_phys N
    #
    # J^T = [[hu, 0, dy_du*hu],
    #        [-0.5 hv, (s3/2)hv, dy_dv*hv],
    #        [0, 0, (Y-yf)*ht]]
    #
    # g_y = p_t / ((Y-yf)*ht) = dNt / (Y-yf)  but p_t = dN/dxi_t = dNsh_t
    # wait dN/dxi_t = dNsh[:,:,2], p_t = that, g_y = p_t / ((Y-yf)*ht) = dNt/(Y-yf)? 
    # dNt defined as dNsh/ht so dN/dtt = dNt, and dN/dxi_t = dN/dtt * ht = dNsh_2
    # g_y = dNsh_2 / ((Y-yf)*ht) = dNt / (Y-yf)

    Gy = dNsh[None, :, :, 2] / ((Y - yf) * ht)[:, :, None]

    # back-sub for g_x1, g_x2 from first two rows of J^T g = p
    # hu * g_x1 + (dy_du*hu) * g_y = p_u = dNsh_0
    # -0.5 hv g_x1 + (s3/2)hv g_x2 + dy_dv*hv g_y = p_v = dNsh_1
    p_u = dNsh[None, :, :, 0]   # (ne,G,8)
    p_v = dNsh[None, :, :, 1]
    # hu * Gx1 + (dy_du*hu) * Gy = p_u  =>  Gx1 = p_u/hu - dy_du * Gy
    Gx1 = p_u / hu - dy_du[:, :, None] * Gy
    # -0.5 hv Gx1 + (√3/2) hv Gx2 + dy_dv*hv*Gy = p_v
    Gx2 = (p_v + 0.5 * hv * Gx1 - dy_dv[:, :, None] * (hv * Gy)) / (
        (math.sqrt(3) / 2) * hv)

    wS = wgt[None, :] * detJ / y
    wM = wgt[None, :] * detJ / y ** 3
    Se = (np.einsum("eg,ega,egb->eab", wS, Gx1, Gx1)
          + np.einsum("eg,ega,egb->eab", wS, Gx2, Gx2)
          + np.einsum("eg,ega,egb->eab", wS, Gy, Gy))
    Me = np.einsum("eg,ga,gb->eab", wM, Nsh, Nsh)
    ae = np.einsum("eg,ga->ea", wM, Nsh)

    rows = np.repeat(conn, 8, axis=1).ravel()
    cols = np.tile(conn, (1, 8)).ravel()
    np.add.at(S, (rows, cols), Se.ravel())
    np.add.at(M, (rows, cols), Me.ravel())
    np.add.at(avec, conn.ravel(), ae.ravel())

    # top face τ=1: t = ∫_P v dx1 dx2 = ∫ v (√3/2) du dv
    # Q1: each top corner of top-layer element gets (√3/2)*hu*hv/4
    top_w = (math.sqrt(3) / 2) * hu * hv / 4.0
    top_corner_ids = [c_idx for c_idx, c in enumerate(corners) if c[2] == 1]
    for (i, j, k), cn in zip(els, conn):
        if k == N3 - 1:
            for c_idx in top_corner_ids:
                tvec[cn[c_idx]] += top_w

    return S, M, avec, tvec


def constrained_min_eig(A, M, ell, k=2):
    n = A.shape[0]
    v = ell.copy()
    v[0] += np.sign(ell[0] if ell[0] != 0 else 1.0) * np.linalg.norm(ell)
    beta = 2.0 / np.dot(v, v)

    def hah(B):
        w = B @ v
        B2 = (B - beta * np.outer(v, w) - beta * np.outer(w, v)
              + (beta ** 2 * np.dot(v, w)) * np.outer(v, v))
        return B2[1:, 1:]

    return eigh(hah(A), hah(M), eigvals_only=True, subset_by_index=[0, k - 1])


def mode_bound_ok(Y):
    # dual of hexagonal lattice: min |mu| = 2/sqrt(3)
    mu_min = 2.0 / math.sqrt(3)
    return 4 * math.pi ** 2 * Y ** 2 / mu_min ** 2


def run(N1, N2, N3, Y, lams, verbose=True):
    if verbose:
        print(f"--- Eisenstein–Picard FEM M0: mesh {N1}x{N2}x{N3}, Y={Y} ---")
        print(f"  |T|=sqrt(3)/6={AREA_T:.6f}  mode-bound 4pi^2 Y^2/|mu|^2 "
              f"= {mode_bound_ok(Y):.2f} (need >1)")
    S, M, avec, tvec = build_matrices(N1, N2, N3, Y)
    one = np.ones(S.shape[0])
    volK_fem = one @ (M @ one)
    t1 = tvec @ one
    # Volume self-check:
    #   EGM core:  vol(K_Y) = vol(F) − |T|/(2 Y²)   (paper domain)
    #   this mesh: vol(K_comp) = ∫_{R_comp} ½(yf^{-2}−Y^{-2})  (FEM mass target)
    vols = _vol_checks(Y)
    if verbose:
        print(f"  EGM vol(K_Y)={vols['vK_egm']:.8f}  "
              f"(vol(F)={vols['vF']:.8f} − tail={vols['tail']:.8f})")
        if vols["vK_comp"] is not None:
            rel = abs(volK_fem - vols["vK_comp"]) / max(vols["vK_comp"], 1e-30)
            print(f"  check vol(K_comp): fem {volK_fem:.8f}  "
                  f"quad {vols['vK_comp']:.8f}  (rel err {rel:.2e})")
        else:
            print(f"  check vol(K): fem {volK_fem:.8f}  (no quad available)")
        print(f"  check t(1):   fem {t1:.8f}   exact |T|={AREA_T:.8f}")
        print(f"  check |S@1|:  {np.abs(S @ one).max():.2e}")
        neu = eigh(S, M, eigvals_only=True, subset_by_index=[0, 2])
        print(f"  Neumann eigs of core: {neu}")

    # beta = (1-s)/(|T| Y^2); Picard special case |T|=1/2 => 2(1-s)/Y^2
    print(f"  {'lambda':>8} {'s':>7} {'mu1':>10} {'mu2':>10}")
    mus = []
    for lam in lams:
        s = math.sqrt(1.0 - lam)
        beta = (1.0 - s) / (AREA_T * Y ** 2)
        A = S - lam * M - beta * np.outer(tvec, tvec)
        ell = avec + tvec / ((1.0 + s) * Y ** 2)
        m = constrained_min_eig(A, M, ell)
        mus.append(m[0])
        print(f"  {lam:8.4f} {s:7.4f} {m[0]:10.5f} {m[1]:10.5f}")
    mmin = min(mus)
    print(f"  margin: min mu = {mmin:.5f}  at lambda = "
          f"{lams[int(np.argmin(mus))]}")
    return mmin


if __name__ == "__main__":
    lams = np.array([0.05, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 0.999])
    print("Independent exclusion M0 — PSL(2,Z[omega]) level 1")
    print("=" * 60)
    # modest mesh first (parallelogram [0,1]^2 denser than Picard half-rect)
    m1 = run(8, 8, 4, 1.25, lams)
    print()
    m2 = run(12, 12, 6, 1.25, lams)
    print()
    print("=" * 60)
    print(f"SUMMARY min mu @ 8^2x4: {m1:.4f}; @ 12^2x6: {m2:.4f}")
    if m2 > 1.5:
        print("GO: float margin healthy for Neumann relaxation (like Picard M0)")
    elif m2 > 0:
        print("MARGINAL: positive but thin — tighten mesh/Y before CR cert")
    else:
        print("NO-GO: nonpositive discrete mu — criterion or geometry bug")
