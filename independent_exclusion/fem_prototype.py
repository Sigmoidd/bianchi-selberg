"""FEM prototype for the Target-B exclusion criterion (see DESIGN.md).

For each lambda = 1 - s^2 in (0,1), computes

    mu(lambda) = min { A_s(v) / M_K(v) : v in V_h, L_s(v) = 0 }

where A_s(v) = Q_K(v) - lambda*M_K(v) - 2(1-s) t(v)^2 / Y^2 and
L_s(v) = a(v) + t(v)/((1+s) Y^2), on the compact core
K = { (x1,x2,y) : x1 in [-1/2,1/2], x2 in [0,1/2], sqrt(1-|z|^2) <= y <= Y }.

Criterion needs mu(lambda) > 0 for all lambda in (0,1). Conforming FEM
bounds the true infimum from ABOVE, so a positive mu here is evidence and
margin measurement, not proof (see DESIGN.md section 4).
"""

import numpy as np
from scipy.linalg import eigh

VOL_F = 0.30532186472  # 2 zeta_K(2)/pi^2, Picard fundamental domain (check value)


def build_matrices(N1, N2, N3, Y):
    """Assemble dense S (weighted stiffness), M (weighted mass),
    a (volume functional), t (top-face trace functional)."""
    x1g = np.linspace(-0.5, 0.5, N1 + 1)
    x2g = np.linspace(0.0, 0.5, N2 + 1)
    tg = np.linspace(0.0, 1.0, N3 + 1)
    h1, h2, h3 = 1.0 / N1, 0.5 / N2, 1.0 / N3
    nn = (N1 + 1) * (N2 + 1) * (N3 + 1)

    def nid(i, j, k):
        return (i * (N2 + 1) + j) * (N3 + 1) + k

    # reference Q1: corner a=(a1,a2,a3) in {0,1}^3
    corners = np.array([[a1, a2, a3] for a1 in (0, 1) for a2 in (0, 1)
                        for a3 in (0, 1)])          # (8,3)
    gp1 = np.array([0.5 - 0.5 / np.sqrt(3), 0.5 + 0.5 / np.sqrt(3)])
    gpts = np.array([[u, v, w] for u in gp1 for v in gp1 for w in gp1])  # (8,3)
    wgt = np.full(8, 0.125)

    def shp(xi):  # xi (G,3) -> N (G,8), dN (G,8,3)
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

    Nsh, dNsh = shp(gpts)  # (8,8), (8,8,3)

    S = np.zeros((nn, nn))
    M = np.zeros((nn, nn))
    avec = np.zeros(nn)
    tvec = np.zeros(nn)

    # element list
    els = [(i, j, k) for i in range(N1) for j in range(N2) for k in range(N3)]
    conn = np.array([[nid(i + c[0], j + c[1], k + c[2]) for c in corners]
                     for (i, j, k) in els])         # (ne,8)
    ne = len(els)

    i_arr = np.array([e[0] for e in els])
    j_arr = np.array([e[1] for e in els])
    k_arr = np.array([e[2] for e in els])

    # gauss-point physical data, shape (ne, G)
    x1 = x1g[i_arr][:, None] + gpts[None, :, 0] * h1
    x2 = x2g[j_arr][:, None] + gpts[None, :, 1] * h2
    tt = tg[k_arr][:, None] + gpts[None, :, 2] * h3
    yf = np.sqrt(1.0 - x1 ** 2 - x2 ** 2)
    y = yf * (1.0 - tt) + Y * tt

    dy_dx1 = (1.0 - tt) * (-x1 / yf)
    dy_dx2 = (1.0 - tt) * (-x2 / yf)
    Ja = dy_dx1 * h1          # dy/dxi1
    Jb = dy_dx2 * h2          # dy/dxi2
    Jc = (Y - yf) * h3        # dy/dxi3
    detJ = h1 * h2 * Jc       # (ne,G)
    assert detJ.min() > 0, "degenerate element (need Y>1)"

    # physical gradients of shape functions: (ne,G,8,3)
    dN0 = dNsh[None, :, :, 0]
    dN1 = dNsh[None, :, :, 1]
    dN2 = dNsh[None, :, :, 2]
    GX0 = dN0 / h1 - (Ja / (h1 * Jc))[:, :, None] * dN2
    GX1 = dN1 / h2 - (Jb / (h2 * Jc))[:, :, None] * dN2
    GX2 = dN2 / Jc[:, :, None]

    wS = (wgt[None, :] * detJ / y)          # (ne,G)
    wM = (wgt[None, :] * detJ / y ** 3)

    Se = (np.einsum('eg,ega,egb->eab', wS, GX0, GX0)
          + np.einsum('eg,ega,egb->eab', wS, GX1, GX1)
          + np.einsum('eg,ega,egb->eab', wS, GX2, GX2))
    Me = np.einsum('eg,ga,gb->eab', wM, Nsh, Nsh)
    ae = np.einsum('eg,ga->ea', wM, Nsh)

    rows = np.repeat(conn, 8, axis=1).ravel()
    cols = np.tile(conn, (1, 8)).ravel()
    np.add.at(S, (rows, cols), Se.ravel())
    np.add.at(M, (rows, cols), Me.ravel())
    np.add.at(avec, conn.ravel(), ae.ravel())

    # top-face trace functional t(v) = int_T v(.,Y) dx  (Euclidean),
    # exact for Q1: each top corner of a top-layer element gets h1*h2/4
    top_corner_ids = [c_idx for c_idx, c in enumerate(corners) if c[2] == 1]
    for (i, j, k), cn in zip(els, conn):
        if k == N3 - 1:
            for c_idx in top_corner_ids:
                tvec[cn[c_idx]] += h1 * h2 / 4.0

    return S, M, avec, tvec


def constrained_min_eig(A, M, ell, k=2):
    """Smallest k eigenvalues of x'Ax/x'Mx subject to ell'x = 0,
    via Householder reflection H mapping ell -> +-|ell| e_1."""
    n = A.shape[0]
    v = ell.copy()
    v[0] += np.sign(ell[0] if ell[0] != 0 else 1.0) * np.linalg.norm(ell)
    beta = 2.0 / np.dot(v, v)

    def hah(B):
        w = B @ v
        vBv = np.dot(v, w)
        B2 = B - beta * np.outer(v, w) - beta * np.outer(w, v) \
             + (beta ** 2 * vBv) * np.outer(v, v)
        return B2[1:, 1:]

    Am, Mm = hah(A), hah(M)
    vals = eigh(Am, Mm, eigvals_only=True, subset_by_index=[0, k - 1])
    return vals


def run(N1, N2, N3, Y, lams):
    print(f"--- mesh {N1}x{N2}x{N3}, Y={Y} ---")
    S, M, avec, tvec = build_matrices(N1, N2, N3, Y)
    one = np.ones(S.shape[0])

    volK = one @ M @ one
    volK_exact = VOL_F - 1.0 / (4.0 * Y * Y)
    print(f"check vol(K):  fem {volK:.6f}   exact {volK_exact:.6f}")
    print(f"check t(1):    fem {tvec @ one:.6f}   exact 0.500000")
    print(f"check |S@1|:   {np.abs(S @ one).max():.2e}")

    # context: pure Neumann spectrum of the core
    neu = eigh(S, M, eigvals_only=True, subset_by_index=[0, 2])
    print(f"Neumann eigs of core: {neu}")

    print(f"{'lambda':>8} {'s':>7} {'mu1':>10} {'mu2':>10}")
    mus = []
    for lam in lams:
        s = np.sqrt(1.0 - lam)
        A = S - lam * M - (2.0 * (1.0 - s) / Y ** 2) * np.outer(tvec, tvec)
        ell = avec + tvec / ((1.0 + s) * Y ** 2)
        m = constrained_min_eig(A, M, ell)
        mus.append(m[0])
        print(f"{lam:8.4f} {s:7.4f} {m[0]:10.5f} {m[1]:10.5f}")
    print(f"margin: min mu = {min(mus):.5f}  at lambda = {lams[np.argmin(mus)]}")
    return min(mus)


if __name__ == "__main__":
    lams = np.array([0.05, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 0.999])
    for Y in (1.25, 1.5):
        run(16, 8, 6, Y, lams)
        print()
    # refinement check at the chosen Y
    run(24, 12, 9, 1.25, lams)
