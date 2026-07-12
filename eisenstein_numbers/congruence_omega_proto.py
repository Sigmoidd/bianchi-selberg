"""Stage 8 float prototype — multi-copy CR for Γ₀(𝔭) ⊂ PSL(2, ℤ[ω]).

Reference-Cell Principle: NC = N(𝔭)+1 isometric copies of the level-1
P₃ core; gluing by residue right action (residue_omega) + geometric
face maps for generators T1, Tw, U, S.

Two-cusp pencil (cusps_omega.py):
  A = Q − λ M − β_∞ t_∞ t_∞ᵀ − β_0 t_0 t_0ᵀ
  on { a + κ_c (t_∞+t_0) = 0 }.

Status: **float evidence only** — not an interval certificate.
Self-identifications (perm[c]=c) are Neumann-relaxed (as Gaussian).

Face pairings: EGM P_3 edge dictionary in `face_pairings_p3.py`
  (RIGHT↔LEFT / LOW↔UP / U-rotate / S-floor).

Usage:
  python -u congruence_omega_proto.py           # N=3 default
  python -u congruence_omega_proto.py 7         # N=7
  python -u congruence_omega_proto.py 3 4 2     # N, N_tri, N3
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).resolve().parent))

from reference_cell import AREA_T, Y_DEFAULT, build_P3_cell
from residue_omega import (
    ring_lambda3, ring_split7, ring_split13, gluing_perms,
)
from cusps_omega import cusp_table_prime
from cr_omega import _orient_vol, TET_QP, TET_QW
from face_pairings_p3 import (
    build_pair_maps as discover_pair_maps,
    exterior_faces,
    validate_pair_maps,
)

Y = Y_DEFAULT
SQRT3 = math.sqrt(3.0)


def _ring_for_norm(q: int):
    if q == 3:
        return ring_lambda3()
    if q == 7:
        return ring_split7()
    if q == 13:
        return ring_split13()
    raise ValueError(f"unsupported N={q}; use 3,7,13")


def build_pair_maps(mesh, tf, nfr, bfaces=None, tol=None):
    """EGM face pairings via short words in T1,Tw,U,S (face_pairings_p3)."""
    pm, bf, meta = discover_pair_maps(mesh, tf, nfr, bfaces=bfaces, tol=tol)
    ok, msgs = validate_pair_maps(pm, bf)
    print(f"    pair dictionary tol={meta['tol']:.4f}  "
          f"method={meta.get('method')}  edge_hist={meta.get('edge_hist')}")
    print(f"    stats={meta['stats']}  validate={'OK' if ok else 'FAIL'}")
    for m in msgs[:6]:
        print(f"      {m}")
    for name, m in pm.items():
        print(f"    pair_map {name}: {len(m)} directed face matches")
    return pm


def ref_elements(mesh):
    """Per-tet Se, Me, ae and face table for one P₃ copy."""
    X, tets = mesh["X"], mesh["tets"]
    # orient
    tets = np.array(tets, dtype=int, copy=True)
    for e in range(len(tets)):
        P = X[list(tets[e])]
        A = np.hstack([np.ones((4, 1)), P])
        if np.linalg.det(A) < 0:
            tets[e, [0, 1]] = tets[e, [1, 0]]
    nt = len(tets)
    fid = {}
    tf = np.empty((nt, 4), dtype=int)
    for e in range(nt):
        for a in range(4):
            key = tuple(sorted(int(x) for x in np.delete(tets[e], a)))
            if key not in fid:
                fid[key] = len(fid)
            tf[e, a] = fid[key]
    nfr = len(fid)
    Se = np.zeros((nt, 4, 4))
    Me = np.zeros((nt, 4, 4))
    ae = np.zeros((nt, 4))
    Mloc = np.full((4, 4), -1 / 20.0) + np.eye(4) * (9 / 20.0)
    for e in range(nt):
        P = X[list(tets[e])]
        P, v, g = _orient_vol(P)
        gphi = -3.0 * g
        wq = wm = 0.0
        av = np.zeros(4)
        for q, qw in zip(TET_QP, TET_QW):
            lam = np.array([1 - q.sum(), *q])
            y = max(float(lam @ P[:, 2]), 1e-14)
            wq += qw / y
            wm += qw / y ** 3
            av += qw / y ** 3 * (1 - 3 * lam)
        Se[e] = wq * v * (gphi @ gphi.T)
        Me[e] = wm * v * Mloc
        ae[e] = v * av
    # top face areas
    top_pairs = []
    for tri in mesh["top_faces"]:
        key = tuple(sorted(int(x) for x in tri))
        if key not in fid:
            continue
        P = X[list(tri)]
        area = 0.5 * float(np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0])))
        top_pairs.append((fid[key], area))
    mesh_o = dict(mesh)
    mesh_o["tets"] = tets
    return dict(
        mesh=mesh_o, Se=Se, Me=Me, ae=ae, tf=tf, fid=fid, nfr=nfr,
        top_pairs=top_pairs, tets=tets,
    )


def assemble_level(q: int = 3, N_tri: int = 4, N3: int = 2, Y: float = Y):
    """NC copies of P₃ CR with gluing; return Q,M,a,t_inf,t_0 sparse/dense."""
    R = _ring_for_norm(q)
    ct = cusp_table_prime(q, R.label)
    NC = R.index
    glue = gluing_perms(R)
    pts = R.p1_points()
    cusp_class = [0 if p[0] == 0 else 1 for p in pts]
    assert cusp_class.count(0) == 1 and cusp_class.count(1) == q

    mesh = build_P3_cell(N_tri=N_tri, N3=N3, Y=Y, lift=True)
    ref = ref_elements(mesh)
    nfr = ref["nfr"]
    Se, Me, ae, tf = ref["Se"], ref["Me"], ref["ae"], ref["tf"]
    print(f"  ref: tets={len(ref['tets'])} CR faces/copy={nfr}  "
          f"top_faces={len(ref['top_pairs'])}")

    bfaces = exterior_faces(ref["mesh"], tf, nfr)
    print(f"  boundary faces: {len(bfaces)} "
          f"(floor={sum(1 for b in bfaces if b['kind']=='floor')}, "
          f"vert={sum(1 for b in bfaces if b['kind']=='vert')})")
    pair_maps = build_pair_maps(ref["mesh"], tf, nfr, bfaces=bfaces)

    # union-find on (copy, ref dof)
    parent = list(range(NC * nfr))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    n_glued = 0
    for name, perm in glue.items():
        pmap = pair_maps.get(name, {})
        if not pmap:
            continue
        for c in range(NC):
            j = perm[c]
            if j == c:
                continue  # Neumann self
            for d_src, d_dst in pmap.items():
                union(c * nfr + d_src, j * nfr + d_dst)
                n_glued += 1
    print(f"  cross-copy face glues applied: {n_glued}")

    gid = {}
    gmap = np.empty(NC * nfr, dtype=int)
    for x in range(NC * nfr):
        r = find(x)
        if r not in gid:
            gid[r] = len(gid)
        gmap[x] = gid[r]
    ng = len(gid)
    print(f"  global CR dofs: {ng}  (raw {NC * nfr}, merged {NC * nfr - ng})")

    rows, cols, qv, mv = [], [], [], []
    avec = np.zeros(ng)
    t_inf = np.zeros(ng)
    t_0 = np.zeros(ng)
    for c in range(NC):
        gd = gmap[c * nfr:(c + 1) * nfr]
        ix = gd[tf]
        rows.append(np.repeat(ix, 4, axis=1).ravel())
        cols.append(np.tile(ix, (1, 4)).ravel())
        qv.append(Se.ravel())
        mv.append(Me.ravel())
        np.add.at(avec, ix.ravel(), ae.ravel())
        tv = t_inf if cusp_class[c] == 0 else t_0
        for d, area in ref["top_pairs"]:
            tv[gd[d]] += area

    Q = coo_matrix(
        (np.concatenate(qv), (np.concatenate(rows), np.concatenate(cols))),
        (ng, ng),
    ).tocsr()
    M = coo_matrix(
        (np.concatenate(mv), (np.concatenate(rows), np.concatenate(cols))),
        (ng, ng),
    ).tocsr()

    one = np.ones(ng)
    print(f"  check 1'M1={float(one @ (M @ one)):.6f}  "
          f"|Q1|={np.abs(Q @ one).max():.2e}")
    print(f"  t_∞(1)={float(t_inf @ one):.6f} (want {ct.area_inf:.6f})  "
          f"t_0(1)={float(t_0 @ one):.6f} (want {ct.area_0:.6f})")
    ncomp, _ = connected_components((M + M.T) != 0, directed=False)
    print(f"  connected components (M graph): {ncomp}")

    return dict(
        Q=Q, M=M, a=avec, t_inf=t_inf, t_0=t_0, ng=ng, NC=NC, q=q,
        ct=ct, cusp_class=cusp_class, ncomp=ncomp, Y=Y,
    )


def constrained_mu_two_cusp(Q, M, a, t_inf, t_0, ct, Y, lam):
    """Min of A on L=0 for two-cusp pencil (dense if small)."""
    s = math.sqrt(1.0 - lam)
    b_inf = ct.beta_inf(s, Y)
    b_0 = ct.beta_0(s, Y)
    kap = 1.0 / ((1.0 + s) * Y * Y)
    Qd = Q.toarray() if hasattr(Q, "toarray") else Q
    Md = M.toarray() if hasattr(M, "toarray") else M
    A = (
        Qd - lam * Md
        - b_inf * np.outer(t_inf, t_inf)
        - b_0 * np.outer(t_0, t_0)
    )
    ell = a + kap * (t_inf + t_0)
    v = ell.copy()
    v[0] += np.sign(ell[0] if ell[0] != 0 else 1.0) * np.linalg.norm(ell)
    b = 2.0 / (v @ v)

    def hah(B):
        w = B @ v
        return (
            B - b * np.outer(v, w) - b * np.outer(w, v)
            + b ** 2 * (v @ w) * np.outer(v, v)
        )[1:, 1:]

    vals = eigh(hah(A), hah(Md), eigvals_only=True, subset_by_index=[0, 1])
    return float(vals[0]), float(vals[1])


def run(q=3, N_tri=4, N3=2, Y=Y_DEFAULT):
    print("Eisenstein Γ₀ float prototype (Stage 8)")
    print("=" * 60)
    print(f"  N(𝔭)={q}  mesh N_tri={N_tri} N3={N3}  Y={Y}")
    ct = cusp_table_prime(q)
    print(f"  cusps: |T_∞|={ct.area_inf:.6f} |T_0|={ct.area_0:.6f}  "
          f"modes OK={ct.ok_modes(Y)}")
    data = assemble_level(q=q, N_tri=N_tri, N3=N3, Y=Y)
    ng = data["ng"]
    if ng > 8000:
        print("  skip dense mu (dofs too large)")
        return data, []
    print("\n  two-cusp constrained μ:")
    print(f"  {'lambda':>8} {'s':>7} {'mu1':>10} {'mu2':>10}")
    lams = [0.05, 0.5, 0.9, 0.99]
    mus = []
    for lam in lams:
        s = math.sqrt(1 - lam)
        m1, m2 = constrained_mu_two_cusp(
            data["Q"], data["M"], data["a"],
            data["t_inf"], data["t_0"], data["ct"], Y, lam,
        )
        mus.append(m1)
        print(f"  {lam:8.3f} {s:7.4f} {m1:10.5f} {m2:10.5f}")
    mmin = min(mus)
    print(f"\n  min μ = {mmin:.5f}  "
          f"{'PASS (float)' if mmin > 0 else 'FAIL / nonpositive'}")
    if data["ncomp"] != 1:
        print(f"  NOTE: M-graph has {data['ncomp']} components "
              f"(gluing incomplete — geometric pair_maps partial)")
    return data, mus


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    q = int(argv[0]) if len(argv) >= 1 else 3
    # defaults: coarse for scan; N=3 needs ~6×3 for μ>0 near λ=1
    N_tri = int(argv[1]) if len(argv) >= 2 else (6 if q == 3 else 4)
    N3 = int(argv[2]) if len(argv) >= 3 else (3 if q == 3 else 2)
    data, mus = run(q=q, N_tri=N_tri, N3=N3)
    if mus and min(mus) <= 0:
        return 1
    if data.get("ncomp", 1) != 1:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
