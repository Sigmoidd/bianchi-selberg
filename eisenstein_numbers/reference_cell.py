"""Stage 2 — reference cell for Eisenstein–Picard core.

Primary mesh (default): EGM truncated core over P_3
  K_Y = { (z,y) : z ∈ P_3,  y_f(z) ≤ y ≤ Y },  y_f = √(1−|z|²)
  max|z| on P_3 = 1/√3 ⇒ y_f ≥ √(2/3) ≈ 0.816  (no Ford degeneracy)
  vol(K_Y) = vol(F) − |T|/(2 Y²)   [exact; geometry_fund.vol_KY_exact]
  |T| = area(P_3) = √3/6

Legacy mesh: R_comp parallelogram (u,v)∈[0,1]×[0,1/3] — hits |z|=1 at z=1,
so y_f→0; kept only for M0 regression via build_Rcomp_cell.

Side pairings for the true quotient (vs Neumann): GEOMETRY.md §4.
"""
from __future__ import annotations

import math
import numpy as np

SQRT3 = math.sqrt(3.0)
AREA_T = SQRT3 / 6.0  # area(P_3) = √3/6
Y_DEFAULT = 1.25
VMAX = 1.0 / 3.0  # legacy R_comp only
YF_MIN_P3 = math.sqrt(2.0 / 3.0)  # min y_f on P_3

# Kuhn split of unit cube -> 6 tets (legacy R_comp)
KUHN = [
    [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)],
    [(0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 1)],
    [(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)],
    [(0, 0, 0), (0, 1, 0), (0, 1, 1), (1, 1, 1)],
    [(0, 0, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1)],
    [(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)],
]


def exact_vol_KY(Y=Y_DEFAULT):
    """Exact EGM truncated-core volume vol(F)−|T|/(2Y²)."""
    try:
        from geometry_fund import vol_KY_exact
        return vol_KY_exact(Y)
    except Exception:
        vF = 0.16915693440160895
        tail = AREA_T / (2.0 * Y * Y)
        return vF - tail, 1e-12


def exact_vol_K_comp(Y=Y_DEFAULT, n=500):
    """Quadrature of hyp. volume of legacy R_comp prism (not EGM K_Y)."""
    try:
        from geometry_fund import vol_K_comp_quad
        return vol_K_comp_quad(Y, n=n)
    except Exception:
        return None


def z_from_uv(u, v):
    return u - 0.5 * v, (SQRT3 / 2.0) * v


def y_floor(x1, x2):
    r2 = x1 * x1 + x2 * x2
    return math.sqrt(max(1.0 - r2, 0.0))


def y_floor_arr(x1, x2):
    r2 = x1 ** 2 + x2 ** 2
    return np.sqrt(np.maximum(1.0 - r2, 0.0))


# ---------------------------------------------------------------------------
# P_3 vertices (EGM / DP20)
# ---------------------------------------------------------------------------
# T_up:  (0,0), (1/2, 1/(2√3)), (0, 1/√3)
# T_low: (0,0), (1/2, 1/(2√3)), (1/2, -1/(2√3))
_P3_A = (0.0, 0.0)
_P3_B = (0.5, 0.5 / SQRT3)
_P3_C = (0.0, 1.0 / SQRT3)
_P3_D = (0.5, -0.5 / SQRT3)


def _subdivide_triangle(V0, V1, V2, N):
    """Uniform barycentric subdivision into N² small triangles.

    Returns (points list of (x,y), tris as index triples).
    """
    pts = []
    idx = {}
    for i in range(N + 1):
        for j in range(N + 1 - i):
            a = 1.0 - (i + j) / N
            b = i / N
            c = j / N
            x = a * V0[0] + b * V1[0] + c * V2[0]
            y = a * V0[1] + b * V1[1] + c * V2[1]
            idx[(i, j)] = len(pts)
            pts.append((x, y))
    tris = []
    for i in range(N):
        for j in range(N - i):
            # up-pointing
            tris.append((idx[(i, j)], idx[(i + 1, j)], idx[(i, j + 1)]))
            # down-pointing (if room)
            if i + j + 1 < N:
                tris.append(
                    (idx[(i + 1, j)], idx[(i + 1, j + 1)], idx[(i, j + 1)])
                )
    return pts, tris


def _merge_base_meshes(pieces, tol=1e-12):
    """Merge (pts, tris) pieces with vertex welding by rounded coordinates."""
    key_to_id = {}
    pts_out = []
    tris_out = []

    def kid(x, y):
        return (round(x / tol) * tol, round(y / tol) * tol)

    for pts, tris in pieces:
        local = []
        for x, y in pts:
            k = kid(x, y)
            if k not in key_to_id:
                key_to_id[k] = len(pts_out)
                pts_out.append((x, y))
            local.append(key_to_id[k])
        for a, b, c in tris:
            tris_out.append((local[a], local[b], local[c]))
    return pts_out, tris_out


def build_P3_cell(N_tri=6, N3=4, Y=Y_DEFAULT, lift=True):
    """Tet mesh of EGM core K_Y over P_3.

    Parameters
    ----------
    N_tri : edge subdivisions per base triangle (N_tri² small tris each)
    N3    : vertical layers
    Y     : truncation height
    lift  : Lemma-G style floor lift (sag of sphere over each base edge)

    Returns mesh dict compatible with cr_omega.assemble_cr.
    """
    pieces = [
        _subdivide_triangle(_P3_A, _P3_B, _P3_C, N_tri),
        _subdivide_triangle(_P3_A, _P3_B, _P3_D, N_tri),
    ]
    base_pts, base_tris = _merge_base_meshes(pieces)
    nb = len(base_pts)
    xy = np.array(base_pts)

    # sphere floor + optional lift
    base_yf = np.array([y_floor(x, y) for x, y in base_pts])
    assert base_yf.min() >= YF_MIN_P3 - 1e-9, base_yf.min()

    lift_b = np.zeros(nb)
    if lift:
        # per base triangle: diam^2 / (8 yf_min^3)
        cell_lift = np.zeros(len(base_tris))
        for t, (a, b, c) in enumerate(base_tris):
            corners = [xy[a], xy[b], xy[c]]
            d2 = 0.0
            for p in range(3):
                for q in range(p + 1, 3):
                    d2 = max(d2, float(np.sum((corners[p] - corners[q]) ** 2)))
            yf_min = min(base_yf[a], base_yf[b], base_yf[c])
            cell_lift[t] = 0.125 * d2 / (yf_min ** 3) * (1 + 1e-9)
        # node lift = max over incident tris
        for t, (a, b, c) in enumerate(base_tris):
            for n in (a, b, c):
                lift_b[n] = max(lift_b[n], cell_lift[t])

    floor_y = base_yf + lift_b
    delta_bar = float(lift_b.max())

    # 3D nodes: base_id * (N3+1) + k
    nn = nb * (N3 + 1)
    X = np.empty((nn, 3))
    is_top = np.zeros(nn, dtype=bool)
    is_floor = np.zeros(nn, dtype=bool)
    lift_node = np.zeros(nn)

    def nid(b, k):
        return b * (N3 + 1) + k

    for b in range(nb):
        for k in range(N3 + 1):
            tt = k / N3
            n = nid(b, k)
            y = floor_y[b] * (1 - tt) + Y * tt
            X[n] = (xy[b, 0], xy[b, 1], y)
            is_floor[n] = k == 0
            is_top[n] = k == N3
            if k == 0:
                lift_node[n] = lift_b[b]

    # prism → 3 tets per base tri × layer.
    # Conformity: sort base vertices by base-index so every extruded quad
    # gets the same diagonal on both sides (standard ordered-prism split).
    tets = []
    for tri in base_tris:
        order = tuple(sorted(tri))  # base ids ascending
        # map original top/bottom after sorting
        # order = (p,q,r) with p<q<r as base indices
        p, q, r = order
        for k in range(N3):
            p0, q0, r0 = nid(p, k), nid(q, k), nid(r, k)
            p1, q1, r1 = nid(p, k + 1), nid(q, k + 1), nid(r, k + 1)
            tets.append([p0, q0, r0, r1])
            tets.append([p0, q0, q1, r1])
            tets.append([p0, p1, q1, r1])
    tets = np.array(tets, dtype=int)

    # top faces = base tris at k=N3
    top_faces = [
        (nid(a, N3), nid(b, N3), nid(c, N3)) for a, b, c in base_tris
    ]

    # mean planar edge length (diagnostic)
    h_xy = 0.0
    cnt = 0
    for a, b, c in base_tris:
        for p, q in ((a, b), (b, c), (c, a)):
            h_xy += float(np.linalg.norm(xy[p] - xy[q]))
            cnt += 1
    h_xy = h_xy / max(cnt, 1)

    return dict(
        X=X,
        tets=tets,
        is_top=is_top,
        is_floor=is_floor,
        lift_node=lift_node,
        top_faces=top_faces,
        Y=Y,
        N=(N_tri, N_tri, N3),
        N_tri=N_tri,
        N3=N3,
        AREA_T=AREA_T,
        VMAX=VMAX,
        delta_bar=delta_bar,
        y_floor_min=float(base_yf.min()),
        domain="P3",
        n_base=nb,
        n_base_tris=len(base_tris),
        h_xy=h_xy,
    )


def build_Rcomp_cell(N1, N2, N3, Y=Y_DEFAULT, lift=True):
    """Legacy parallelogram mesh (hits |z|=1). Prefer build_P3_cell."""
    ug = np.linspace(0.0, 1.0, N1 + 1)
    vg = np.linspace(0.0, VMAX, N2 + 1)
    hu, hv = 1.0 / N1, VMAX / N2

    base = np.empty((N1 + 1, N2 + 1))
    for i in range(N1 + 1):
        for j in range(N2 + 1):
            x1, x2 = z_from_uv(ug[i], vg[j])
            base[i, j] = y_floor(x1, x2)

    lift_n = np.zeros((N1 + 1, N2 + 1))
    if lift:
        cell_lift = np.empty((N1, N2))
        for i in range(N1):
            for j in range(N2):
                corners = [
                    z_from_uv(ug[i + a], vg[j + b]) for a in (0, 1) for b in (0, 1)
                ]
                d2 = 0.0
                for p in range(4):
                    for q in range(p + 1, 4):
                        dx = corners[p][0] - corners[q][0]
                        dy = corners[p][1] - corners[q][1]
                        d2 = max(d2, dx * dx + dy * dy)
                yf_min = min(
                    base[i, j], base[i + 1, j], base[i, j + 1], base[i + 1, j + 1]
                )
                yf_min = max(yf_min, 0.15)
                raw = 0.125 * d2 / (yf_min ** 3) * (1 + 1e-9)
                cell_lift[i, j] = min(raw, 0.08)
        for i in range(N1 + 1):
            for j in range(N2 + 1):
                adj = cell_lift[
                    max(i - 1, 0) : min(i + 1, N1),
                    max(j - 1, 0) : min(j + 1, N2),
                ]
                lift_n[i, j] = adj.max() if adj.size else 0.0

    floor_y = base + lift_n
    delta_bar = float(lift_n.max())

    nn = (N1 + 1) * (N2 + 1) * (N3 + 1)
    X = np.empty((nn, 3))
    is_top = np.zeros(nn, dtype=bool)
    is_floor = np.zeros(nn, dtype=bool)
    lift_node = np.zeros(nn)

    def nid(i, j, k):
        return (i * (N2 + 1) + j) * (N3 + 1) + k

    for i in range(N1 + 1):
        for j in range(N2 + 1):
            x1, x2 = z_from_uv(ug[i], vg[j])
            for k in range(N3 + 1):
                tt = k / N3
                n = nid(i, j, k)
                y = floor_y[i, j] * (1 - tt) + Y * tt
                X[n] = (x1, x2, y)
                is_floor[n] = k == 0
                is_top[n] = k == N3
                if k == 0:
                    lift_node[n] = lift_n[i, j]

    tets = []
    for i in range(N1):
        for j in range(N2):
            for k in range(N3):
                for tt in KUHN:
                    tets.append(
                        [nid(i + a, j + b, k + c) for (a, b, c) in tt]
                    )
    tets = np.array(tets, dtype=int)

    top_faces = []
    for i in range(N1):
        for j in range(N2):
            c00 = nid(i, j, N3)
            c10 = nid(i + 1, j, N3)
            c01 = nid(i, j + 1, N3)
            c11 = nid(i + 1, j + 1, N3)
            top_faces.append((c00, c10, c11))
            top_faces.append((c00, c11, c01))

    return dict(
        X=X,
        tets=tets,
        is_top=is_top,
        is_floor=is_floor,
        lift_node=lift_node,
        top_faces=top_faces,
        Y=Y,
        N=(N1, N2, N3),
        AREA_T=AREA_T,
        VMAX=VMAX,
        delta_bar=delta_bar,
        y_floor_min=float(base.min()),
        domain="Rcomp",
        hu=hu,
        hv=hv,
    )


def build_reference_cell(N1=6, N2=6, N3=4, Y=Y_DEFAULT, lift=True, domain="P3"):
    """Build the computational reference cell.

    domain='P3'   (default): EGM planar section — correct volume & y_f ≥ √(2/3)
    domain='Rcomp': legacy parallelogram (N1×N2×N3); degenerates at |z|=1

    For domain='P3', N1 is used as N_tri (edge subdivisions); N2 ignored;
    N3 = vertical layers.
    """
    if domain == "P3":
        return build_P3_cell(N_tri=N1, N3=N3, Y=Y, lift=lift)
    if domain == "Rcomp":
        return build_Rcomp_cell(N1, N2, N3, Y=Y, lift=lift)
    raise ValueError(f"unknown domain {domain!r}")
