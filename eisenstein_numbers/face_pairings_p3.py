"""Face-pairing dictionary for the EGM P_3 truncated core (multi-copy gluing).

The planar domain P_3 ([EGM98, §7.3], [DP20, §2.3]) has exterior edges

  LEFT  : x = 0,              0 ≤ y ≤ 1/√3
  RIGHT : x = 1/2,       |y| ≤ 1/(2√3)
  UP    : y = (1−x)/√3,       0 ≤ x ≤ 1/2
  LOW   : y = −x/√3,          0 ≤ x ≤ 1/2

plus the hemispherical floor |z|²+t²=1 and free top t=Y.

Poincaré face pairings used for Γ-gluing (names match residue generators
in residue_omega / Swan–EGM):

  T1 : RIGHT ↔ LEFT     parameter map (x,y,t) ↦ (x−1/2, y, t) on RIGHT
                        (equivalently opposite-side ID of the P_3 strip;
                        residue class of translation by 1)
  Tw : LOW   ↔ UP       opposite slant sides (hexagonal lattice direction)
  U  : cycle via z ↦ ωz when image of a vertical face centre lands on
       another vertical face (order-3 unit about 0)
  S  : FLOOR ↔ FLOOR    inversion in the unit sphere (self-paired)

These are the load-bearing maps for multi-copy CR assembly.  Discovery
by short words is available as a cross-check (`discover_pair_maps`) but
the **dictionary maps** (`build_pair_maps`) are the ones used in production.

References: GEOMETRY.md §4, T0_AREA.md, independent_exclusion/CONGRUENCE.md §3.
"""
from __future__ import annotations

import math

import numpy as np

from geometry_fund import GEN, OMEGA, mat_act

SQRT3 = math.sqrt(3.0)
S = 1.0 / SQRT3  # 1/√3


# ---------------------------------------------------------------------------
# Exact edge predicates (planar)
# ---------------------------------------------------------------------------

def _dist_point_seg(px, py, ax, ay, bx, by):
    abx, aby = bx - ax, by - ay
    t = ((px - ax) * abx + (py - ay) * aby) / max(abx * abx + aby * aby, 1e-30)
    t = max(0.0, min(1.0, t))
    qx, qy = ax + t * abx, ay + t * aby
    return math.hypot(px - qx, py - qy)


def classify_edge(x: float, y: float, eps: float = 0.04) -> str | None:
    """Which exterior edge of P_3 is (x,y) on? (by distance to segment)."""
    # segments
    edges = {
        "LEFT": (0.0, 0.0, 0.0, S),
        "RIGHT": (0.5, -0.5 * S, 0.5, 0.5 * S),
        "UP": (0.0, S, 0.5, 0.5 * S),
        "LOW": (0.0, 0.0, 0.5, -0.5 * S),
    }
    best, bd = None, 1e9
    for name, (ax, ay, bx, by) in edges.items():
        d = _dist_point_seg(x, y, ax, ay, bx, by)
        if d < bd:
            bd, best = d, name
    if bd <= eps:
        return best
    return None


# ---------------------------------------------------------------------------
# Boundary faces on mesh
# ---------------------------------------------------------------------------

def exterior_faces(mesh, tf, nfr, edge_eps: float = 0.045):
    """Boundary faces with kind and exact edge tag."""
    X, tets = mesh["X"], mesh["tets"]
    count = np.zeros(nfr, dtype=int)
    owner = np.full(nfr, -1, dtype=int)
    for e in range(len(tets)):
        for a in range(4):
            f = int(tf[e, a])
            count[f] += 1
            owner[f] = e
    bfaces = []
    for f in range(nfr):
        if count[f] != 1:
            continue
        e = owner[f]
        tet = tets[e]
        nodes = None
        for a in range(4):
            if int(tf[e, a]) == f:
                nodes = tuple(sorted(int(x) for x in np.delete(tet, a)))
                break
        P = X[list(nodes)]
        cen = P.mean(axis=0)
        x, y = float(cen[0]), float(cen[1])
        if mesh["is_floor"][list(nodes)].all():
            kind, edge = "floor", "FLOOR"
        elif mesh["is_top"][list(nodes)].all():
            kind, edge = "top", "TOP"
        else:
            kind = "vert"
            edge = classify_edge(x, y, eps=edge_eps) or "OTHER"
        bfaces.append(
            dict(fid=f, nodes=nodes, cen=cen, kind=kind, edge=edge, P=P)
        )
    return bfaces


def _match_by_height(src_list, dst_list, tol_t, also_y=True):
    """Greedy pair by hyperbolic height t (=x₃); optionally planar y."""
    m = {}
    used = set()
    # sort src by height for stable matching
    src_sorted = sorted(src_list, key=lambda b: b["cen"][2])
    for bf in src_sorted:
        best, bd = None, 1e9
        for tg in dst_list:
            if tg["fid"] in used:
                continue
            d = (float(bf["cen"][2]) - float(tg["cen"][2])) ** 2
            if also_y:
                d += 0.25 * (float(bf["cen"][1]) - float(tg["cen"][1])) ** 2
            if d < bd:
                bd, best = d, tg
        if best is not None and bd < tol_t * tol_t:
            m[bf["fid"]] = best["fid"]
            used.add(best["fid"])
    return m


# ---------------------------------------------------------------------------
# Dictionary pair maps (production)
# ---------------------------------------------------------------------------

def build_pair_maps(mesh, tf, nfr, bfaces=None, tol=None):
    """Return pair_maps for T1,Tw,U,S from the EGM edge dictionary.

    Also returns (bfaces, meta) for logging.
    """
    if bfaces is None:
        bfaces = exterior_faces(mesh, tf, nfr)
    if tol is None:
        hs = []
        for bf in bfaces[:40]:
            P = bf["P"]
            for i in range(3):
                hs.append(float(np.linalg.norm(P[i] - P[(i + 1) % 3])))
        hmed = float(np.median(hs)) if hs else 0.1
        tol = max(0.08, 0.9 * hmed)

    by_edge = {}
    for bf in bfaces:
        by_edge.setdefault(bf["edge"], []).append(bf)
    floors = by_edge.get("FLOOR", [])
    verts = [bf for bf in bfaces if bf["kind"] == "vert"]

    pair_maps = {"T1": {}, "Tw": {}, "U": {}, "S": {}}

    # T1: RIGHT ↔ LEFT (opposite sides of P_3)
    # Match by height only: edges have different planar y-ranges.
    rights = by_edge.get("RIGHT", [])
    lefts = by_edge.get("LEFT", [])
    m = _match_by_height(rights, lefts, tol_t=max(tol, 0.15), also_y=False)
    for s, t in m.items():
        pair_maps["T1"][s] = t
        pair_maps["T1"][t] = s

    # Tw: LOW ↔ UP (slanted opposite sides) — match by height
    lows = by_edge.get("LOW", [])
    ups = by_edge.get("UP", [])
    m = _match_by_height(lows, ups, tol_t=max(tol, 0.15), also_y=False)
    for s, t in m.items():
        pair_maps["Tw"][s] = t
        pair_maps["Tw"][t] = s

    # U: z ↦ ωz on vertical face centres → nearest vertical face
    for bf in verts:
        z = complex(float(bf["cen"][0]), float(bf["cen"][1])) * OMEGA
        img = np.array([z.real, z.imag, float(bf["cen"][2])])
        best, bd = None, 1e9
        for tg in verts:
            if tg["fid"] == bf["fid"]:
                continue
            d = float(np.sum((tg["cen"] - img) ** 2))
            if d < bd:
                bd, best = d, tg
        if best is not None and bd < (1.2 * tol) ** 2:
            pair_maps["U"][bf["fid"]] = best["fid"]

    # S: sphere inversion on floor
    for bf in floors:
        z = complex(float(bf["cen"][0]), float(bf["cen"][1]))
        t = float(bf["cen"][2])
        try:
            z2, t2 = mat_act(GEN["S"], z, t)
        except Exception:
            continue
        img = np.array([z2.real, z2.imag, t2])
        best, bd = None, 1e9
        for tg in floors:
            d = float(np.sum((tg["cen"] - img) ** 2))
            if d < bd:
                bd, best = d, tg
        if best is not None and bd < (1.5 * tol) ** 2:
            pair_maps["S"][bf["fid"]] = best["fid"]

    # soft involution fill
    for name in pair_maps:
        for s, t in list(pair_maps[name].items()):
            pair_maps[name].setdefault(t, s)

    stats = {k: len(v) for k, v in pair_maps.items()}
    edge_hist = {}
    for bf in bfaces:
        edge_hist[bf["edge"]] = edge_hist.get(bf["edge"], 0) + 1
    meta = dict(tol=tol, stats=stats, edge_hist=edge_hist, method="dictionary")
    return pair_maps, bfaces, meta


def validate_pair_maps(pair_maps, bfaces):
    by_fid = {bf["fid"]: bf for bf in bfaces}
    msgs = []
    for name, m in pair_maps.items():
        for s, t in m.items():
            if s not in by_fid or t not in by_fid:
                msgs.append(f"FAIL unknown fid {s}->{t}")
                continue
            ks, kt = by_fid[s]["kind"], by_fid[t]["kind"]
            if name == "S" and (ks != "floor" or kt != "floor"):
                msgs.append(f"FAIL S non-floor {s}->{t}")
            if name != "S" and (ks != "vert" or kt != "vert"):
                msgs.append(f"FAIL {name} non-vert {s}->{t}")
    ok = not any(x.startswith("FAIL") for x in msgs)
    return ok, msgs


def main():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from reference_cell import build_P3_cell
    from congruence_omega_proto import ref_elements

    print("EGM P_3 face-pairing dictionary")
    print("=" * 56)
    for N_tri, N3 in ((4, 2), (6, 3)):
        mesh = build_P3_cell(N_tri=N_tri, N3=N3, Y=1.25, lift=True)
        ref = ref_elements(mesh)
        pm, bf, meta = build_pair_maps(ref["mesh"], ref["tf"], ref["nfr"])
        ok, msgs = validate_pair_maps(pm, bf)
        print(f"\n  mesh {N_tri}x{N3}: tol={meta['tol']:.4f}")
        print(f"  edge_hist={meta['edge_hist']}")
        print(f"  stats={meta['stats']}  validate={'OK' if ok else 'FAIL'}")
        for m in msgs[:5]:
            print(f"    {m}")
        n_vert = sum(1 for b in bf if b["kind"] == "vert")
        n_floor = sum(1 for b in bf if b["kind"] == "floor")
        cov_v = len(set(pm["T1"]) | set(pm["Tw"]) | set(pm["U"]))
        print(f"  vert {n_vert} covered {cov_v}; floor {n_floor} "
              f"S-covered {len(pm['S'])}")


if __name__ == "__main__":
    main()
