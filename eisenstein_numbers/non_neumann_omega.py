"""Optional: level-1 CR space with Poincaré face pairings (non-Neumann).

Imposes the P₃ edge dictionary (T1, Tw, U, S) by identifying CR face-mean
dofs on paired exterior faces.  TOP faces stay free (cusp ODE / β).

Why this is optional for the theorem
------------------------------------
The Neumann free H¹(K) is a *larger* test space.  Theorem G1 positivity
there already implies positivity on the quotient subspace of Γ-invariant
(face-paired) functions.  This module:

  1. Builds the genuine orbifold FE space for the paper's geometric story.
  2. Checks that the reduced system still has ker(Q)≃constants, t(1)=|T|,
     μ>0, and float G1 windows pass.
  3. Compares Lemma I1 vs sharper κ=CZZ (κ₁) on the same space.

Usage:
  python -u non_neumann_omega.py           # 6×3, both κ modes
  python -u non_neumann_omega.py 6 3
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).resolve().parent))

from constants import KAPPA_1, KAPPA_SC, g1_coeffs  # noqa: E402
from face_pairings_p3 import (  # noqa: E402
    build_pair_maps, exterior_faces, validate_pair_maps,
)
from reference_cell import AREA_T, Y_DEFAULT, build_P3_cell  # noqa: E402
from cr_omega import (  # noqa: E402
    assemble_cr, compute_g1_constants, constrained_mu_dense, window_check,
)


def _union_find(n: int):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
            return True
        return False

    return find, union, parent


def assemble_cr_paired(mesh, pair_maps=None, bfaces=None):
    """CR assembly with exterior face IDs from the pairing dictionary.

    Returns reduced Q,M,a,t and maps raw face → global dof.
    """
    raw = assemble_cr(mesh)
    nfr = raw["nfr"]
    tf = raw["tet_faces"]
    Se = None  # re-use raw sparse by remapping

    if bfaces is None:
        # need face tables: rebuild via exterior_faces needs tf/nfr
        from congruence_omega_proto import ref_elements
        ref = ref_elements(mesh)
        # align: use same mesh orientation as assemble_cr
        bfaces = exterior_faces(raw["mesh"], tf, nfr)
        if pair_maps is None:
            pair_maps, bfaces, meta = build_pair_maps(
                raw["mesh"], tf, nfr, bfaces=bfaces
            )
        else:
            meta = {}
    else:
        if pair_maps is None:
            pair_maps, bfaces, meta = build_pair_maps(
                raw["mesh"], tf, nfr, bfaces=bfaces
            )
        else:
            meta = {}

    vok, msgs = validate_pair_maps(pair_maps, bfaces)
    find, union, _ = _union_find(nfr)
    n_id = 0
    for name, pmap in pair_maps.items():
        for s, t in pmap.items():
            if s < t and union(s, t):
                n_id += 1

    # tops must not be glued to anything non-top
    top_fids = {b["fid"] for b in bfaces if b["kind"] == "top"}
    for b in bfaces:
        if b["kind"] == "top":
            # ensure top stays free: no union already applied to tops
            pass
    glued_top = any(
        find(f) != f or any(find(f) == find(g) and f != g for g in range(nfr))
        for f in top_fids
    )
    # simpler: check no top appears in pair_maps
    hit_top = any(set(m) & top_fids for m in pair_maps.values())

    gid = {}
    gmap = np.empty(nfr, dtype=int)
    for x in range(nfr):
        r = find(x)
        if r not in gid:
            gid[r] = len(gid)
        gmap[x] = gid[r]
    ng = len(gid)

    # reassemble element matrices into reduced dofs
    X, tets = raw["mesh"]["X"], raw["mesh"]["tets"]
    # reuse per-tet Se, Me, ae from raw: reconstruct by local re-assembly
    # Easier path: scatter raw dense contributions via gmap
    Qr = raw["Q"].tocoo()
    Mr = raw["M"].tocoo()
    rows_q = gmap[Qr.row]
    cols_q = gmap[Qr.col]
    rows_m = gmap[Mr.row]
    cols_m = gmap[Mr.col]
    Q = coo_matrix((Qr.data, (rows_q, cols_q)), shape=(ng, ng)).tocsr()
    M = coo_matrix((Mr.data, (rows_m, cols_m)), shape=(ng, ng)).tocsr()
    a = np.zeros(ng)
    t = np.zeros(ng)
    np.add.at(a, gmap, raw["a"])
    np.add.at(t, gmap, raw["t"])

    return dict(
        Q=Q, M=M, a=a, t=t, nfr=ng, nfr_raw=nfr,
        AREA_T=raw["AREA_T"], Y=raw["Y"],
        gmap=gmap, pair_maps=pair_maps, bfaces=bfaces,
        meta=meta, validate_ok=vok, validate_msgs=msgs,
        n_identifications=n_id, hit_top=hit_top,
        # keep geometry for G1 constants from *unpaired* mesh (same)
        wQ=raw["wQ"], wM=raw["wM"], vol=raw["vol"], hT=raw["hT"],
        ymin=raw["ymin"], ymax=raw["ymax"], grads=raw["grads"],
        tet_faces=tf, face_area=raw["face_area"],
        top_tets=raw["top_tets"], floor_tets=raw["floor_tets"],
        mesh=raw["mesh"],
    )


def run_non_neumann(
    N_tri: int = 6,
    N3: int = 3,
    Y: float = Y_DEFAULT,
    nwin: int = 8,
    rho: float = 55.0,
    nu_star: float = 1.05,
    do_mu: bool = True,
) -> dict:
    print("Level-1 non-Neumann CR (face pairings imposed)")
    print("=" * 64)
    print(f"  mesh P3 {N_tri}x{N3}, Y={Y}, rho={rho}")

    mesh = build_P3_cell(N_tri=N_tri, N3=N3, Y=Y, lift=True)
    data = assemble_cr_paired(mesh)
    stats = {k: len(v) for k, v in data["pair_maps"].items()}
    print(f"  dictionary validate={'OK' if data['validate_ok'] else 'FAIL'} "
          f"stats={stats}")
    print(f"  IDs applied={data['n_identifications']}")
    print(f"  CR dofs: raw={data['nfr_raw']} → paired={data['nfr']} "
          f"(merged {data['nfr_raw'] - data['nfr']})")
    print(f"  tops free (no TOP in maps): {not data['hit_top']}")

    one = np.ones(data["nfr"])
    m1 = float(one @ (data["M"] @ one))
    t1 = float(data["t"] @ one)
    q1 = float(np.abs(data["Q"] @ one).max())
    print(f"  1'M1={m1:.6f}  t(1)={t1:.6f} (want {AREA_T:.6f})  |Q1|={q1:.2e}")
    ncomp, _ = connected_components(
        (data["M"] + data["M"].T) != 0, directed=False
    )
    print(f"  M-graph components: {ncomp}")

    ok_geom = (
        data["validate_ok"]
        and not data["hit_top"]
        and abs(t1 - AREA_T) < 1e-6
        and q1 < 1e-8
        and ncomp == 1
    )

    results = {}
    for kappa_name, kappa in (("I1", KAPPA_SC), ("CZZ", KAPPA_1)):
        print(f"\n  --- kappa={kappa_name} ({kappa:.6f}) ---")
        cst = compute_g1_constants(data, kappa=kappa)
        print(
            f"  gamma={cst['gamma']:.4f} tau={cst['tau']:.4f} "
            f"alpha_h={cst['alpha_h']:.4f} S_Q={cst['S_Q']:.3e}"
        )
        Qd = data["Q"].toarray()
        Md = data["M"].toarray()
        a, t = data["a"], data["t"]
        sgrid = np.linspace(0.0, 1.0, nwin + 1)
        n_ok = 0
        min_ratio = float("inf")
        mus = []
        print(f"  {'window':>14} {'c_e/d_e':>9} {'PSD':>5} {'mu':>9}")
        for w in range(nwin):
            co = g1_coeffs(
                sgrid[w], sgrid[w + 1], data["AREA_T"], Y, cst,
                rho=rho, nu_star=nu_star,
            )
            r = window_check(Qd, Md, a, t, co)
            mu = float("nan")
            if do_mu and data["nfr"] <= 6000:
                lam = min(1.0 - sgrid[w] ** 2, 0.999)
                try:
                    mu = constrained_mu_dense(
                        Qd, Md, a, t, data["AREA_T"], Y, lam
                    )
                except Exception:
                    mu = float("nan")
            if r["ok"]:
                n_ok += 1
            min_ratio = min(min_ratio, r["ratio"])
            mus.append(mu)
            print(
                f"  [{sgrid[w]:.2f},{sgrid[w+1]:.2f}] "
                f"{r['ratio']:9.3f} "
                f"{'yes' if r['psd'] else ' no':>5} {mu:9.4f}  "
                f"{'PASS' if r['ok'] else 'fail'}"
            )
        n_mu = sum(1 for m in mus if m == m and m > 0)
        print(f"  G1 PASS: {n_ok}/{nwin}  mu>0: {n_mu}/{nwin}  "
              f"min c_e/d_e={min_ratio:.3f}")
        results[kappa_name] = dict(
            n_ok=n_ok, n_mu=n_mu, min_ratio=min_ratio, cst=cst, mus=mus,
        )

    # Summary
    print("\n" + "=" * 64)
    i1 = results["I1"]
    czz = results["CZZ"]
    ok_float = i1["n_ok"] == nwin and i1["n_mu"] == nwin and ok_geom
    print(f"OVERALL non-Neumann float: {'PASS' if ok_float else 'FAIL'}")
    print(f"  geometry/IDs/t(1)/ker(Q): {'OK' if ok_geom else 'FAIL'}")
    print(f"  G1 I1:  {i1['n_ok']}/{nwin}  min ratio={i1['min_ratio']:.3f}")
    print(f"  G1 CZZ: {czz['n_ok']}/{nwin}  min ratio={czz['min_ratio']:.3f}")
    print("  Note: Neumann cert already implies this subspace; this is the")
    print("  geometric (paired) FE space for the paper writeup.")
    return dict(
        ok=ok_float, data=data, results=results, ok_geom=ok_geom,
    )


def main() -> int:
    N_tri, N3 = 6, 3
    if len(sys.argv) >= 3:
        N_tri, N3 = int(sys.argv[1]), int(sys.argv[2])
    out = run_non_neumann(N_tri=N_tri, N3=N3)
    return 0 if out["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
