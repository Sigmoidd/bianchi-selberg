#!/usr/bin/env python3
"""
V_h^{P1,per} ⊂ H¹(Γ\\H³) and companion J_h CR face-mean jump bridge.

Purpose
-------
Make concrete the morningbrief warning:

  P1 Neumann with *dropped* pairing is lower-bound flavor only.
  Dual / defect upper bounds need a conforming periodic trial space
      V_h^{P1,per}  ⊂  H¹(Γ \\ H³)
  or a companion jump functional J_h on free CR that measures the defect
  of not being periodic.

This module:
  1. Builds the multi-copy CR face-dof gluing for Γ₀(2+i) (N=5) using
     independent_exclusion.congruence_prototype PAIRINGS + glue.
  2. Defines V_h^{P1,per} as the image of the CR space under face-mean
     identification on *cross-copy* pairings (self-gluings stay free =
     Neumann relaxation of interior self-faces, as in CONGRUENCE.md).
  3. Defines J_h(u) = max over pairing faces of |mean(u, F) − mean(u, φ(F))|
     on the *unreduced* multi-copy CR vector (before identification).
  4. Diagnostics: constants have J_h≈0; free-Neumann random vectors have
     J_h = O(1); reduced (periodic) vectors have J_h≡0 by construction.

Language
--------
Engineering scaffold — not a certified dual upper bound, not Rung-4 green.
Does not change hard map. Does not claim Rayleigh on free Neumann is an
upper bound for the quotient spectrum.

Usage:
  python p1_per_jh_bridge.py
  python p1_per_jh_bridge.py --N1 4 --N2 2 --N3 2
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
_IE = _ROOT / "independent_exclusion"
sys.path.insert(0, str(_IE))

from congruence_prototype import (  # noqa: E402
    INDEX,
    LEVEL,
    NP,
    PAIRINGS,
    build_gluing,
    build_reference,
    ref_geometry,
    set_level,
    tri_key,
)


def build_face_tables(
    N1: int = 4,
    N2: int = 2,
    N3: int = 2,
) -> Dict[str, Any]:
    """Reference-cell CR face table + boundary pairing maps (one copy)."""
    X, tets, btri = build_reference(N1, N2, N3)
    nt = len(tets)
    fid: Dict[Tuple[int, ...], int] = {}
    tf = np.empty((nt, 4), dtype=int)
    face_nodes: List[Tuple[int, int, int]] = []
    for e in range(nt):
        for a in range(4):
            key = tuple(sorted(np.delete(tets[e], a).tolist()))
            if key not in fid:
                fid[key] = len(fid)
                face_nodes.append(key)  # type: ignore[arg-type]
            tf[e, a] = fid[key]
    nfr = len(fid)

    def dof_of(tri: Sequence[int]) -> int:
        return fid[tuple(sorted(tri))]

    pair_maps: Dict[str, Dict[int, int]] = {}
    for name, src, dst, mp in PAIRINGS:
        dst_lookup = {tri_key(X, t): dof_of(t) for t in btri[dst]}
        m: Dict[int, int] = {}
        for t in btri[src]:
            key = tri_key(X, t, mapping=mp)
            m[dof_of(t)] = dst_lookup[key]
        pair_maps[name] = m

    # face areas (Euclidean) for weighted means
    face_area = np.zeros(nfr)
    for e in range(nt):
        for a in range(4):
            nodes = list(np.delete(tets[e], a))
            P = X[nodes]
            area = 0.5 * np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))
            face_area[tf[e, a]] = area  # each interior face written twice ok

    return dict(
        X=X,
        tets=tets,
        btri=btri,
        tf=tf,
        nfr=nfr,
        nt=nt,
        face_nodes=face_nodes,
        pair_maps=pair_maps,
        face_area=face_area,
        N=(N1, N2, N3),
    )


def build_periodic_identification(
    nfr: int,
    pair_maps: Dict[str, Dict[int, int]],
    glue: Dict[str, List[int]],
    nc: int = INDEX,
) -> Dict[str, Any]:
    """
    Union-find on (copy, ref_face) for cross-copy pairings only.

    Self-identifications (glue[c]==c) are *not* merged — Neumann relaxation
    of interior self-faces, matching CONGRUENCE.md / assemble_level_p.
    """
    parent = list(range(nc * nfr))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    n_merged_edges = 0
    n_self_relaxed = 0
    for name, pmap in pair_maps.items():
        perm = glue[name]
        for c in range(nc):
            j = perm[c]
            if j == c:
                n_self_relaxed += len(pmap)
                continue
            for d_src, d_dst in pmap.items():
                union(c * nfr + d_src, j * nfr + d_dst)
                n_merged_edges += 1

    gid: Dict[int, int] = {}
    for x in range(nc * nfr):
        r = find(x)
        if r not in gid:
            gid[r] = len(gid)
    gmap = np.array([gid[find(x)] for x in range(nc * nfr)], dtype=int)
    n_raw = nc * nfr
    n_per = len(gid)
    return dict(
        gmap=gmap,
        n_raw=n_raw,
        n_per=n_per,
        n_merged_edges=n_merged_edges,
        n_self_relaxed=n_self_relaxed,
        reduction_ratio=n_raw / max(n_per, 1),
        language=(
            "V_h^{P1,per} ≅ CR face-means modulo cross-copy PAIRINGS gluing; "
            "self-faces free (Neumann). Subspace of multi-copy CR, not of free "
            "single-cell Neumann for dual upper bounds."
        ),
    )


def lift_periodic_to_raw(v_per: np.ndarray, gmap: np.ndarray) -> np.ndarray:
    """Embed periodic dofs into unreduced multi-copy face vector."""
    return v_per[gmap]


def project_raw_to_periodic(v_raw: np.ndarray, gmap: np.ndarray, n_per: int) -> np.ndarray:
    """Average raw face values onto periodic dofs (L2 face-mean projection)."""
    out = np.zeros(n_per, dtype=float)
    cnt = np.zeros(n_per, dtype=float)
    for i, g in enumerate(gmap):
        out[g] += v_raw[i]
        cnt[g] += 1.0
    cnt = np.maximum(cnt, 1.0)
    return out / cnt


def J_h(
    v_raw: np.ndarray,
    nfr: int,
    pair_maps: Dict[str, Dict[int, int]],
    glue: Dict[str, List[int]],
    face_area: np.ndarray,
    nc: int = INDEX,
    only_cross_copy: bool = False,
) -> Dict[str, Any]:
    """
    Companion jump functional on unreduced multi-copy CR face means:

      J_h(u) = max | u_{c,F} − u_{j, φ_γ(F)} |,  j = glue[γ](c).

    only_cross_copy=True: only j≠c terms — these *must* vanish on V_h^{P1,per}
    (the identification subspace). Self-faces (j=c) are Neumann-relaxed and
    need not vanish; report them separately as J_self when only_cross_copy=False.
    """
    per_gen: Dict[str, float] = {}
    max_all = 0.0
    max_cross = 0.0
    max_self = 0.0
    rms_sq = 0.0
    n_terms = 0
    n_cross = 0
    n_self = 0
    for name, pmap in pair_maps.items():
        perm = glue[name]
        jumps: List[float] = []
        for c in range(nc):
            j = perm[c]
            cross = j != c
            if only_cross_copy and not cross:
                continue
            for d_src, d_dst in pmap.items():
                i0 = c * nfr + d_src
                i1 = j * nfr + d_dst
                w = math.sqrt(max(float(face_area[d_src]), 1e-30))
                jmp = abs(float(v_raw[i0]) - float(v_raw[i1]))
                jumps.append(jmp)
                rms_sq += (jmp * w) ** 2
                n_terms += 1
                if cross:
                    max_cross = max(max_cross, jmp)
                    n_cross += 1
                else:
                    max_self = max(max_self, jmp)
                    n_self += 1
        dmax = max(jumps) if jumps else 0.0
        per_gen[name] = dmax
        max_all = max(max_all, dmax)
    return dict(
        J_h_max=max_all,
        J_h_cross=max_cross,
        J_h_self=max_self,
        J_h_rms=math.sqrt(rms_sq / max(n_terms, 1)),
        per_generator=per_gen,
        n_terms=n_terms,
        n_cross=n_cross,
        n_self=n_self,
    )


def diagnose(
    N1: int = 4,
    N2: int = 2,
    N3: int = 2,
    seed: int = 0,
) -> Dict[str, Any]:
    t0 = time.time()
    set_level("(2+i)")
    ref = build_face_tables(N1, N2, N3)
    pts, glue, cusp_class = build_gluing("(2+i)")
    assert len(pts) == INDEX
    ident = build_periodic_identification(ref["nfr"], ref["pair_maps"], glue, INDEX)

    n_raw = ident["n_raw"]
    n_per = ident["n_per"]
    gmap = ident["gmap"]
    rng = np.random.default_rng(seed)

    fa = ref["face_area"]
    pm = ref["pair_maps"]
    nfr = ref["nfr"]

    # 1) constant field → all jumps 0
    v_const = np.ones(n_raw)
    j_const = J_h(v_const, nfr, pm, glue, fa)

    # 2) free Neumann random (unreduced) → cross J typically O(1)
    v_free = rng.standard_normal(n_raw)
    v_free = v_free - v_free.mean()
    nrm = np.linalg.norm(v_free)
    if nrm > 0:
        v_free = v_free / nrm * math.sqrt(n_raw)
    j_free = J_h(v_free, nfr, pm, glue, fa)
    j_free_cross = J_h(v_free, nfr, pm, glue, fa, only_cross_copy=True)

    # 3) periodic: random on V_h^{P1,per}, lift → cross J ≡ 0
    v_per = rng.standard_normal(n_per)
    v_per = v_per - v_per.mean()
    nrm = np.linalg.norm(v_per)
    if nrm > 0:
        v_per = v_per / nrm * math.sqrt(n_per)
    v_lift = lift_periodic_to_raw(v_per, gmap)
    j_per_all = J_h(v_lift, nfr, pm, glue, fa)
    j_per_cross = J_h(v_lift, nfr, pm, glue, fa, only_cross_copy=True)

    # 4) project free → periodic → lift: cross J ≡ 0
    v_proj = project_raw_to_periodic(v_free, gmap, n_per)
    v_proj_lift = lift_periodic_to_raw(v_proj, gmap)
    j_proj_cross = J_h(v_proj_lift, nfr, pm, glue, fa, only_cross_copy=True)
    defect_free_to_per = float(
        np.linalg.norm(v_free - v_proj_lift) / max(np.linalg.norm(v_free), 1e-300)
    )

    out = dict(
        level=LEVEL,
        NP=NP,
        INDEX=INDEX,
        mesh=dict(N1=N1, N2=N2, N3=N3, nfr=nfr, nt=ref["nt"]),
        identification=dict(
            n_raw=n_raw,
            n_per=n_per,
            n_merged_edges=ident["n_merged_edges"],
            n_self_relaxed=ident["n_self_relaxed"],
            reduction_ratio=ident["reduction_ratio"],
            language=ident["language"],
        ),
        J_h_constant=j_const,
        J_h_free_neumann=j_free,
        J_h_free_cross_only=j_free_cross,
        J_h_periodic_lift_all=j_per_all,
        J_h_periodic_cross_only=j_per_cross,
        J_h_projected_cross_only=j_proj_cross,
        free_to_periodic_L2_defect=defect_free_to_per,
        warnings=[
            "Free Neumann Rayleigh is NOT an upper bound for Spec(Γ\\H³).",
            "V_h^{P1,per} kills only cross-copy jumps; self-faces stay Neumann.",
            "J_h_cross is the companion defect of a free trial vs periodic space.",
            "This scaffold does not assemble Q/M on V_h^{P1,per} for dual cert.",
        ],
        checks=dict(
            constant_J_near_zero=j_const["J_h_max"] < 1e-12,
            periodic_cross_J_near_zero=j_per_cross["J_h_max"] < 1e-12,
            projected_cross_J_near_zero=j_proj_cross["J_h_max"] < 1e-12,
            free_cross_J_positive=j_free_cross["J_h_max"] > 1e-6,
        ),
        seconds=time.time() - t0,
    )
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--N1", type=int, default=4)
    p.add_argument("--N2", type=int, default=2)
    p.add_argument("--N3", type=int, default=2)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--json-out", type=str, default="")
    args = p.parse_args(argv)

    print("=== V_h^{P1,per} / J_h bridge (engineering scaffold) ===")
    out = diagnose(N1=args.N1, N2=args.N2, N3=args.N3, seed=args.seed)
    idn = out["identification"]
    print(f"level={out['level']}  NP={out['NP']}  copies={out['INDEX']}")
    print(
        f"mesh N={out['mesh']['N1']}x{out['mesh']['N2']}x{out['mesh']['N3']}  "
        f"nfr/copy={out['mesh']['nfr']}  tets/copy={out['mesh']['nt']}"
    )
    print(
        f"dofs: raw={idn['n_raw']}  periodic={idn['n_per']}  "
        f"reduction={idn['reduction_ratio']:.3f}×  "
        f"merged_edges={idn['n_merged_edges']}  self_relaxed={idn['n_self_relaxed']}"
    )
    print("\nJ_h diagnostics (face-mean jumps):")
    for label, key in [
        ("constant", "J_h_constant"),
        ("free Neumann (all)", "J_h_free_neumann"),
        ("free cross-only", "J_h_free_cross_only"),
        ("periodic all", "J_h_periodic_lift_all"),
        ("periodic cross-only", "J_h_periodic_cross_only"),
        ("proj free→per cross", "J_h_projected_cross_only"),
    ]:
        j = out[key]
        print(
            f"  {label:22s}  J_max={j['J_h_max']:.3e}  "
            f"J_cross={j.get('J_h_cross', float('nan')):.3e}  "
            f"J_self={j.get('J_h_self', float('nan')):.3e}"
        )
    print(f"\n  free→periodic L² defect = {out['free_to_periodic_L2_defect']:.4f}")
    print("  checks:", out["checks"])
    for w in out["warnings"]:
        print(f"  WARN: {w}")
    print(f"  ({out['seconds']:.2f}s)")

    ok = all(out["checks"].values())
    path = Path(args.json_out) if args.json_out else _HERE / "p1_per_jh_result.json"
    path.write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    print(f"  wrote {path}")
    print("  PASS" if ok else "  FAIL checks")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
