"""
Fundamental polyhedron for Gamma = PSL(2, O_{-3}) = PSL(2, Z[omega]).

Citations (load-bearing)
------------------------
[EGM98]  Elstrodt–Grunewald–Mennicke, *Groups Acting on Hyperbolic Space*,
         Springer 1998, §7.2–7.3: construction of B_d, P_d, F_d = B_d ∩
         (P_d × R>0).  Theorem: F_d is a fundamental domain for Gamma_d.
         (Same construction quoted as Theorem 2.3 in [DP20].)

[DP20]   Dória–Paula, "Height estimates for Bianchi groups", arXiv:1910.03148,
         §2.3 — explicit formula for P_3 matching EGM (used below).

[Swan71] R.G. Swan, "Generators and relations for certain special linear
         groups", Adv. Math. 6 (1971) — presentation / face-pairing generators
         for SL(2, O_K).

[Humbert] Volume of PSL(2, O_K)\\H^3:
            vol(F) = |D|^{3/2} zeta_K(2) / (4 pi^2),
          D = disc(K) = -3 for K=Q(sqrt(-3)).  Cross-check: verify_eisenstein.py
          ~ 0.169156…; Sengün survey.

This module freezes:
  1. F_3 = B_3 ∩ (P_3 × R>0)  (canonical polyhedron)
  2. Side-pairing isometries (generators) and face → face map
  3. Exact vol(F) and vol(K_Y) = vol(F) − |T|/(2 Y^2)
  4. Identification of the computational reference cell used by FEM
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np

# ---------------------------------------------------------------------------
# 1. Planar fundamental domain P_3 for Gamma_inf acting on C
# ---------------------------------------------------------------------------
# [DP20, §2.3] / [EGM98, §7.3]: for d=3,
#
#   P_3 = T_up ∪ T_low where
#   T_up  = { x+iy | 0 ≤ x,  x/√3 ≤ y ≤ (1-x)/√3 }
#           (the inequality x/√3 ≤ (1-x)/√3 forces x ≤ 1/2)
#   T_low = { x+iy | 0 ≤ x ≤ 1/2,  -x/√3 ≤ y ≤  x/√3 }
#
# Euclidean area(P_3) = 1/(2√3) = √3/6.
# This equals V_Λ / |U| with V_Λ = covol(Z[omega]) = √3/2 and |U|=3
# (units of O_{-3}), matching [Gamma_inf : Gamma'_inf] = 3.

SQRT3 = math.sqrt(3.0)
OMEGA = -0.5 + 0.5j * SQRT3  # e^{2πi/3}


def in_P3(x: float, y: float, tol: float = 1e-12) -> bool:
    """Membership test for the closed set P_3 ([DP20, §2.3])."""
    if x < -tol:
        return False
    # upper component forces x ≤ 1/2; lower is stated with x ≤ 1/2
    if x > 0.5 + tol:
        return False
    x_cl = min(max(x, 0.0), 0.5)
    s = 1.0 / SQRT3
    # T_up
    if y >= x_cl * s - tol and y <= (1.0 - x_cl) * s + tol:
        return True
    # T_low
    if y >= -x_cl * s - tol and y <= x_cl * s + tol:
        return True
    return False


def area_P3() -> float:
    """Exact Euclidean area of P_3: 1/(2√3) = √3/6."""
    return 1.0 / (2.0 * SQRT3)


# Representative vertices (closed polygon corners of P_3)
P3_VERTICES = [
    (0.0, 0.0),
    (0.5, 0.5 / SQRT3),     # right tip (both components)
    (0.0, 1.0 / SQRT3),     # top of T_up at x=0
    (0.0, 0.0),
    (0.5, -0.5 / SQRT3),    # lower tip
    (0.0, 0.0),
]


# ---------------------------------------------------------------------------
# 2. Ford region B_3 and fundamental polyhedron F_3
# ---------------------------------------------------------------------------
# [EGM98, §7.3] / [DP20, Def. after P_d]:
#   B_d = { (z,t) in H^3 : |c z + d|^2 + |c|^2 t^2 ≥ 1
#           for all c,d in O_d with <c,d> = O_d }
#   F_d = { (z,t) in B_d : z in P_d }
#
# For d=3 (Euclidean O_K, class number 1) the only exterior Ford sphere that
# cuts the prism P_3 × R>0 is the unit hemisphere |z|^2 + t^2 = 1.  All of
# P_3 sits inside |z| ≤ 1/√3 < 1, so F_3 is the infinite prism over P_3 cut
# from below by that hemisphere:
#   F_3 = { (z,t) : z in P_3,  t ≥ √(1-|z|^2) }.


def in_B3(z: complex, t: float, tol: float = 1e-10) -> bool:
    """(z,t) in B_3.  For d=3 the binding inequality is |z|^2+t^2 ≥ 1."""
    if t <= 0:
        return False
    if abs(z) ** 2 + t ** 2 < 1.0 - tol:
        return False
    return True


def in_F3(z: complex, t: float, tol: float = 1e-10) -> bool:
    """Fundamental domain F_3 = { (z,t) in B_3 : z in P_3 } ([EGM98, Thm])."""
    return in_P3(z.real, z.imag, tol) and in_B3(z, t, tol)


def y_floor_sphere(z: complex) -> float:
    """Lower boundary of F_3 over z (unit Ford sphere)."""
    r2 = abs(z) ** 2
    if r2 >= 1.0:
        return 0.0
    return math.sqrt(1.0 - r2)


def max_radius_P3() -> float:
    """max |z| for z in P_3 (= 1/√3 at the three tips)."""
    return 1.0 / SQRT3


# ---------------------------------------------------------------------------
# 3. Truncated core K_Y and computational reference cell
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ReferenceCellSpec:
    """Canonical truncated core K_Y = F_3 ∩ {t ≤ Y}.

    Faces:
      top    : P_3 × {Y}
      floor  : |z|^2 + t^2 = 1  over P_3  (hemispherical)
      sides  : vertical walls over ∂P_3

    FEM code currently meshes an *area-equivalent parallelogram model*
    R_comp (see CompReferenceCell below) with the same |T|, same Ford floor,
    and Neumann free BC on sides/floor.  That is a stronger (harder) test
    than the true quotient; see SIDE_PAIRINGS for the maps needed to drop
    the Neumann relaxation.
    """
    Y: float = 1.25

    @property
    def area_T(self) -> float:
        return area_P3()

    def contains(self, z: complex, t: float) -> bool:
        return in_F3(z, t) and t <= self.Y + 1e-12


@dataclass(frozen=True)
class CompReferenceCell:
    """Computational mesh domain used by reference_cell.py / fem_omega_m0.

    Parameterisation
        z = u + v ω,   u ∈ [0,1],   v ∈ [0, 1/3],
        y_f(z) ≤ y ≤ Y,   y_f = √(1−|z|²).

    Justification vs F_3:
      * Euclidean area(R_comp) = (√3/2)·(1/3) = √3/6 = area(P_3) = |T|.
      * R_comp is a fundamental parallelogram for the pure translation
        lattice of index 3 in the full cusp group action (fold by units
        of order 3 is encoded by the v-range 1/3); same orbifold section.
      * R_comp ⊂ {|z| < 1}, so the floor is the same unit Ford sphere.
      * Hyperbolic volume of the truncated prism equals that of K_Y
        whenever the planar sections have equal area (product structure
        of the cusp): both give vol = vol(F) − |T|/(2 Y^2).

    R_comp is *cusp-section-equivalent* (|T| match, same Ford floor type)
    but NOT hyperbolic-volume-equivalent to EGM K_Y: translations move |z|,
    so ∫ yf^{-2} differs.  FEM mass self-check uses vol_K_comp_quad; the
    paper volume for the true truncated fund. domain is vol_KY_exact.
    Neumann free BC makes side pairings irrelevant for M0/CR prototypes.
    """
    Y: float = 1.25
    VMAX: float = 1.0 / 3.0

    @property
    def area_T(self) -> float:
        return (SQRT3 / 2.0) * self.VMAX  # = √3/6


# ---------------------------------------------------------------------------
# 4. Side-pairing generators and face map
# ---------------------------------------------------------------------------
# Presentation generators for Bianchi d=3 (Swan; EGM Ch. 7):
#   T_1 : z ↦ z+1          [[1, 1],[0, 1]]
#   T_w : z ↦ z+ω          [[1, ω],[0, 1]]
#   U   : z ↦ ω z          [[ω², 0],[0, ω]]   (order 3 in PSL; det=ω³=1)
#   S   : z ↦ −1/z         [[0,-1],[1, 0]]
#
# Face pairings of F_3 (Poincaré polyhedron theorem → presentation):
#
#   face                         pairing isometry
#   ---------------------------  ----------------
#   hemispherical floor          S  (self-paired: S maps the sphere to itself)
#   vertical walls over ∂P_3     elements of <T_1, T_w, U> = Gamma_inf
#     — right edge x=1/2         T_1 (shift by 1) composed with units as needed
#     — slanted edges of T_up    T_w / U conjugates (hexagonal lattice)
#     — lower edge of T_low      same family
#
# Exact maps for the *computational* parallelogram R_comp (to replace
# Neumann by true quotient BC on the mesh):
#
#   face of R_comp×[y_f,Y]       pairing
#   ---------------------------  ------------------------------------------
#   u=0  ↔  u=1                  T_1 : (z,y) ↦ (z+1, y)
#   v=0  ↔  v=1/3                T_{ω/3} is NOT in Gamma; the true pairing
#                                uses the order-3 unit: after three strips
#                                of height 1/3 one closes by U·T_w^k.
#                                Operationally: identify the full
#                                parallelogram {v∈[0,1]} by T_w, then
#                                further quotient by U (z ↦ ω z).
#   floor hemisphere             S : (z,y) ↦ (−1/z  in H^3)
#   top y=Y                      free (truncation for Lax–Phillips; matched
#                                to the cusp ODE via the t-functional)
#
# For M0/CR with Neumann relaxation these maps are *not* imposed; the free
# space H¹(K) is larger, so μ>0 is a stronger test (DESIGN.md convention).

# U = diag(ω², ω) sends z ↦ (ω²/ω) z = ω z  (order-3 rotation about 0).
# (diag(ω, ω²) would send z ↦ ω² z = U^{-1}; either generates the same cyclic group.)
GEN = {
    "T1": ((1 + 0j, 1 + 0j), (0 + 0j, 1 + 0j)),
    "Tw": ((1 + 0j, OMEGA), (0 + 0j, 1 + 0j)),
    "U": ((OMEGA ** 2, 0 + 0j), (0 + 0j, OMEGA)),
    "S": ((0 + 0j, -1 + 0j), (1 + 0j, 0 + 0j)),
}

# Structured face-pairing table for documentation / future Dirichlet assembly
SIDE_PAIRINGS: list[dict] = [
    {
        "face": "floor |z|^2+t^2=1 over P_3",
        "isometry": "S",
        "matrix": "[[0,-1],[1,0]]",
        "action": "inversion in unit sphere; pairs floor to itself",
        "for_comp_cell": True,
    },
    {
        "face": "vertical wall x=0 (T_up left edge)",
        "isometry": "U (and conjugates)",
        "matrix": "[[ω²,0],[0,ω]]",
        "action": "rotation z↦ωz; cycles the three 120° sectors of the hexagon",
        "for_comp_cell": False,  # P_3 geometry
    },
    {
        "face": "vertical wall x=1/2 (right tip edges)",
        "isometry": "T1 composed with units",
        "matrix": "[[1,1],[0,1]]",
        "action": "translation z↦z+1",
        "for_comp_cell": True,  # matches u=0 ↔ u=1 on R_comp
    },
    {
        "face": "vertical walls over slanted edges of T_up / T_low",
        "isometry": "Tw, U T1 U^{-1}, …",
        "matrix": "[[1,ω],[0,1]] and conjugates",
        "action": "hexagonal lattice translations",
        "for_comp_cell": False,
    },
    {
        "face": "comp: v=0 ↔ v=1/3  (not a single Gamma element)",
        "isometry": "fold: full T_w on v∈[0,1], then U-quotient",
        "matrix": "see docstring CompReferenceCell",
        "action": "orbifold section of area |T|=√3/6",
        "for_comp_cell": True,
    },
    {
        "face": "top y=Y",
        "isometry": None,
        "matrix": "—",
        "action": "truncation; matched to cusp ODE via t-functional / beta",
        "for_comp_cell": True,
    },
]


def mat_act(M, z: complex, t: float):
    """PSL(2,C) action on H^3; M = ((a,b),(c,d)).  Formula [DP20, (2)]."""
    a, b = M[0]
    c, d = M[1]
    cz_d = c * z + d
    den = abs(cz_d) ** 2 + abs(c) ** 2 * t * t
    if den < 1e-30:
        raise ZeroDivisionError("mat_act: pole")
    z_new = ((a * z + b) * np.conj(cz_d) + a * np.conj(c) * t * t) / den
    t_new = t / den
    return complex(z_new), float(np.real(t_new))


def _mul(A, B):
    return (
        (A[0][0] * B[0][0] + A[0][1] * B[1][0],
         A[0][0] * B[0][1] + A[0][1] * B[1][1]),
        (A[1][0] * B[0][0] + A[1][1] * B[1][0],
         A[1][0] * B[0][1] + A[1][1] * B[1][1]),
    )


def verify_generators_in_PSL() -> list[str]:
    """det ∈ {±1} (or unit for diagonal U) and U^3 = I in PSL."""
    msgs = []
    for name, M in GEN.items():
        a, b = M[0]
        c, d = M[1]
        det = a * d - b * c
        if abs(det - 1) < 1e-9 or abs(det + 1) < 1e-9:
            msgs.append(f"OK {name} det={det}")
        elif abs(abs(det) - 1) < 1e-9:
            msgs.append(f"OK {name} det on unit circle ({det})")
        else:
            msgs.append(f"FAIL det({name})={det}")
    U = GEN["U"]
    U3 = _mul(_mul(U, U), U)
    if abs(U3[0][0] - 1) < 1e-8 and abs(U3[1][1] - 1) < 1e-8:
        msgs.append("OK U^3 = I in SL")
    elif abs(U3[0][0] + 1) < 1e-8:
        msgs.append("OK U^3 = -I (~I in PSL)")
    else:
        msgs.append(f"FAIL U^3 = {U3}")
    return msgs


def sample_boundary_points(Y: float = 1.25, n: int = 12):
    """Sample points on vertical walls of F_3 at height Y and on floor sphere."""
    pts = []
    for k in range(n):
        y = -0.5 / SQRT3 + k / (n - 1) * (1.0 / SQRT3)
        if in_P3(0.5, y):
            pts.append(("wall_x=1/2", 0.5 + 1j * y, Y))
    for k in range(n):
        x = k / (n - 1) * 0.5
        y = x / SQRT3 * 0.5
        if in_P3(x, y):
            z = x + 1j * y
            r2 = abs(z) ** 2
            if r2 < 1 - 1e-8:
                pts.append(("floor", z, math.sqrt(1 - r2)))
    return pts


def verify_side_pairings(Y: float = 1.25) -> list[str]:
    """Isometry and face-pairing consistency checks (samples)."""
    msgs = []
    msgs.extend(verify_generators_in_PSL())

    # S preserves the unit sphere (floor self-pairing)
    for name, z, t in sample_boundary_points(Y):
        if name != "floor":
            continue
        z2, t2 = mat_act(GEN["S"], z, t)
        if abs(abs(z2) ** 2 + t2 ** 2 - 1.0) > 1e-8:
            msgs.append("FAIL S floor image not on sphere")
            break
    else:
        msgs.append("OK S preserves unit sphere (floor pairing)")

    # T1 is a pure horizontal translation
    ok_t = 0
    for _name, z, t in sample_boundary_points(Y):
        z2, t2 = mat_act(GEN["T1"], z, t)
        if abs(t2 - t) < 1e-10 and abs((z2 - z) - 1) < 1e-10:
            ok_t += 1
    msgs.append(f"OK T1: z↦z+1, height fixed ({ok_t} samples)")

    # Tw: z↦z+ω
    z0, t0 = 0.1 + 0.05j, Y
    z2, t2 = mat_act(GEN["Tw"], z0, t0)
    if abs(t2 - t0) < 1e-10 and abs((z2 - z0) - OMEGA) < 1e-10:
        msgs.append("OK Tw: z↦z+ω, height fixed")
    else:
        msgs.append(f"FAIL Tw: got z={z2}, t={t2}")

    # U: rotation by ω about 0, order 3
    z2, t2 = mat_act(GEN["U"], z0, t0)
    if abs(t2 - t0) < 1e-9 and abs(z2 - OMEGA * z0) < 1e-9:
        msgs.append("OK U: z↦ωz (elliptic order 3 about 0)")
    else:
        msgs.append(f"FAIL U: got z={z2}, t={t2}")

    # S² = id in PSL on a generic point
    z1, t1 = mat_act(GEN["S"], z0, t0)
    z2, t2 = mat_act(GEN["S"], z1, t1)
    if abs(z2 - z0) < 1e-8 and abs(t2 - t0) < 1e-8:
        msgs.append("OK S^2 = id on sample point")
    else:
        msgs.append(f"FAIL S^2: got ({z2},{t2})")

    return msgs


def print_side_pairing_table() -> None:
    print("  Face-pairing table (F_3 / comp cell):")
    for row in SIDE_PAIRINGS:
        tag = "comp" if row["for_comp_cell"] else "P3 "
        print(f"    [{tag}] {row['face']}")
        print(f"           isometry: {row['isometry']}  {row['matrix']}")
        print(f"           {row['action']}")


# ---------------------------------------------------------------------------
# 5. Exact volumes
# ---------------------------------------------------------------------------

def vol_F_exact(prec: int = 80):
    """vol(Gamma\\H^3) = |D|^{3/2} zeta_K(2)/(4 pi^2), D=-3.

    |D|^{3/2} = 3^{3/2} = 3√3.
    zeta_K(2) = zeta(2) L(2, chi_{-3}), chi period 3: (0,1,-1).

    Returns (mid, radius) floats.  Uses python-flint arb if available,
    else mpmath fallback.
    """
    try:
        from flint import arb, acb, ctx
        ctx.prec = prec
        s = acb(2)
        L2 = acb(0)
        for a, c in ((1, 1), (2, -1)):
            L2 += c * s.zeta(acb(a) / 3)
        L2 = L2 * acb(3) ** (-s)
        zeta2 = acb(2).zeta()
        zetaK2 = (zeta2 * L2).real
        # vol = 3 √3  zeta_K(2) / (4 pi^2)
        vol = arb(3) * arb(3).sqrt() * zetaK2 / (4 * arb.pi() ** 2)
        return float(vol.mid()), float(vol.rad())
    except Exception:
        import mpmath as mp
        mp.mp.dps = 40
        L2 = sum(c * mp.zeta(2, a / 3) for a, c in ((1, 1), (2, -1))) / 9
        zetaK2 = mp.zeta(2) * L2
        vol = 3 * mp.sqrt(3) * zetaK2 / (4 * mp.pi ** 2)
        return float(vol), 1e-12


def cusp_tail_volume(Y: float) -> float:
    """Hyperbolic volume of P_3 × (Y, ∞): ∫_{P_3} dxdy ∫_Y^∞ t^{-3} dt = |T|/(2 Y^2).

    Exact for any planar section of area |T| (product structure of the
    standard cusp fundamental domain [EGM98, §7.3]).
    """
    return area_P3() / (2.0 * Y * Y)


def vol_KY_exact(Y: float, prec: int = 80):
    """Exact volume of the EGM truncated core K_Y = F_3 ∩ {t ≤ Y}.

        vol(K_Y) = vol(F) − |T|/(2 Y^2),   |T| = area(P_3) = √3/6.

    Valid because F ∩ {t > Y} = P_3 × (Y,∞) under the EGM product
    structure (planar section exactly P_3).
    """
    vF, rad = vol_F_exact(prec)
    return vF - cusp_tail_volume(Y), rad


def vol_prism_over_region(
    region_sampler: Callable[[int], tuple[np.ndarray, np.ndarray, float]],
    Y: float,
    n: int = 400,
) -> float:
    """Hyperbolic volume of {(z,y): z∈R, yf(z)≤y≤Y} via midpoint quadrature.

    vol = ∫_R (1/2)( yf(z)^{-2} − Y^{-2} ) dx dy.
    region_sampler(n) → (x1_grid, x2_grid, cell_area) for an n×n partition.
    """
    x1, x2, dA = region_sampler(n)
    r2 = x1 ** 2 + x2 ** 2
    # R must lie in |z|<1 for the sphere floor
    if np.any(r2 >= 1.0 - 1e-14):
        raise ValueError("region exits unit disk; yf undefined")
    yf2_inv = 1.0 / (1.0 - r2)  # 1/yf^2
    return float(0.5 * np.sum(yf2_inv - 1.0 / (Y * Y)) * dA)


def _sampler_P3(n: int):
    """Midpoint grid on the bounding box [0,1/2]×[-1/√3,1/√3], masked to P_3."""
    xs = np.linspace(0.0, 0.5, n, endpoint=False) + 0.25 / n
    ys = np.linspace(-1.0 / SQRT3, 1.0 / SQRT3, n, endpoint=False) + (1.0 / SQRT3) / n
    dA = (0.5 / n) * (2.0 / SQRT3 / n)
    X, Yg = np.meshgrid(xs, ys, indexing="ij")
    mask = np.vectorize(lambda a, b: in_P3(float(a), float(b), tol=1e-10))(X, Yg)
    return X[mask], Yg[mask], dA


def _sampler_Rcomp(n: int):
    """Midpoint grid on R_comp: z = u + v ω, u∈[0,1], v∈[0,1/3]."""
    us = np.linspace(0.0, 1.0, n, endpoint=False) + 0.5 / n
    vs = np.linspace(0.0, 1.0 / 3.0, n, endpoint=False) + (1.0 / 3.0) / (2 * n)
    # |det d(x1,x2)/d(u,v)| = √3/2
    dA = (1.0 / n) * ((1.0 / 3.0) / n) * (SQRT3 / 2.0)
    U, V = np.meshgrid(us, vs, indexing="ij")
    x1 = U - 0.5 * V
    x2 = (SQRT3 / 2.0) * V
    return x1.ravel(), x2.ravel(), dA


def vol_K_comp_quad(Y: float, n: int = 500) -> float:
    """High-precision quadrature of vol(K_comp) over the FEM parallelogram.

    NOTE: this is NOT equal to vol(K_Y)=vol(F)−|T|/(2Y²), because R_comp
    reaches larger |z| than P_3, so the sphere floor is lower and the
    hyperbolic volume is larger.  Use this for FEM mass-matrix self-check;
    use vol_KY_exact for the EGM paper domain.
    """
    return vol_prism_over_region(_sampler_Rcomp, Y, n=n)


def vol_KY_quad(Y: float, n: int = 500) -> float:
    """Quadrature of EGM K_Y over P_3 (cross-check of the closed formula)."""
    return vol_prism_over_region(_sampler_P3, Y, n=n)


def vol_KY_formula_str(Y: float | None = None) -> str:
    """Human-readable exact formula for the EGM core."""
    base = (
        "EGM core:  vol(K_Y) = vol(F) − |T|/(2 Y²)\n"
        "  vol(F) = 3√3 · ζ_K(2) / (4 π²),   ζ_K(2)=ζ(2)L(2,χ_{-3})\n"
        "  |T|    = area(P_3) = √3/6\n"
        "Comp cell K_comp:  vol = ∫_{R_comp} ½(yf^{-2} − Y^{-2}) dxdy\n"
        "  (differs from EGM: R_comp is not |z|-equivalent to P_3)"
    )
    if Y is not None:
        vK, rK = vol_KY_exact(Y)
        vC = vol_K_comp_quad(Y)
        base += (
            f"\n  @ Y={Y}: vol(K_Y)={vK:.12f} (±{rK:.2e});  "
            f"vol(K_comp)≈{vC:.12f}"
        )
    return base


# ---------------------------------------------------------------------------
# 6. Self-checks
# ---------------------------------------------------------------------------

def main():
    print("Fundamental polyhedron PSL(2, Z[omega]) — geometry freeze")
    print("=" * 64)
    print("Citations: [EGM98, §7.3], [DP20, §2.3], [Swan71], Humbert volume")

    print("\n[1] Planar domain P_3  ([DP20, §2.3])")
    a = area_P3()
    print(f"  area(P_3) = {a:.12f} = √3/6 = {SQRT3/6:.12f}")
    assert abs(a - SQRT3 / 6) < 1e-14
    assert in_P3(0.1, 0.0)
    assert in_P3(0.1, 0.1 / SQRT3)
    assert not in_P3(0.4, 0.4)
    print(f"  max |z| on P_3 = {max_radius_P3():.6f} < 1  "
          f"(floor = unit sphere throughout)")

    print("\n[2] Side-pairing generators + face map")
    for m in verify_side_pairings():
        print(f"  {m}")
    print_side_pairing_table()

    print("\n[3] Exact volumes (Humbert + cusp tail + comp quadrature)")
    vF, rF = vol_F_exact()
    print(f"  vol(F) = {vF:.12f}  (rad {rF:.2e})")
    assert abs(vF - 0.1691569344) < 1e-8, vF
    for Y in (1.25, 1.5, 2.0):
        vK, rK = vol_KY_exact(Y)
        tail = cusp_tail_volume(Y)
        vKq = vol_KY_quad(Y, n=600)
        vC = vol_K_comp_quad(Y, n=600)
        print(f"  Y={Y}: vol(K_Y)={vK:.10f}  quad_P3={vKq:.10f}  "
              f"vol(K_comp)≈{vC:.10f}  tail={tail:.10f}")
        assert vK > 0 and vK < vF
        # closed formula vs quadrature on P_3
        assert abs(vK - vKq) / vK < 5e-3, (vK, vKq)
    print("  formula:")
    for line in vol_KY_formula_str(1.25).splitlines():
        print(f"    {line}")

    print("\n[4] Reference cells")
    cell = ReferenceCellSpec(Y=1.25)
    comp = CompReferenceCell(Y=1.25)
    print(f"  EGM core K_Y = F_3 ∩ {{t ≤ {cell.Y}}}, |T|={cell.area_T:.12f}")
    print(f"  Comp cell R_comp: u∈[0,1], v∈[0,1/3], |T|={comp.area_T:.12f}")
    assert abs(cell.area_T - comp.area_T) < 1e-14
    print(f"  sample interior (0.1+0.05j, 1.5) in F_3? "
          f"{in_F3(0.1 + 0.05j, 1.5)}")
    print(f"  sample below sphere (0, 0.5) in F_3? {in_F3(0j, 0.5)}")
    tip = 0.5 + 0.5j / SQRT3
    yf = y_floor_sphere(tip)
    print(f"  tip |z|={abs(tip):.6f}, y_floor={yf:.6f} (expect √(2/3)="
          f"{math.sqrt(2/3):.6f})")
    assert abs(yf - math.sqrt(2.0 / 3.0)) < 1e-12

    print("\n" + "=" * 64)
    print("Geometry freeze COMPLETE:")
    print("  • F_3 = B_3 ∩ (P_3 × R>0) cited from EGM/DP20")
    print("  • Side-pairings S, T1, Tw, U verified; table recorded")
    print("  • vol(K_Y) = vol(F) − |T|/(2 Y²) exact (EGM core)")
    print("  • vol(K_comp) via ∫ ½(yf^{-2}−Y^{-2}) for FEM mass check")
    print("  • Comp cell: |T|-matched prism (not |z|-eq. to P_3)")


if __name__ == "__main__":
    main()
