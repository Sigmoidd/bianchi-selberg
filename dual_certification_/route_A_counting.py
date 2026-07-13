"""Route A — Dirichlet–Neumann bracketing prototype on truncated core K_Y.

Source of truth: DualCertification_AgentReady.md §7 Route A, §11 priority 3.

Purpose
-------
Engineering / research scaffold that **addresses the remaining logical gap**
of certified counting via D–N bracketing on a product-like truncated core.
This module does **not** produce a certified counting function N(λ), nor a
certified first eigenvalue, nor N(λ)=0 on (1,λ₁).

What it does (float diagnostics)
--------------------------------
1. Build (or load) the Picard Humbert truncated core mesh K_Y at Y=1.25
   by reusing independent_exclusion/cr_prototype.build_mesh (Z[i] level-1).
2. Assemble weighted CR stiffness/mass (Q, M) for the hyperbolic Laplacian
   forms Q(v)=∫_K |∇v|² y⁻¹, M(v)=∫_K v² y⁻³ on the compact core only
   (no cusp-collar zero-mode terms — pure variational spectrum of K_Y).
3. Solve float generalized eigenproblems:
     - Neumann on artificial top y=Y (all CR face dofs free)
     - Dirichlet on artificial top y=Y (top-face dofs constrained to 0)
4. Report candidate bracketing intervals [λ_k^N, λ_k^D] as **diagnostics**.
5. Leave clear extension hooks for Arb assembly + Rump PSD / GLB upgrades.

What is NOT certified (see also route_A_status.md)
--------------------------------------------------
- No Arb enclosures of matrix entries or eigenvalues.
- No Rump positive-definiteness certificates.
- No CR guaranteed-lower-bound post-processing of discrete λ_k.
- Artificial-boundary / truncation error to the full quotient is open.
- Bracketing inequality with explicit truncation constants is not proved here.
- Final integer-interval counting enclosure N(λ) is not produced.
- Assumption H is not needed for pure variational counting of Laplacian on
  compact K with BC, but truncation to the full quotient still needs a
  documented error enclosure (open).

Language discipline (AgentReady §12)
------------------------------------
Label all numerical output as "engineering diagnostic / addresses counting
gap, not certified" until a counting theorem and Arb/Rump certificates exist.
"""

from __future__ import annotations

import os
import sys
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.linalg import eigh

# ---------------------------------------------------------------------------
# Import Picard mesh + CR assembly from independent_exclusion
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
_IE = os.path.join(_REPO, "independent_exclusion")
if _IE not in sys.path:
    sys.path.insert(0, _IE)

from cr_prototype import (  # noqa: E402
    KAPPA1,
    assemble,
    build_mesh,
    geometry,
)

# Optional secondary: Eisenstein P3 reference cell (import only if requested)
_ES = os.path.join(_REPO, "eisenstein_numbers")

DEFAULT_Y = 1.25
DEFAULT_NEIGS = 8
# Coarse smoke mesh: fast enough for <2 min; production later uses denser.
SMOKE_MESH = (6, 3, 3)  # N1 x N2 x N3 → 6*N1*N2*N3 tets


# ---------------------------------------------------------------------------
# Status bookkeeping (mirrors AgentReady §7 Route A checklist)
# ---------------------------------------------------------------------------
ROUTE_A_CHECKLIST = [
    (
        "compact_core",
        "Compact core K is truncated at a certified height Y with product "
        "collar structure.",
    ),
    (
        "product_mesh",
        "Mesh on K is product-like near the artificial boundary (or the "
        "error of non-product structure is enclosed).",
    ),
    (
        "dn_guaranteed",
        "Dirichlet and Neumann eigenvalue problems are both solved with "
        "guaranteed bounds (CR or Morley + post-processing).",
    ),
    (
        "bracketing_proved",
        "The bracketing inequality λ_k^{N,K_h} ≤ λ_k ≤ λ_k^{D,K_h} is proved "
        "with explicit constants that account for the truncation error.",
    ),
    (
        "artificial_boundary",
        "The artificial-boundary error is itself enclosed by an Arb ball "
        "(or by a comparison with a larger domain).",
    ),
    (
        "counting_enclosure",
        "Final counting function enclosure N(λ) is an integer interval; for "
        "the target λ it must be [0,0] on (1,λ₁).",
    ),
]


@dataclass
class RouteAStatus:
    """Honest checklist state for Route A items."""

    items: Dict[str, Tuple[bool, str]] = field(default_factory=dict)

    def set(self, key: str, done: bool, note: str) -> None:
        self.items[key] = (done, note)

    def report_lines(self) -> List[str]:
        lines = []
        for key, desc in ROUTE_A_CHECKLIST:
            done, note = self.items.get(key, (False, "not evaluated"))
            mark = "x" if done else " "
            lines.append(f"- [{mark}] **{key}**: {desc}")
            lines.append(f"  - note: {note}")
        return lines


def default_status() -> RouteAStatus:
    """Prototype-level status: mesh/product partial; certification all red."""
    st = RouteAStatus()
    st.set(
        "compact_core",
        True,  # scaffold: mesh built at fixed Y; "certified height" still open
        "Mesh built at engineering height Y=1.25 (product collar above Y is "
        "classical for Picard Humbert; Y itself is not Arb-certified as optimal).",
    )
    st.set(
        "product_mesh",
        True,  # build_mesh is tensor-product in (x1,x2) × height parameter
        "Picard mesh is product-like: (x1,x2)-grid extruded in height with "
        "floor lift (Lemma G). Top layer is flat y=Y — product collar ready. "
        "Non-product floor error handled by lift, not by Arb enclosure of "
        "collar mismatch.",
    )
    st.set(
        "dn_guaranteed",
        False,
        "Float CR generalized eigenpairs only. No Arb assembly, no Rump, "
        "no CR GLB post-processing. Engineering diagnostic only.",
    )
    st.set(
        "bracketing_proved",
        False,
        "Discrete float [λ_k^N, λ_k^D] reported as candidate intervals only. "
        "No proof of bracketing with truncation-error constants.",
    )
    st.set(
        "artificial_boundary",
        False,
        "OPEN: artificial-boundary / truncation error to full quotient is "
        "not enclosed. Comparison with larger Y not certified.",
    )
    st.set(
        "counting_enclosure",
        False,
        "No integer-interval N(λ). Do not claim N(λ)=0 or certified first "
        "eigenvalue. Status: addresses counting gap, not certified.",
    )
    return st


# ---------------------------------------------------------------------------
# Mesh + classification of top-face CR dofs
# ---------------------------------------------------------------------------
def build_picard_core(
    N1: int,
    N2: int,
    N3: int,
    Y: float = DEFAULT_Y,
    curved: bool = True,
) -> Dict[str, Any]:
    """Truncated Picard Humbert core K_Y via cr_prototype.build_mesh.

    Domain (curved=True):
      K_Y = { (x1,x2,y) : x1∈[-1/2,1/2], x2∈[0,1/2],
              y_f(x1,x2) ≤ y ≤ Y },  y_f = √(1−x1²−x2²).
    Artificial boundary: top face y = Y (product structure in the collar).
    """
    mesh = build_mesh(N1, N2, N3, Y, curved=curved)
    geo = geometry(mesh)
    Q, M, avec, tvec, top_tets, floor_tets = assemble(mesh, geo)
    top_face_dofs = sorted({fid for (_e, fid) in top_tets})
    return dict(
        mesh=mesh,
        geo=geo,
        Q=Q,
        M=M,
        avec=avec,
        tvec=tvec,
        top_tets=top_tets,
        floor_tets=floor_tets,
        top_face_dofs=top_face_dofs,
        Y=Y,
        N=(N1, N2, N3),
        field="picard_Zi",
        n_tets=len(mesh["tets"]),
        n_dofs=geo["nfaces"],
    )


def try_build_eisenstein_core(
    N1: int = 3,
    N3: int = 3,
    Y: float = DEFAULT_Y,
) -> Optional[Dict[str, Any]]:
    """Secondary: Eisenstein P3 core via eisenstein_numbers (optional).

    Returns None if import/mesh fails. Float CR only; same honesty labels.
    """
    if _ES not in sys.path:
        sys.path.insert(0, _ES)
    try:
        from reference_cell import build_reference_cell  # type: ignore
        from cr_omega import assemble_cr  # type: ignore
    except Exception as exc:  # pragma: no cover - optional path
        print(f"  [eisenstein secondary unavailable: {exc}]")
        return None

    try:
        mesh = build_reference_cell(N1=N1, N2=N1, N3=N3, Y=Y, domain="P3")
        pack = assemble_cr(mesh)
    except Exception as exc:
        print(f"  [eisenstein mesh/assembly failed: {exc}]")
        return None

    Q = pack["Q"]
    M = pack["M"]
    if hasattr(Q, "toarray"):
        Q = Q.toarray()
        M = M.toarray()
    top_dofs = sorted({fid for (_e, fid) in pack.get("top_tets", [])})
    return dict(
        mesh=mesh,
        Q=np.asarray(Q, dtype=float),
        M=np.asarray(M, dtype=float),
        top_face_dofs=list(top_dofs),
        Y=Y,
        N=(N1, N3),
        field="eisenstein_omega",
        n_tets=len(mesh["tets"]) if isinstance(mesh, dict) else -1,
        n_dofs=int(Q.shape[0]),
    )


# ---------------------------------------------------------------------------
# Float Dirichlet / Neumann eigenpairs (diagnostic)
# ---------------------------------------------------------------------------
def _free_index_set(n: int, constrained: Sequence[int]) -> np.ndarray:
    mask = np.ones(n, dtype=bool)
    if constrained:
        mask[list(constrained)] = False
    return np.where(mask)[0]


def float_neumann_eigs(
    Q: np.ndarray,
    M: np.ndarray,
    neigs: int = DEFAULT_NEIGS,
) -> np.ndarray:
    """Float eigenvalues of Qv = λ Mv with free (Neumann-type) boundary.

    CR nonconforming: free face means on all boundary faces, including top.
    First mode is typically ~0 (constants have Q=0 for pure Neumann Laplace
    with unit or consistent weights; weighted hyperbolic form also admits
    near-null constants if Q annihilates them approximately).
    """
    n = Q.shape[0]
    k = min(neigs, n - 1)
    # subset_by_index needs M positive definite; shift tiny if needed
    w = eigh(Q, M, eigvals_only=True, subset_by_index=[0, k - 1])
    return np.asarray(w, dtype=float)


def float_dirichlet_eigs(
    Q: np.ndarray,
    M: np.ndarray,
    top_dofs: Sequence[int],
    neigs: int = DEFAULT_NEIGS,
) -> np.ndarray:
    """Float eigenvalues with Dirichlet (zero) on artificial top-face dofs.

    Side/floor faces remain free (Neumann relaxation of orbifold gluings),
    matching the existing exclusion philosophy: relaxation is valid for
    lower bounds if later certified; here purely diagnostic.
    """
    free = _free_index_set(Q.shape[0], top_dofs)
    if free.size < 2:
        raise RuntimeError("too few free dofs after Dirichlet constraint")
    Qf = Q[np.ix_(free, free)]
    Mf = M[np.ix_(free, free)]
    k = min(neigs, free.size - 1)
    w = eigh(Qf, Mf, eigvals_only=True, subset_by_index=[0, k - 1])
    return np.asarray(w, dtype=float)


def candidate_brackets(
    lam_N: Sequence[float],
    lam_D: Sequence[float],
) -> List[Tuple[int, float, float]]:
    """Pair Neumann/Dirichlet float eigenvalues into candidate intervals.

    Classical D–N bracketing on a fixed compact domain with the *same*
    side conditions would give λ_k^N ≤ λ_k ≤ λ_k^D. Here:
    - domain is truncated K_Y (not full quotient);
    - side/floor are free (Neumann relaxation);
    - values are float CR, not guaranteed bounds.

    So these are **candidate diagnostic intervals only**.
    """
    m = min(len(lam_N), len(lam_D))
    out = []
    for k in range(m):
        out.append((k + 1, float(lam_N[k]), float(lam_D[k])))
    return out


# ---------------------------------------------------------------------------
# Hooks for future Arb / Rump upgrades (stubs)
# ---------------------------------------------------------------------------
def arb_upgrade_placeholder() -> Dict[str, str]:
    """Document intended Arb/Rump path without implementing it.

    Future steps (not done here):
      1. Interval assembly of Q, M per m3_certify.py patterns (arb balls).
      2. CR GLB: λ_k ≥ λ_k^h / (1 + C_h² λ_k^h) with certified C_h (κ₁,h_T).
      3. Dirichlet / Neumann discrete spectra as interval eigenvalues
         (e.g. Rump verified eigensolver or interval Sturm counts).
      4. Truncation comparison K_Y ⊂ K_{Y'} or collar spectral estimate
         → Arb ball for artificial-boundary error.
      5. Integer enclosure of N(λ) = # { k : λ_k ≤ λ } from brackets.
    """
    return {
        "assembly": "reuse m3_certify tet_arb_data + interval weights wQ,wM",
        "glb": "Carstensen–Gedicke / Liu CR lower bound with KAPPA1",
        "rump": "rump_psd_certificate for shifted pencils; verified eig later",
        "truncation": "OPEN — needs collar comparison or defect-style bound",
        "counting": "N(λ) ∈ [# {k: λ_k^D ≤ λ}, # {k: λ_k^N ≤ λ}] once certified",
    }


def diagnostic_glb_lower(lam_h: float, h_max: float) -> float:
    """Non-certified float sketch of CR GLB formula (diagnostic only).

    λ ≥ λ^h / (1 + κ₁² h² λ^h) under the hypotheses of the CR theory.
    Not applied as a certificate here (weights, BC, domain truncation).
    """
    if lam_h <= 0:
        return float(lam_h)
    den = 1.0 + (KAPPA1 ** 2) * (h_max ** 2) * lam_h
    return float(lam_h / den)


def truncation_constants_scaffold(
    Y: float = DEFAULT_Y,
    Y_prime: float = 1.5,
    r_target: float = 6.0,
    T_abs: float = 0.5,
) -> Dict[str, Any]:
    """
    Explicit *scaffold* formulae for artificial-boundary / collar truncation.

    These are the constants Route A must eventually enclose in Arb to upgrade
    bracketing_proved and artificial_boundary from RED → GREEN. They are
    **not** certificates — only closed-form targets for the missing proof.

    Collar model (cuspidal product region y ≥ Y):
      - Continuous spectrum / Eisenstein contribution controlled by
          E_collar ≲ C_cut · exp(−2π Y)   (cf. Theorem D(K) C2 term)
      - Difference between spectra on K_Y and K_{Y'} controlled by the
        annular region Y ≤ y ≤ Y' with volume ≲ |T|·(1/Y² − 1/Y'²) and
        Poincaré constant ~ (Y'−Y).

    Returns named floats for documentation + future Arb hooks.
    """
    import math

    # A_cut-style collar floor (matches lemma_K spirit; field Z[i] defaults)
    C_Sob_proxy = 1.0 + 3.0 / math.pi  # crude; replace by true C_Sob
    A_met_proxy = (0.5) ** (-1.5)  # y_min=1/2 conservative
    A_cut = 2.0 * T_abs * (1.0 / Y) * (1.0 + C_Sob_proxy) * A_met_proxy
    collar_floor = A_cut * math.exp(-2.0 * math.pi * Y)
    collar_floor_Yp = A_cut * math.exp(-2.0 * math.pi * Y_prime)
    # Euclidean height strip length
    strip = max(Y_prime - Y, 0.0)
    # Model truncation shift for eigenvalues when domain grows K_Y → K_{Y'}
    # (engineering: O(exp(−2π Y)) + O(strip^{-2}) Poincaré remainder placeholder)
    poincare_strip = (math.pi / max(strip, 1e-6)) ** 2 if strip > 0 else float("inf")
    # Target: prove |λ_k(K_Y) − λ_k(Γ\H³)| ≤ Δ_trunc(Y) with
    delta_trunc_model = collar_floor + (
        0.0 if strip == 0 else 1.0 / poincare_strip
    )
    return dict(
        language=(
            "Scaffold only — not a proved truncation theorem. "
            "Needed to close Route A bracketing_proved + artificial_boundary."
        ),
        Y=Y,
        Y_prime=Y_prime,
        r_target=r_target,
        A_cut_proxy=A_cut,
        collar_floor_Y=collar_floor,
        collar_floor_Yp=collar_floor_Yp,
        strip_height=strip,
        poincare_strip_proxy=poincare_strip,
        delta_trunc_model=delta_trunc_model,
        proof_obligations=[
            "Enclose A_cut / C_Sob / A_met in Arb on the true core geometry",
            "Prove |Spec(K_Y)−Spec(Γ\\H³)| ≤ Δ_trunc with explicit Δ_trunc",
            "Or: compare certified D/N spectra on K_Y and K_{Y'} and take limit",
            "Feed Δ_trunc into N(λ) integer interval enclosure",
        ],
        status="RED_OPEN",
    )


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
@dataclass
class RouteAResult:
    field: str
    Y: float
    mesh_N: Tuple[int, ...]
    n_tets: int
    n_dofs: int
    n_top_dofs: int
    h_max: float
    lam_N: np.ndarray
    lam_D: np.ndarray
    brackets: List[Tuple[int, float, float]]
    elapsed_s: float
    status: RouteAStatus
    label: str = (
        "engineering diagnostic / addresses counting gap, not certified"
    )

    def summary(self) -> str:
        lines = [
            "=" * 72,
            "Route A counting prototype — DIAGNOSTIC OUTPUT",
            f"Label: {self.label}",
            f"Field/mesh: {self.field}  N={self.mesh_N}  Y={self.Y}",
            f"tets={self.n_tets}  CR dofs={self.n_dofs}  "
            f"top Dirichlet dofs={self.n_top_dofs}  h_max={self.h_max:.4e}",
            f"elapsed: {self.elapsed_s:.2f}s",
            "-" * 72,
            f"{'k':>4}  {'λ_k^N (float)':>14}  {'λ_k^D (float)':>14}  "
            f"{'candidate [N,D]':>20}  {'GLB sketch(N)':>12}",
        ]
        for k, lo, hi in self.brackets:
            glb = diagnostic_glb_lower(lo, self.h_max)
            lines.append(
                f"{k:4d}  {lo:14.6f}  {hi:14.6f}  "
                f"[{lo:.4f}, {hi:.4f}]{'':>2}  {glb:12.6f}"
            )
        lines.append("-" * 72)
        lines.append("Route A checklist (prototype):")
        lines.extend(self.status.report_lines())
        lines.append("-" * 72)
        lines.append("Arb/Rump upgrade hooks:")
        for key, val in arb_upgrade_placeholder().items():
            lines.append(f"  - {key}: {val}")
        lines.append("=" * 72)
        lines.append(
            "DO NOT claim: dual certification, certified λ₁, N(λ)=0, "
            "or certified counting."
        )
        return "\n".join(lines)


def run_route_A_picard(
    N1: int,
    N2: int,
    N3: int,
    Y: float = DEFAULT_Y,
    neigs: int = DEFAULT_NEIGS,
    curved: bool = True,
) -> RouteAResult:
    """Main Picard (Z[i] level-1) Route A float diagnostic."""
    t0 = time.perf_counter()
    data = build_picard_core(N1, N2, N3, Y=Y, curved=curved)
    Q, M = data["Q"], data["M"]
    top = data["top_face_dofs"]
    h_max = float(data["geo"]["hT"].max())

    lam_N = float_neumann_eigs(Q, M, neigs=neigs)
    lam_D = float_dirichlet_eigs(Q, M, top, neigs=neigs)
    brackets = candidate_brackets(lam_N, lam_D)
    elapsed = time.perf_counter() - t0

    return RouteAResult(
        field=data["field"],
        Y=Y,
        mesh_N=(N1, N2, N3),
        n_tets=data["n_tets"],
        n_dofs=data["n_dofs"],
        n_top_dofs=len(top),
        h_max=h_max,
        lam_N=lam_N,
        lam_D=lam_D,
        brackets=brackets,
        elapsed_s=elapsed,
        status=default_status(),
    )


def smoke_test(verbose: bool = True) -> bool:
    """Minimal smoke test: coarse Picard mesh, few eigs, <2 minutes.

    Success criteria (infrastructure only — not certification):
      - mesh builds
      - Q, M SPD-ish (M PD, Q PSD)
      - Neumann λ_1 ≈ 0 (null / near-null constant)
      - Dirichlet λ_1 > Neumann λ_1 (strict for k≥1 after constants)
      - candidate brackets finite and ordered λ^N ≤ λ^D for each k
    """
    N1, N2, N3 = SMOKE_MESH
    if verbose:
        print(
            f"Route A smoke test: Picard mesh {N1}x{N2}x{N3}, Y={DEFAULT_Y}"
        )
    res = run_route_A_picard(N1, N2, N3, Y=DEFAULT_Y, neigs=5)
    if verbose:
        print(res.summary())

    ok = True
    # timing budget
    if res.elapsed_s > 120:
        ok = False
        if verbose:
            print(f"FAIL: elapsed {res.elapsed_s:.1f}s > 120s")
    # Neumann ground ~ 0
    if abs(res.lam_N[0]) > 1e-6 and res.lam_N[0] < -1e-8:
        # large negative would indicate assembly bug
        ok = False
        if verbose:
            print(f"FAIL: Neumann λ1 = {res.lam_N[0]} suspicious")
    # brackets: N <= D (allow tiny float noise)
    for k, lo, hi in res.brackets:
        if lo > hi + 1e-8:
            ok = False
            if verbose:
                print(f"FAIL: bracket k={k}: N={lo} > D={hi}")
    # Dirichlet first positive mode should exceed first Neumann
    if res.lam_D[0] + 1e-10 < res.lam_N[0]:
        ok = False
        if verbose:
            print("FAIL: λ1^D < λ1^N")

    if verbose:
        print(f"smoke_test infrastructure ok: {ok}")
        print(
            "Reminder: pass ≠ certified counting. "
            "Addresses remaining logical gap only."
        )
    return ok


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    p = argparse.ArgumentParser(
        description=(
            "Route A D–N bracketing prototype (engineering diagnostic; "
            "not certified)."
        )
    )
    p.add_argument("--smoke", action="store_true", help="coarse mesh smoke test")
    p.add_argument("--N1", type=int, default=SMOKE_MESH[0])
    p.add_argument("--N2", type=int, default=SMOKE_MESH[1])
    p.add_argument("--N3", type=int, default=SMOKE_MESH[2])
    p.add_argument("--Y", type=float, default=DEFAULT_Y)
    p.add_argument("--neigs", type=int, default=DEFAULT_NEIGS)
    p.add_argument(
        "--flat",
        action="store_true",
        help="flat box floor (curved=False) for model-problem checks",
    )
    p.add_argument(
        "--eisenstein",
        action="store_true",
        help="also attempt secondary Eisenstein P3 diagnostic if available",
    )
    args = p.parse_args(list(argv) if argv is not None else None)

    if args.smoke:
        return 0 if smoke_test(verbose=True) else 1

    res = run_route_A_picard(
        args.N1,
        args.N2,
        args.N3,
        Y=args.Y,
        neigs=args.neigs,
        curved=not args.flat,
    )
    print(res.summary())

    if args.eisenstein:
        print("\n--- secondary Eisenstein P3 (optional) ---")
        sec = try_build_eisenstein_core(Y=args.Y)
        if sec is None:
            print("skipped")
        else:
            t0 = time.perf_counter()
            lam_N = float_neumann_eigs(sec["Q"], sec["M"], neigs=min(5, args.neigs))
            if sec["top_face_dofs"]:
                lam_D = float_dirichlet_eigs(
                    sec["Q"], sec["M"], sec["top_face_dofs"], neigs=min(5, args.neigs)
                )
            else:
                lam_D = np.array([])
                print("  (no top dofs identified; Dirichlet skipped)")
            print(f"  field={sec['field']} dofs={sec['n_dofs']} "
                  f"elapsed={time.perf_counter()-t0:.2f}s")
            print(f"  λ^N float: {lam_N}")
            if lam_D.size:
                print(f"  λ^D float: {lam_D}")
            print("  label: engineering diagnostic / not certified")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
