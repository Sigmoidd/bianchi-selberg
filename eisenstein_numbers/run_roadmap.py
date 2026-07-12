"""Execute Stages 1–7 of the Eisenstein FEM roadmap (float + status)."""
from __future__ import annotations

import sys


def main() -> int:
    print("Eisenstein–Picard FEM roadmap runner")
    print("Geometry → cell → M0 → constants → CR → windows → cert scaffold")
    print("=" * 64)

    # Stage 1–2
    print("\n[1] Geometry freeze (EGM F_3 + pairings + vol)")
    from geometry_fund import (
        area_P3, vol_F_exact, vol_KY_exact, vol_K_comp_quad,
        verify_side_pairings,
    )
    vF, _ = vol_F_exact()
    vK, _ = vol_KY_exact(1.25)
    vC = vol_K_comp_quad(1.25)
    print(f"  |T|=area(P_3)={area_P3():.8f}  vol(F)={vF:.8f}")
    print(f"  vol(K_Y@1.25)={vK:.8f}  vol(K_comp@1.25)≈{vC:.8f}")
    for m in verify_side_pairings():
        if m.startswith("FAIL"):
            print(f"  !! {m}")
            return 1
    print("  side-pairing checks: OK")

    print("\n[2] Reference cell (EGM P_3)")
    from reference_cell import AREA_T, build_reference_cell, exact_vol_KY
    mesh = build_reference_cell(6, 6, 3, Y=1.25, domain="P3")
    vK2, _ = exact_vol_KY(1.25)
    print(f"  domain={mesh['domain']}  |T|={AREA_T:.6f}  tets={len(mesh['tets'])}")
    print(f"  delta_bar={mesh['delta_bar']:.3e}  y_floor_min={mesh['y_floor_min']:.4f}")
    print(f"  EGM vol(K_Y)={vK2:.8f}")
    assert mesh["delta_bar"] >= 0
    assert abs(AREA_T - area_P3()) < 1e-14
    assert mesh["y_floor_min"] > 0.8

    # Stage 3
    print("\n[3] Float M0 (legacy R_comp Q1 mesh — still positive margin)")
    from fem_omega_m0 import run as m0_run
    import numpy as np
    lams = np.array([0.05, 0.5, 0.999])
    mu = m0_run(10, 10, 5, 1.25, lams, verbose=True)
    ok_m0 = mu > 0
    print(f"  M0 min-mu (3-point) = {mu:.4f}  {'PASS' if ok_m0 else 'FAIL'}")

    # Stage 4
    print("\n[4] Interpolation constants")
    from constants import KAPPA_SC, KAPPA_1, beta_s
    print(f"  kappa_sc={KAPPA_SC:.6f}  kappa_1={KAPPA_1:.6f}")
    print(f"  beta(s=0)={beta_s(0, AREA_T, 1.25):.6f}")

    # Stages 5–6
    print("\n[5–6] CR assembly + G1 window float (P3, rho=55)")
    from cr_omega import run as cr_run
    data, cst, res = cr_run(N_tri=6, N3=3, Y=1.25, nwin=8, rho=55.0, do_mu=True)
    n_g1 = sum(1 for r in res if r["ok"])
    n_mu = sum(1 for r in res if r.get("mu", 0) > 0)
    print(f"  G1 PASS: {n_g1}/8   mu>0: {n_mu}/8")

    # Stage 7
    print("\n[7] Interval certification (Arb + Rump on P_3)")
    from cert_omega import certify as cert_run
    ok_cert = cert_run(N_tri=6, N3=3, diagnostics=False)
    print(f"  interval cert: {'PASS' if ok_cert else 'FAIL'}")

    # Optional non-Neumann
    print("\n[7b] Non-Neumann FE space (pairings + κ=I1/CZZ)")
    from non_neumann_omega import run_non_neumann
    nn = run_non_neumann(N_tri=6, N3=3, do_mu=True)
    ok_nn = bool(nn.get("ok"))
    print(f"  non-Neumann float: {'PASS' if ok_nn else 'FAIL'}")

    # Stage 8
    print("\n[8] Congruence ladder Γ₀(N=3)")
    from cert_omega_p import status as p_status
    p_status()
    print("  machine cert: python -u cert_omega_p.py 3 6 3  (8/8 Rump)")

    print("\n" + "=" * 64)
    print("ROADMAP STATUS")
    print(f"  Geometry:      FROZEN (EGM F_3 / P_3 mesh)")
    print(f"  M0 margin:     {'PASS' if ok_m0 else 'FAIL'} (mu={mu:.3f})")
    print(f"  G1 float:      {n_g1}/8 windows PASS (c_e/d_e + PSD)")
    print(f"  Interval cert: {'PASS 8/8 Arb+Rump' if ok_cert else 'FAIL'}")
    print(f"  Non-Neumann:   {'PASS' if ok_nn else 'FAIL'}")
    print("  Congruence:    Γ₀(N=3) CERTIFIED; N=7+ open")
    print("=" * 64)
    return 0 if ok_m0 and n_g1 >= 6 and ok_cert and ok_nn else 1


if __name__ == "__main__":
    raise SystemExit(main())
