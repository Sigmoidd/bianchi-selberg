"""Report serialization with a separate, explicit theorem-proof gate."""
import platform
import flint


def report_payload(evaluation):
    E = evaluation
    return dict(
        schema_version=1,
        report_kind=("mechanical-screen" if not E.include_elliptic else
                     "legacy-regression" if E.group.inventory_status == "legacy" else
                     "proved-inventory-evaluation"),
        spectral_certificate=False,
        group=E.group.name,
        d=E.group.field.d,
        D=E.group.field.D,
        parameters=dict(k=E.k, delta_hex=E.delta.hex(), R=E.R, precision_bits=E.precision_bits),
        arithmetic=dict(python=platform.python_version(), python_flint=flint.__version__),
        inventory_status=E.group.inventory_status,
        inventory_provenance=E.group.inventory_provenance,
        normalization_status=[e.normalization_status for e in E.group.elliptic_classes],
        systole_trace=list(E.group.systole_trace),
        systole_ball=E.group.systole().str(40),
        bound_ball=E.bound.str(40),
        lower_endpoint_ball=E.bound.lower().str(40),
        upper_endpoint_ball=E.bound.upper().str(40),
        terms={name: value.str(40) for name, value in E.terms.items()},
        notes=["Ball strings retain their radii; decimal endpoints are not exact scalars.",
               "Historical reports preserve the existing elliptic normalization.",
               "Mechanical screens omit elliptics and do not imply a spectral gap."],
    )


def certificate_payload(evaluation):
    E = evaluation
    if not E.include_elliptic or E.group.inventory_status != "self-contained":
        raise ValueError("a new certificate requires a self-contained complete inventory")
    if any(e.normalization_status != "proved" for e in E.group.elliptic_classes):
        raise ValueError("all elliptic normalizations must be proved")
    proof = E.group.verify_inventory()
    from groups.systoles import verify_systole
    from flint import arb
    verify_systole(E.group)
    if not (arb(E.delta) > 0 and (2*E.k*arb(E.delta)).upper() <= E.group.systole().lower()):
        raise ValueError("certificate support is outside the proved systole")
    if not E.bound.upper() < 1:
        raise ValueError("this test function does not prove B < 1")
    payload = report_payload(E)
    payload.update(report_kind="spectral-certificate", spectral_certificate=True)
    payload.update(
        inventory_proof=proof,
        support_ball=(2*E.k*arb(E.delta)).str(40),
        elliptic_coefficient_ball=E.group.analytic_data()["C_ell"].str(40),
        conclusion="No discrete Laplace eigenvalue in (0,1) at level 1.",
        proof_dependencies=[
            "Friedman, arXiv:math/0612807v1, Theorem 4.1.1 (standard trace formula)",
            "docs/D2_INVENTORY_PROOF.md (arithmetic and orbital derivation)",
            "Minkowski ideal-class bound and elementary local valuation theory",
            "docs/SYSTOLES.md; docs/ANALYTIC_DERIVATIONS.md (analytic reductions)",
            "Arb interval arithmetic and certified quadrature via python-flint",
        ],
        notes=["Ball strings retain their radii; decimal endpoints are not exact scalars.",
               "Completeness is proved in the derivation note; finite arithmetic witnesses replay exactly.",
               "This is a mathematical proof with executable checks, not a proof-assistant formalization.",
               "The two inverse order-3 element classes are counted separately.",
               "Existing Eisenstein coefficients and frozen historical reports are unchanged."],
    )
    return payload
