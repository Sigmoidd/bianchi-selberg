"""Report serialization with a separate, explicit theorem-proof gate."""
import platform
import flint


def report_payload(evaluation):
    E = evaluation
    return dict(
        schema_version=1,
        report_kind="legacy-regression" if E.include_elliptic else "mechanical-screen",
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
    if not E.bound.upper() < 1:
        raise ValueError("this test function does not prove B < 1")
    payload = report_payload(E)
    payload.update(report_kind="spectral-certificate", spectral_certificate=True)
    return payload
