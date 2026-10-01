"""Uniform trace assembly; cusp and scattering formulas come from a backend."""
from dataclasses import dataclass
import math
from flint import arb, ctx
from groups import get_group
from core.bspline import g0, gpp0, validate_k
from core.testfunctions import h_isig
from core.backends import get_backend, validate_geometry
# Compatibility exports for the historical bianchi_omega_arb module.
from core.numerics import integ, setp, PREC, RTOL


@dataclass(frozen=True)
class Evaluation:
    group: object
    k: int
    delta: float
    R: float
    include_elliptic: bool
    precision_bits: int
    bound: object
    terms: dict


def evaluate(kind, k=2, frac=0.999, R=40, verbose=True, include_elliptic=True,
             *, group_registry=None, analytic_registry=None):
    validate_k(k)
    if not (arb(frac) > 0 and arb(frac) <= 1):
        raise ValueError("support fraction must lie in (0, 1]")
    if not math.isfinite(float(R)) or R < 7:
        raise ValueError("shifted-line tail bound requires finite R >= 7")
    setp()
    group = get_group(kind, registry=group_registry)
    backend = get_backend(group, registry=analytic_registry)
    backend.validate(group)
    if include_elliptic:
        group.require_inventory(registry=group_registry)
    geometry = backend.geometry(group)
    validate_geometry(group, geometry)
    # Float delta is an exact binary input. Follow rounding with an Arb check.
    d = math.nextafter(float(geometry.systole.lower())*float(frac)/(2*k), -math.inf)
    if not (arb(d) > 0 and (2*k*arb(d)).upper() <= geometry.systole.lower()):
        raise ValueError("support 2k*delta is not provably within the systole")
    I = geometry.volume/(2*arb.pi())*(-gpp0(k,d))
    coefficient = sum((C.coefficient() for C in group.elliptic_classes if not C.cuspidal), arb(0))
    NCE = g0(k,d)*coefficient if include_elliptic else arb(0)
    extra = backend.terms(group, geometry, k, d, R, include_elliptic)
    required = {"CE", "Ch0", "PARg0", "PSI", "PHIINT"}
    if not required <= extra.keys() or any(name in extra for name in ("I", "NCE", "hi")):
        raise ValueError("analytic backend returned an invalid trace-term block")
    if any(not isinstance(value, arb) or not value.is_finite() for value in extra.values()):
        raise ValueError("analytic backend terms must be finite Arb enclosures")
    hi = h_isig(k,d)
    B = I+NCE+extra["CE"]+extra["Ch0"]+extra["PARg0"]+extra["PSI"]+extra["PHIINT"]-hi
    terms = dict(I=I, NCE=NCE, **extra, hi=hi)
    if verbose:
        for name in ("I", "NCE", "CE", "Ch0", "PARg0", "PSI", "PHIINT", "hi"):
            print(f"    {name:12s} = {terms[name]}")
        print(f"    {'B':12s} = {B}")
    return Evaluation(group, k, d, R, include_elliptic, ctx.prec, B, terms)
