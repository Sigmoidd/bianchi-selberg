#!/usr/bin/env python3
"""Validated real K_{i r}(x) values for the dual-certification pipeline.

The production Hejhal matrix must never use an analytic upper bound as a
surrogate value.  This module therefore fails closed when python-flint is not
available or when Arb does not return a finite enclosure containing zero in
the (mathematically vanishing) imaginary component.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, dataclass
from functools import lru_cache
import math
from typing import Iterator

try:
    from flint import acb, arb, ctx
except Exception as exc:  # pragma: no cover - fail-closed import path
    raise RuntimeError("python-flint is required for validated K_{ir}") from exc


@dataclass(frozen=True)
class KirEnclosure:
    r: float
    x: float
    midpoint: float
    radius: float
    lower: float
    upper: float
    precision_bits: int

    def to_dict(self) -> dict:
        return asdict(self)


@contextmanager
def _precision(bits: int) -> Iterator[None]:
    old = int(ctx.prec)
    ctx.prec = max(int(bits), 64)
    try:
        yield
    finally:
        ctx.prec = old


@lru_cache(maxsize=131072)
def _kir_cached(r_text: str, x_text: str, bits: int) -> KirEnclosure:
    with _precision(bits):
        rr = arb(r_text)
        xx = arb(x_text)
        if not rr.is_finite() or not xx.is_finite() or not bool(xx > 0):
            raise ValueError(f"K_ir requires finite r and x>0; got r={rr}, x={xx}")
        z = acb(xx).bessel_k(acb(0, rr))
        if not z.is_finite():
            raise ArithmeticError(f"Arb failed to enclose K_{{i r}}(x): {z}")
        # K_{i r}(x) is real for real r and x>0.  Arb should independently
        # enclose zero in the imaginary part; otherwise reject the result.
        if not z.imag.contains(0):
            raise ArithmeticError(f"unexpected non-real K_{{i r}}(x) enclosure: {z}")
        mid = float(z.real.mid())
        real_rad = float(z.real.rad())
        imag_rad = abs(float(z.imag.mid())) + float(z.imag.rad())
        rad = max(real_rad, imag_rad, 0.0)
        if not (math.isfinite(mid) and math.isfinite(rad)):
            raise ArithmeticError(f"nonfinite K_{{i r}}(x) midpoint/radius: {z}")
        return KirEnclosure(
            r=float(r_text), x=float(x_text), midpoint=mid, radius=rad,
            lower=mid - rad, upper=mid + rad, precision_bits=int(bits),
        )


def kir_enclosure(r: float, x: float, bits: int = 160) -> KirEnclosure:
    """Return an Arb-backed enclosure using round-trippable decimal inputs."""
    return _kir_cached(repr(float(r)), repr(float(x)), int(bits))


def kir_mid_rad(r: float, x: float, bits: int = 160) -> tuple[float, float]:
    q = kir_enclosure(r, x, bits)
    return q.midpoint, q.radius


def kir_abs_upper(r: float, x: float, bits: int = 160) -> float:
    q = kir_enclosure(r, x, bits)
    return max(abs(q.lower), abs(q.upper))


def clear_cache() -> None:
    _kir_cached.cache_clear()
