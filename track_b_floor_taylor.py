#!/usr/bin/env python3
"""Validated multivariate Taylor models for the Track-B floor certificate.

Polynomials use normalized variables xi in [-1,1]^3 and a configurable
total-degree truncation.  A model is P(xi)+E with |E| <= error.  Polynomial
coefficients remain complex Arb balls and are summed before any absolute
value is applied.  Products accumulate every discarded coefficient by
multi-index before bounding it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache
import math
from typing import Any

from flint import acb, arb


Index = tuple[int, int, int]
Degree = int | tuple[int, int, int]


def _series_order(degree: Degree) -> int:
    """Univariate composition order needed to fill the polynomial space."""
    return degree if isinstance(degree, int) else sum(degree)


def _kept_index(index: Index, degree: Degree) -> bool:
    if isinstance(degree, int):
        return sum(index) <= degree
    return all(index[axis] <= degree[axis] for axis in range(3))


@dataclass
class TaylorArithmeticAudit:
    """Diagnostic attribution; it does not alter any enclosing bound."""

    discarded_by_multi_index: dict[Index, arb] = field(default_factory=dict)
    maximum_finite_cancellation_condition: arb = field(default_factory=lambda: arb(1))
    maximum_cancellation_multi_index: Index | None = None
    unresolved_cancellation_count: int = 0

    def discarded(self, terms: dict[Index, acb]) -> None:
        for key, value in terms.items():
            self.discarded_by_multi_index[key] = (
                self.discarded_by_multi_index.get(key, arb(0)) + _mag(value)
            ).upper()

    def addition(self, key: Index, left: acb, right: acb, result: acb) -> None:
        numerator = (_mag(left) + _mag(right)).upper()
        denominator = abs(result)
        if bool(denominator.lower() > 0):
            condition = (numerator / denominator.lower()).upper()
            if bool(condition > self.maximum_finite_cancellation_condition):
                self.maximum_finite_cancellation_condition = condition
                self.maximum_cancellation_multi_index = key
        elif bool(numerator > 0):
            self.unresolved_cancellation_count += 1

    def summary(self, top: int = 12) -> dict[str, Any]:
        ranked = sorted(
            self.discarded_by_multi_index.items(),
            key=lambda item: float(item[1]), reverse=True,
        )[:top]
        return {
            "dominant_omitted_multi_indices": [
                {"multi_index": list(key), "aggregate_discarded_upper": str(value)}
                for key, value in ranked
            ],
            "maximum_finite_cancellation_condition": str(
                self.maximum_finite_cancellation_condition.upper()
            ),
            "maximum_cancellation_multi_index": (
                None if self.maximum_cancellation_multi_index is None
                else list(self.maximum_cancellation_multi_index)
            ),
            "unresolved_cancellation_count": self.unresolved_cancellation_count,
            "scope": (
                "all Taylor polynomial additions and discarded product coefficients; "
                "diagnostic only, certified remainder is unchanged"
            ),
        }


_ACTIVE_ARITHMETIC_AUDIT: TaylorArithmeticAudit | None = None


def begin_arithmetic_audit() -> TaylorArithmeticAudit:
    global _ACTIVE_ARITHMETIC_AUDIT
    if _ACTIVE_ARITHMETIC_AUDIT is not None:
        raise RuntimeError("nested Taylor arithmetic audit")
    _ACTIVE_ARITHMETIC_AUDIT = TaylorArithmeticAudit()
    return _ACTIVE_ARITHMETIC_AUDIT


def end_arithmetic_audit() -> TaylorArithmeticAudit:
    global _ACTIVE_ARITHMETIC_AUDIT
    if _ACTIVE_ARITHMETIC_AUDIT is None:
        raise RuntimeError("no active Taylor arithmetic audit")
    result = _ACTIVE_ARITHMETIC_AUDIT
    _ACTIVE_ARITHMETIC_AUDIT = None
    return result


def _mag(q: acb) -> arb:
    return abs(q).upper()


def _point_mid(q: acb) -> acb:
    return acb(q.real.mid(), q.imag.mid())


def _positive_power_upper(x: arb, n: int) -> arb:
    if n == 0:
        return arb(1)
    xu = x.upper()
    if not bool(xu >= 0):
        raise ArithmeticError(f"expected nonnegative power base, got {x}")
    return xu ** n


def _factorial(n: int) -> arb:
    return arb(math.factorial(n))


@dataclass
class TaylorModel:
    degree: Degree
    coeff: dict[Index, acb]
    error: arb

    @classmethod
    def constant(cls, degree: Degree, value: Any) -> "TaylorModel":
        return cls(degree, {(0, 0, 0): acb(value)}, arb(0))

    @classmethod
    def variable(
        cls, degree: Degree, center: arb, halfwidth: arb, index: int
    ) -> "TaylorModel":
        key = [0, 0, 0]
        key[index] = 1
        return cls(
            degree,
            {(0, 0, 0): acb(center), tuple(key): acb(halfwidth)},
            arb(0),
        )

    @staticmethod
    def coerce(other: Any, degree: Degree) -> "TaylorModel":
        if isinstance(other, TaylorModel):
            if other.degree != degree:
                raise ValueError("Taylor-model degree mismatch")
            return other
        return TaylorModel.constant(degree, other)

    def copy(self) -> "TaylorModel":
        return TaylorModel(self.degree, dict(self.coeff), arb(self.error))

    def constant_coefficient(self) -> acb:
        return self.coeff.get((0, 0, 0), acb(0))

    def polynomial_sup(self) -> arb:
        return sum((_mag(q) for q in self.coeff.values()), arb(0)).upper()

    def total_sup(self) -> arb:
        return (self.polynomial_sup() + self.error).upper()

    def evaluate(self, xi: tuple[arb, arb, arb]) -> acb:
        total = acb(0)
        for powers, value in self.coeff.items():
            factor = arb(1)
            for axis in range(3):
                factor *= xi[axis] ** powers[axis]
            total += value * factor
        return acb(
            arb(total.real, total.real.rad() + self.error),
            arb(total.imag, total.imag.rad() + self.error),
        )

    def __add__(self, other: Any) -> "TaylorModel":
        b = self.coerce(other, self.degree)
        out = dict(self.coeff)
        for key, value in b.coeff.items():
            left = out.get(key, acb(0))
            result = left + value
            if _ACTIVE_ARITHMETIC_AUDIT is not None and key in out:
                _ACTIVE_ARITHMETIC_AUDIT.addition(key, left, value, result)
            out[key] = result
        out = {k: v for k, v in out.items() if not (v.real.is_zero() and v.imag.is_zero())}
        return TaylorModel(self.degree, out, (self.error + b.error).upper())

    __radd__ = __add__

    def __neg__(self) -> "TaylorModel":
        return TaylorModel(self.degree, {k: -v for k, v in self.coeff.items()}, self.error)

    def __sub__(self, other: Any) -> "TaylorModel":
        return self + (-self.coerce(other, self.degree))

    def __rsub__(self, other: Any) -> "TaylorModel":
        return self.coerce(other, self.degree) - self

    def __mul__(self, other: Any) -> "TaylorModel":
        b = self.coerce(other, self.degree)
        kept: dict[Index, acb] = {}
        discarded: dict[Index, acb] = {}
        for ka, va in self.coeff.items():
            for kb, vb in b.coeff.items():
                key = (ka[0] + kb[0], ka[1] + kb[1], ka[2] + kb[2])
                target = kept if _kept_index(key, self.degree) else discarded
                target[key] = target.get(key, acb(0)) + va * vb
        tail = sum((_mag(q) for q in discarded.values()), arb(0))
        if _ACTIVE_ARITHMETIC_AUDIT is not None:
            _ACTIVE_ARITHMETIC_AUDIT.discarded(discarded)
        error = (
            tail
            + self.polynomial_sup() * b.error
            + b.polynomial_sup() * self.error
            + self.error * b.error
        ).upper()
        kept = {
            k: v for k, v in kept.items()
            if not (v.real.is_zero() and v.imag.is_zero())
        }
        return TaylorModel(self.degree, kept, error)

    __rmul__ = __mul__

    def __pow__(self, n: int) -> "TaylorModel":
        if n < 0:
            return (self ** (-n)).reciprocal()
        out = TaylorModel.constant(self.degree, 1)
        base = self
        k = n
        while k:
            if k & 1:
                out = out * base
            base = base * base
            k >>= 1
        return out

    def __truediv__(self, other: Any) -> "TaylorModel":
        return self * self.coerce(other, self.degree).reciprocal()

    def __rtruediv__(self, other: Any) -> "TaylorModel":
        return self.coerce(other, self.degree) / self

    def with_extra_error(self, extra: arb) -> "TaylorModel":
        return TaylorModel(self.degree, dict(self.coeff), (self.error + extra).upper())

    def centered(self) -> tuple[acb, "TaylorModel", arb]:
        center = _point_mid(self.constant_coefficient())
        z = self - center
        return center, z, z.total_sup()

    def real_range(self) -> arb:
        center, _z, radius = self.centered()
        if not center.imag.is_zero():
            raise ArithmeticError("real Taylor model has nonreal center")
        return arb(center.real, radius)

    @staticmethod
    def from_series(
        argument: "TaylorModel", coefficients: list[acb], analytic_tail: arb
    ) -> "TaylorModel":
        center, z, _radius = argument.centered()
        del center
        out = TaylorModel.constant(argument.degree, coefficients[0])
        power = TaylorModel.constant(argument.degree, 1)
        for n in range(1, len(coefficients)):
            power = power * z
            out = out + coefficients[n] * power
        return out.with_extra_error(analytic_tail)

    def exp(self) -> "TaylorModel":
        center, _z, radius = self.centered()
        if not center.imag.is_zero():
            raise ArithmeticError("complex exponential model is not implemented")
        c = center.real
        order = _series_order(self.degree)
        coefficients = [acb(c.exp() / _factorial(n)) for n in range(order + 1)]
        upper = self.real_range().upper()
        tail = (
            upper.exp()
            * _positive_power_upper(radius, order + 1)
            / _factorial(order + 1)
        ).upper()
        return self.from_series(self, coefficients, tail)

    def reciprocal(self) -> "TaylorModel":
        center, _z, radius = self.centered()
        if not center.imag.is_zero():
            raise ArithmeticError("complex reciprocal model is not implemented")
        c = center.real
        bounds = self.real_range()
        if not bool(bounds > 0):
            raise ArithmeticError(f"reciprocal model does not prove positivity: {bounds}")
        order = _series_order(self.degree)
        coefficients = [acb(((-1) ** n) / (c ** (n + 1))) for n in range(order + 1)]
        tail = (
            _positive_power_upper(radius, order + 1)
            / (bounds.lower() ** (order + 2))
        ).upper()
        return self.from_series(self, coefficients, tail)

    def sqrt(self) -> "TaylorModel":
        center, _z, radius = self.centered()
        if not center.imag.is_zero():
            raise ArithmeticError("complex square-root model is not implemented")
        c = center.real
        bounds = self.real_range()
        if not bool(bounds > 0):
            raise ArithmeticError(f"square-root model does not prove positivity: {bounds}")
        coefficients: list[acb] = []
        binomial = arb(1)
        order = _series_order(self.degree)
        for n in range(order + 1):
            if n:
                binomial *= (arb(1) / 2 - (n - 1)) / n
            coefficients.append(acb(binomial * c ** (arb(1) / 2 - n)))
        n = order + 1
        binomial *= (arb(1) / 2 - (n - 1)) / n
        tail = (
            abs(binomial).upper()
            * bounds.lower() ** (arb(1) / 2 - n)
            * _positive_power_upper(radius, n)
        ).upper()
        return self.from_series(self, coefficients, tail)

    def sin_cos(self) -> tuple["TaylorModel", "TaylorModel"]:
        center, _z, radius = self.centered()
        if not center.imag.is_zero():
            raise ArithmeticError("complex trigonometric model is not implemented")
        c = center.real
        sin_coeff: list[acb] = []
        cos_coeff: list[acb] = []
        order = _series_order(self.degree)
        for n in range(order + 1):
            phase = n % 4
            sd = (c.sin(), c.cos(), -c.sin(), -c.cos())[phase]
            cd = (c.cos(), -c.sin(), -c.cos(), c.sin())[phase]
            sin_coeff.append(acb(sd / _factorial(n)))
            cos_coeff.append(acb(cd / _factorial(n)))
        tail = (
            _positive_power_upper(radius, order + 1)
            / _factorial(order + 1)
        ).upper()
        return (
            self.from_series(self, sin_coeff, tail),
            self.from_series(self, cos_coeff, tail),
        )

    def sin(self) -> "TaylorModel":
        return self.sin_cos()[0]

    def cos(self) -> "TaylorModel":
        return self.sin_cos()[1]


@dataclass
class TaylorJet:
    value: TaylorModel
    gradient: list[TaylorModel]

    @classmethod
    def constant(cls, degree: Degree, value: Any) -> "TaylorJet":
        zero = TaylorModel.constant(degree, 0)
        return cls(TaylorModel.constant(degree, value), [zero, zero, zero])

    @classmethod
    def variable(
        cls, degree: Degree, center: arb, halfwidth: arb, index: int
    ) -> "TaylorJet":
        gradient = [TaylorModel.constant(degree, 0) for _ in range(3)]
        gradient[index] = TaylorModel.constant(degree, 1)
        return cls(TaylorModel.variable(degree, center, halfwidth, index), gradient)

    @staticmethod
    def coerce(other: Any, degree: Degree) -> "TaylorJet":
        if isinstance(other, TaylorJet):
            if other.value.degree != degree:
                raise ValueError("Taylor-jet degree mismatch")
            return other
        return TaylorJet.constant(degree, other)

    def __add__(self, other: Any) -> "TaylorJet":
        b = self.coerce(other, self.value.degree)
        return TaylorJet(
            self.value + b.value,
            [self.gradient[i] + b.gradient[i] for i in range(3)],
        )

    __radd__ = __add__

    def __neg__(self) -> "TaylorJet":
        return TaylorJet(-self.value, [-q for q in self.gradient])

    def __sub__(self, other: Any) -> "TaylorJet":
        return self + (-self.coerce(other, self.value.degree))

    def __rsub__(self, other: Any) -> "TaylorJet":
        return self.coerce(other, self.value.degree) - self

    def __mul__(self, other: Any) -> "TaylorJet":
        b = self.coerce(other, self.value.degree)
        return TaylorJet(
            self.value * b.value,
            [self.gradient[i] * b.value + self.value * b.gradient[i] for i in range(3)],
        )

    __rmul__ = __mul__

    def __pow__(self, n: int) -> "TaylorJet":
        if n < 0:
            return (self ** (-n)).reciprocal()
        out = TaylorJet.constant(self.value.degree, 1)
        base = self
        k = n
        while k:
            if k & 1:
                out = out * base
            base = base * base
            k >>= 1
        return out

    def reciprocal(self) -> "TaylorJet":
        inverse = self.value.reciprocal()
        return TaylorJet(
            inverse,
            [-self.gradient[i] * inverse * inverse for i in range(3)],
        )

    def __truediv__(self, other: Any) -> "TaylorJet":
        return self * self.coerce(other, self.value.degree).reciprocal()

    def __rtruediv__(self, other: Any) -> "TaylorJet":
        return self.coerce(other, self.value.degree) / self

    def exp(self) -> "TaylorJet":
        value = self.value.exp()
        return TaylorJet(value, [value * q for q in self.gradient])

    def sqrt(self) -> "TaylorJet":
        value = self.value.sqrt()
        return TaylorJet(value, [q / (2 * value) for q in self.gradient])

    def sin_cos(self) -> tuple["TaylorJet", "TaylorJet"]:
        sine, cosine = self.value.sin_cos()
        return (
            TaylorJet(sine, [cosine * q for q in self.gradient]),
            TaylorJet(cosine, [-sine * q for q in self.gradient]),
        )

    def sin(self) -> "TaylorJet":
        return self.sin_cos()[0]

    def cos(self) -> "TaylorJet":
        return self.sin_cos()[1]


@dataclass
class BesselAudit:
    direct_count: int = 0
    fallback_count: int = 0
    majorant_count: int = 0
    real_k_cache: dict[tuple[str, int], arb] = field(default_factory=dict)
    real_order_counts: dict[int, int] = field(default_factory=dict)


def _radial_center_coefficients(ev: Any, mode_index: int, magnitude: arb,
                                center: arb, count: int,
                                audit: BesselAudit) -> list[acb]:
    omega = 2 * ev.pi * magnitude
    arg = omega * center
    kval = acb(arg).bessel_k(ev.order)
    km = acb(arg).bessel_k(ev.order - 1)
    kp = acb(arg).bessel_k(ev.order + 1)
    audit.direct_count += 3
    if not (kval.is_finite() and km.is_finite() and kp.is_finite()):
        audit.fallback_count += 1
        raise ArithmeticError("non-finite direct Arb Bessel value in Taylor center")
    radial = center * kval
    radial_y = kval - center * omega * (km + kp) / 2
    coefficients = [radial, radial_y]
    lam = 1 + ev.r * ev.r
    c2 = center * center
    omega2 = omega * omega
    for n in range(count - 2):
        am1 = coefficients[n - 1] if n >= 1 else acb(0)
        am2 = coefficients[n - 2] if n >= 2 else acb(0)
        numerator = (
            (n + 1) * center * (2 * n - 1) * coefficients[n + 1]
            + (n * n - 2 * n + lam - omega2 * c2) * coefficients[n]
            - 2 * omega2 * center * am1
            - omega2 * am2
        )
        coefficients.append(-numerator / (c2 * (n + 2) * (n + 1)))
    return coefficients


def _k_derivative_majorant(order: int, argument_lower: arb,
                            audit: BesselAudit) -> arb:
    total = arb(0)
    for j in range(order + 1):
        shifted = abs(-order + 2 * j)
        key = (str(argument_lower), shifted)
        bound = audit.real_k_cache.get(key)
        if bound is None:
            value = acb(argument_lower).bessel_k(shifted)
            audit.majorant_count += 1
            audit.real_order_counts[shifted] = audit.real_order_counts.get(shifted, 0) + 1
            if not value.is_finite():
                audit.fallback_count += 1
                raise ArithmeticError("non-finite real-order Bessel tail majorant")
            bound = abs(value).upper()
            audit.real_k_cache[key] = bound
        total += math.comb(order, j) * bound
    return (total / (2 ** order)).upper()


def radial_derivative_majorant(ev: Any, magnitude: arb, y_bounds: arb,
                               order: int, audit: BesselAudit) -> arb:
    omega = 2 * ev.pi * magnitude
    argument_lower = omega.lower() * y_bounds.lower()
    kd = _k_derivative_majorant(order, argument_lower, audit)
    first = y_bounds.upper() * _positive_power_upper(omega, order) * kd
    if order == 0:
        return first.upper()
    kdm1 = _k_derivative_majorant(order - 1, argument_lower, audit)
    second = order * _positive_power_upper(omega, order - 1) * kdm1
    return (first + second).upper()


def radial_taylor_jet(ev: Any, mode_index: int, magnitude: arb,
                      y: TaylorJet, audit: BesselAudit) -> TaylorJet:
    degree = y.value.degree
    order = _series_order(degree)
    center_acb, _z, radius = y.value.centered()
    if not center_acb.imag.is_zero():
        raise ArithmeticError("radial Taylor center is nonreal")
    center = center_acb.real
    y_bounds = y.value.real_range()
    if not bool(y_bounds > 0):
        raise ArithmeticError(f"radial Taylor height not positive: {y_bounds}")
    # a_n are Taylor coefficients r^(n)(center)/n! from the radial ODE.
    a = _radial_center_coefficients(
        ev, mode_index, magnitude, center, order + 2, audit
    )
    value_tail = (
        radial_derivative_majorant(ev, magnitude, y_bounds, order + 1, audit)
        * _positive_power_upper(radius, order + 1)
        / _factorial(order + 1)
    ).upper()
    derivative_tail = (
        radial_derivative_majorant(ev, magnitude, y_bounds, order + 2, audit)
        * _positive_power_upper(radius, order + 1)
        / _factorial(order + 1)
    ).upper()
    value = TaylorModel.from_series(y.value, a[:order + 1], value_tail)
    derivative_coefficients = [acb((n + 1) * a[n + 1]) for n in range(order + 1)]
    derivative = TaylorModel.from_series(
        y.value, derivative_coefficients, derivative_tail
    )
    return TaylorJet(value, [derivative * q for q in y.gradient])


def _double_factorial(n: int) -> int:
    out = 1
    for k in range(n, 0, -2):
        out *= k
    return out


@lru_cache(maxsize=None)
def monomial_legendre_expansion(n: int) -> tuple[tuple[int, Fraction], ...]:
    out = []
    for k in range(n // 2 + 1):
        ell = n - 2 * k
        numerator = (2 * ell + 1) * math.factorial(n)
        denominator = (
            (2 ** k) * math.factorial(k) * _double_factorial(2 * n - 2 * k + 1)
        )
        out.append((ell, Fraction(numerator, denominator)))
    return tuple(out)


def legendre_coefficients(model: TaylorModel) -> dict[Index, acb]:
    out: dict[Index, acb] = {}
    for powers, value in model.coeff.items():
        for a, ca in monomial_legendre_expansion(powers[0]):
            for b, cb in monomial_legendre_expansion(powers[1]):
                for c, cc in monomial_legendre_expansion(powers[2]):
                    key = (a, b, c)
                    scale = ca * cb * cc
                    factor = arb(scale.numerator) / scale.denominator
                    out[key] = out.get(key, acb(0)) + factor * value
    return out


def polynomial_l2_squared_upper(model: TaylorModel, h1: arb, h2: arb,
                                hs: arb) -> tuple[arb, int]:
    legendre = legendre_coefficients(model)
    normalized = arb(0)
    for (a, b, c), value in legendre.items():
        weight = arb(8) / ((2 * a + 1) * (2 * b + 1) * (2 * c + 1))
        normalized += weight * (_mag(value) ** 2)
    return (h1 * h2 * hs * normalized).upper(), len(legendre)
