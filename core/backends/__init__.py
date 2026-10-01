"""Analytic backend contract for geometry, cusps, and scattering.

A backend supplies verified volume/systole and the aggregate cusp/scattering
terms. Test functions, support checks, elliptic sums and certificate gates
remain common across fields and levels.
"""
from dataclasses import dataclass
from typing import Protocol
from flint import arb
from groups.identity import GroupKey


@dataclass(frozen=True)
class Geometry:
    group_key: GroupKey
    volume: arb
    systole: arb


class AnalyticBackend(Protocol):
    backend_id: str
    documentation: tuple[str, ...]

    def validate(self, group) -> None: ...
    def geometry(self, group) -> Geometry: ...
    def terms(self, group, geometry, k, delta, R, include_elliptic) -> dict[str, arb]: ...


class AnalyticRegistry:
    def __init__(self):
        self._backends = {}

    def register(self, backend: AnalyticBackend):
        if not backend.backend_id or not backend.documentation:
            raise ValueError("analytic backend requires an id and derivation documentation")
        if backend.backend_id in self._backends:
            raise ValueError("analytic backend already registered")
        self._backends[backend.backend_id] = backend

    def get(self, backend_id):
        try:
            return self._backends[backend_id]
        except KeyError:
            raise ValueError(f"no analytic backend registered for {backend_id}") from None


def validate_geometry(group, geometry):
    if geometry.group_key != group.key:
        raise ValueError("analytic geometry is bound to a different group identity")
    if (not isinstance(geometry.volume, arb) or not geometry.volume.is_finite()
            or not geometry.volume > 0):
        raise ValueError("analytic backend must prove positive volume")
    if (not isinstance(geometry.systole, arb) or not geometry.systole.is_finite()
            or not geometry.systole > 0):
        raise ValueError("analytic backend must prove a positive systole")


_DEFAULT = None


def default_registry():
    global _DEFAULT
    if _DEFAULT is None:
        from core.backends.level_one import LevelOneBackend
        registry = AnalyticRegistry()
        registry.register(LevelOneBackend())
        _DEFAULT = registry
    return _DEFAULT


def get_backend(group, registry=None):
    return (registry or default_registry()).get(group.trace_backend_id)
