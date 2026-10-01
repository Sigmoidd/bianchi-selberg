"""Group factories and group-bound arithmetic proof backends.

Registering a backend is registering trusted proof code, not accepting
untrusted status strings or treating local witnesses as completeness.
"""
from dataclasses import dataclass
from typing import Callable, TYPE_CHECKING

from fields.quadratic import get_field
from groups.identity import GroupKey
from groups.matrix import MatrixOps, IDENTITY
from groups.arithmetic import QuadraticRing

if TYPE_CHECKING:
    from groups.data import GroupData


@dataclass(frozen=True)
class InventoryBackend:
    key: GroupKey
    proof_id: str
    documentation: tuple[str, ...]
    replay: Callable[["GroupData"], dict]
    dependencies: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    def verify(self, group):
        if group.key != self.key or group.inventory_proof_id != self.proof_id:
            raise ValueError("inventory backend is bound to a different group or proof id")
        if group.inventory_status != "self-contained":
            raise ValueError("inventory does not claim the self-contained proof standard")
        if not self.proof_id or not self.documentation:
            raise ValueError("inventory backend requires a proof id and derivation documentation")
        _verify_element_witnesses(group)
        result = self.replay(group)
        if not isinstance(result, dict) or result.get("proof_id") != self.proof_id:
            raise ValueError("arithmetic replay returned the wrong proof manifest")
        return result


def _verify_element_witnesses(group):
    """Common local identities; the backend separately proves completeness."""
    ops = MatrixOps(QuadraticRing(group.field.d))
    seen = set()
    for C in group.elliptic_classes:
        if C.normalization_status != "proved":
            raise ArithmeticError("class normalization is not proved")
        if not C.label or C.label in seen:
            raise ArithmeticError("class labels must be unique and nonempty")
        seen.add(C.label)
        R, T = C.representative, C.primitive_translation
        if R is None or not group.key.contains(R):
            raise ArithmeticError("elliptic witness is missing or outside the specified group")
        if (C.m not in (2, 3) or ops.canon(R) == ops.canon(IDENTITY)
                or ops.canon(ops.power(R, C.m)) != ops.canon(IDENTITY)):
            raise ArithmeticError("chosen SL lift does not have the claimed elliptic order")
        if C.cuspidal:
            # Its parabolic centralizer and contribution belong to the
            # analytic cusp backend, rather than a loxodromic unit witness.
            continue
        if T is None or not group.key.contains(T):
            raise ArithmeticError("translation witness is missing or outside the specified group")
        if ops.mul(R, T) != ops.mul(T, R):
            raise ArithmeticError("primitive translation does not commute with the elliptic lift")
        if C.flip is not None:
            X = C.flip
            if C.m != 2 or not group.key.contains(X) or ops.power(X, 2) != ops.neg(IDENTITY):
                raise ArithmeticError("invalid endpoint-flip witness")
            if ops.mul(X, R) != ops.neg(ops.mul(R, X)):
                raise ArithmeticError("endpoint flip does not centralize in PSL")


class GroupRegistry:
    def __init__(self):
        self._factories = {}
        self._proofs = {}

    def register(self, key, factory, inventory=None, *, replace=False):
        if key in self._factories and not replace:
            raise ValueError("group already registered")
        if inventory is not None and inventory.key != key:
            raise ValueError("inventory backend key differs from its group factory")
        if replace:
            self._proofs = {k: v for k, v in self._proofs.items() if k[0] != key}
        self._factories[key] = factory
        if inventory is not None:
            self._proofs[key, inventory.proof_id] = inventory

    def get(self, kind):
        from groups.data import GroupData
        if isinstance(kind, GroupData):
            return kind
        key = kind if isinstance(kind, GroupKey) else GroupKey(get_field(kind).d)
        if key not in self._factories:
            raise ValueError(f"no group data registered for {key}; supply its group factory and proof backend")
        group = self._factories[key]()
        if group.key != key:
            raise ValueError("group factory returned a different group identity")
        return group

    def inventory_backend(self, group):
        try:
            return self._proofs[group.key, group.inventory_proof_id]
        except KeyError:
            raise ValueError("no self-contained inventory verifier is registered for this group and proof id") from None

    def verify(self, group):
        return self.inventory_backend(group).verify(group)

    def keys(self):
        return tuple(self._factories)


_DEFAULT = None


def default_registry():
    global _DEFAULT
    if _DEFAULT is None:
        from groups.builtins import register_builtins
        registry = GroupRegistry()
        register_builtins(registry)
        _DEFAULT = registry
    return _DEFAULT
