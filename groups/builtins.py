"""Adapters for existing results and incomplete full-group screens."""
from fields.quadratic import DISCRIMINANTS, QuadraticField
from groups.data import EllipticClass, GroupData, PICARD, EISENSTEIN
from groups.identity import GroupKey
from groups.registry import InventoryBackend


def d2_group():
    from groups.d2_inventory import expected_classes
    provenance = "docs/D2_INVENTORY_PROOF.md; groups/d2_inventory.py (d2-arithmetic-v1)"
    classes = tuple(EllipticClass(**C, cuspidal=False, provenance=provenance,
                                   normalization_status="proved") for C in expected_classes())
    return GroupData(QuadraticField(2), "PSL2(Z[sqrt(-2)])", (0, 1), classes,
                     "self-contained", provenance, inventory_proof_id="d2-arithmetic-v1")


def d2_backend():
    from groups.d2_inventory import verify_group_records
    return InventoryBackend(
        GroupKey(2), "d2-arithmetic-v1", ("docs/D2_INVENTORY_PROOF.md",), verify_group_records,
        dependencies=("Minkowski ideal-class bound and elementary local valuation theory",),
        notes=("The two inverse order-3 element classes are counted separately.",))


def d7_group():
    from groups.d7_inventory import expected_classes
    provenance = "docs/D7_INVENTORY_PROOF.md; groups/d7_inventory.py (d7-arithmetic-v1)"
    classes = tuple(EllipticClass(**C, cuspidal=False, provenance=provenance,
                                   normalization_status="proved") for C in expected_classes())
    return GroupData(QuadraticField(7), "PSL2(O_-7)", (0, 1), classes,
                     "self-contained", provenance, inventory_proof_id="d7-arithmetic-v1")


def d7_backend():
    from groups.d7_inventory import verify_group_records
    return InventoryBackend(
        GroupKey(7), "d7-arithmetic-v1", ("docs/D7_INVENTORY_PROOF.md",), verify_group_records,
        dependencies=("Minkowski ideal-class bound and elementary local valuation theory",),
        notes=("Inverse order-3 elements are conjugate; involution centralizer has no endpoint flip.",))


def d11_group():
    from groups.d11_inventory import expected_classes
    provenance = "docs/D11_INVENTORY_PROOF.md; groups/d11_inventory.py (d11-arithmetic-v1)"
    classes = tuple(EllipticClass(**C, cuspidal=False, provenance=provenance,
                                  normalization_status="proved") for C in expected_classes())
    return GroupData(QuadraticField(11), "PSL2(O_-11)", (0, 1), classes,
                     "self-contained", provenance, inventory_proof_id="d11-arithmetic-v1")


def d11_backend():
    from groups.d11_inventory import verify_group_records
    return InventoryBackend(
        GroupKey(11), "d11-arithmetic-v1", ("docs/D11_INVENTORY_PROOF.md",), verify_group_records,
        dependencies=("Minkowski ideal-class bound and elementary local valuation theory",),
        notes=("The two inverse order-3 classes are distinct; the involution has an endpoint flip.",))


def d19_group():
    from groups.d19_inventory import expected_classes
    provenance = "docs/D19_INVENTORY_PROOF.md; groups/d19_inventory.py (d19-arithmetic-v1)"
    classes = tuple(EllipticClass(**C, cuspidal=False, provenance=provenance,
                                  normalization_status="proved") for C in expected_classes())
    return GroupData(QuadraticField(19), "PSL2(O_-19)", (0, 1), classes,
                     "self-contained", provenance, inventory_proof_id="d19-arithmetic-v1")


def d19_backend():
    from groups.d19_inventory import verify_group_records
    return InventoryBackend(
        GroupKey(19), "d19-arithmetic-v1", ("docs/D19_INVENTORY_PROOF.md",), verify_group_records,
        dependencies=("Minkowski ideal-class bound and elementary local valuation theory",),
        notes=("Inverse order-3 classes merge; the involution has an endpoint flip.",))


def incomplete_group(d):
    return GroupData(QuadraticField(d), f"PSL2(O_-{d})", (0, 1) if d <= 19 else (3, 0),
                     (), "incomplete", "docs/INVENTORY_PROOF.md (open completeness obligations)")


def register_builtins(registry):
    registry.register(GroupKey(1), lambda: PICARD)
    registry.register(GroupKey(3), lambda: EISENSTEIN)
    registry.register(GroupKey(2), d2_group, d2_backend())
    registry.register(GroupKey(7), d7_group, d7_backend())
    registry.register(GroupKey(11), d11_group, d11_backend())
    registry.register(GroupKey(19), d19_group, d19_backend())
    for d in DISCRIMINANTS:
        if d not in (1, 2, 3, 7, 11, 19):
            registry.register(GroupKey(d), lambda d=d: incomplete_group(d))
