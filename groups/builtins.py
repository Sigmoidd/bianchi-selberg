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


def incomplete_group(d):
    return GroupData(QuadraticField(d), f"PSL2(O_-{d})", (0, 1) if d <= 19 else (3, 0),
                     (), "incomplete", "docs/INVENTORY_PROOF.md (open completeness obligations)")


def register_builtins(registry):
    registry.register(GroupKey(1), lambda: PICARD)
    registry.register(GroupKey(3), lambda: EISENSTEIN)
    registry.register(GroupKey(2), d2_group, d2_backend())
    for d in DISCRIMINANTS:
        if d not in (1, 2, 3):
            registry.register(GroupKey(d), lambda d=d: incomplete_group(d))
