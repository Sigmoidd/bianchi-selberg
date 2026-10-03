from dataclasses import dataclass
from fractions import Fraction

from flint import arb

from fields.quadratic import QuadraticField, get_field


@dataclass(frozen=True)
class EllipticClass:
    label: str
    m: int
    finite_centralizer_order: int
    # N(T0) = a + b sqrt(c), stored exactly.
    norm: tuple[int, int, int]
    cuspidal: bool
    provenance: str
    normalization_status: str = "historical"

    def coefficient(self):
        if self.cuspidal:
            raise ValueError("cuspidal elliptics use the separate cusp formula")
        if self.m not in (2, 3) or self.finite_centralizer_order < 1:
            raise ValueError("invalid order-2/order-3 elliptic data")
        a, b, c = self.norm
        N = arb(a)+b*arb(c).sqrt()
        if not N > 1:
            raise ValueError("primitive loxodromic norm must exceed 1")
        sin2 = arb(1) if self.m == 2 else arb(3)/4
        return N.log()/(4*self.finite_centralizer_order*sin2)


@dataclass(frozen=True)
class GroupData:
    field: QuadraticField
    name: str
    # Coordinates of tau in the integral basis documented in SYSTOLES.md.
    systole_trace: tuple[int, int]
    elliptic_classes: tuple[EllipticClass, ...]
    inventory_status: str
    inventory_provenance: str
    cusp_count: int = 1
    ce_g0: tuple[int, int, int] = (0, 1, 1)  # a/b * log(c)
    ce_integral: Fraction = Fraction(0)
    ce_kernel: Fraction = Fraction(1)

    @property
    def GG(self):
        return self.field.units//2

    def systole(self):
        a, b = self.systole_trace
        d = self.field.d
        x = Fraction(a) if self.field.D % 4 == 0 else Fraction(2*a+b, 2)
        y2 = Fraction(d*b*b) if self.field.D % 4 == 0 else Fraction(d*b*b, 4)
        A = x*x+y2
        rad = (A-4)**2+16*y2
        v = (_arb_fraction(A)+_arb_fraction(rad).sqrt())/4
        return (v+(v*v-1).sqrt()).log()

    def require_inventory(self):
        if self.inventory_status not in ("legacy", "self-contained"):
            raise ValueError(f"{self.name}: self-contained elliptic inventory is incomplete; "
                             "only a mechanical screen is available")

    def analytic_data(self, require_inventory=True):
        from fields.quadratic import volume, eta
        if require_inventory:
            self.require_inventory()
        a, b, c = self.ce_g0
        return dict(vol=volume(self.field), systole=self.systole(),
                    C_ell=sum((e.coefficient() for e in self.elliptic_classes
                               if not e.cuspidal), arb(0)),
                    CEg0=arb(a)/b*arb(c).log(),
                    CEint=_arb_fraction(self.ce_integral),
                    ckern=_arb_fraction(self.ce_kernel), GG=self.GG,
                    eta=eta(self.field), D=-self.field.D)


def _arb_fraction(x):
    x = Fraction(x)
    return arb(x.numerator)/x.denominator


PICARD = GroupData(
    QuadraticField(1), "PSL2(Z[i])", (0, 1),
    (EllipticClass("Picard order 3", 3, 3, (7, 4, 3), False,
                   "verify_matthies.py; old_RIGOR_GAPS.md"),),
    "legacy", "Historical class count; new-field proof standard is stricter",
    ce_g0=(5, 16, 2), ce_integral=Fraction(1, 4))
EISENSTEIN = GroupData(
    QuadraticField(3), "PSL2(Z[omega])", (-1, 1),
    (EllipticClass("Eisenstein order 2 (frozen)", 2, 2, (7, 4, 3), False,
                   "old_RIGOR_GAPS.md; docs/NORMALIZATION_ISSUE.md",
                   "unresolved-flips-frozen"),),
    "legacy", "Historical EGM classification citation; normalization issue remains open",
    ce_g0=(2, 9, 3), ce_integral=Fraction(1, 3), ce_kernel=Fraction(1, 2))


def get_group(kind):
    if isinstance(kind, GroupData):
        return kind
    F = get_field(kind)
    if F.d == 1:
        return PICARD
    if F.d == 3:
        return EISENSTEIN
    witness = (0, 1) if F.d == 2 else ((0, 1) if F.d <= 19 else (3, 0))
    return GroupData(F, f"PSL2(O_-{F.d})", witness, (), "incomplete",
                     "docs/INVENTORY_PROOF.md (open completeness obligations)")
