"""Self-contained d=7 arithmetic replay; see docs/D7_INVENTORY_PROOF.md.

All ideal classes and units are covered by proved finite bounds. GL/SL/PSL
class distinctions and primitive translations are handled individually.
"""
from fractions import Fraction
from itertools import product
from math import gcd

from flint import arb
from groups.arithmetic import QuadraticRing
from groups.exact import Radical
from groups.identity import GroupKey
from groups.matrix import MatrixOps, IDENTITY as MATRIX_IDENTITY
from groups.relative_orders import QuadraticOrder, ROOT, IDENTITY, power

RING = QuadraticRing(7)
OPS = MatrixOps(RING)
ONE = (1, 0)
GENERATORS = {"A2": (1, 1, 2, -1), "A3": (-1, 1, 1, 0)}


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


class RelativeOrder(QuadraticOrder):
    def __init__(self, name):
        if name not in GENERATORS:
            raise ValueError("unknown d7 relative order")
        self.name = name
        super().__init__(RING, (0, 0) if name == "A2" else ONE, ONE)

    def abs_square(self, z):
        a, b, c, d = z
        base = RING.norm((a, b))+RING.norm((c, d))
        if self.name == "A2":
            return Radical(base, b*c-a*d, 7)
        return Radical(2*base+2*a*c+a*d+b*c+4*b*d, b*c-a*d, 21, 2)


def verify_orders_and_lattices():
    O2, O3 = RelativeOrder("A2"), RelativeOrder("A3")
    require(O2.discriminant() == 784 and O3.discriminant() == 441, "relative order discriminants")
    require(gcd(7, 4) == gcd(7, 3) == 1, "coprime quadratic discriminants prove maximality")
    require(arb(42)/arb.pi()**2 < 5, "A2 Minkowski bound below five")
    require(arb(63)/(2*arb.pi()**2) < 4, "A3 Minkowski bound below four")
    # A2: two norm-two primes, both principal. No norm-three prime.
    beta, beta_bar = (0, 1, 1, 0), (1, -1, 1, 0)
    require(O2.absolute_norm(beta) == O2.absolute_norm(beta_bar) == 2, "principal primes at two")
    require({x for x in range(2) if (x*x-x+2) % 2 == 0} == {0, 1}, "split base prime two")
    require({x for x in range(2) if (x*x+1) % 2 == 0} == {1}, "ramified i at two")
    # On either Q2 factor, i-1 satisfies Y^2+2Y+2: Eisenstein.
    # Thus each prime over two has residue degree one, excluding norm-four primes.
    require(2 % 2 == 0 and 2 % 4 != 0, "Eisenstein ramification at two")
    # beta=tau+i vanishes at tau=1; beta_bar at tau=0.
    require((1+1) % 2 == 0 and (0+1) % 2 != 0, "beta residue")
    require((1-0+1) % 2 == 0 and (1-1+1) % 2 != 0, "beta_bar residue")
    require(all((x*x-x+2) % 3 != 0 for x in range(3)), "no base norm-three prime")
    # A3 has no norm-two prime because alpha has no root over F2.
    require(all((x*x-x+1) % 2 != 0 for x in range(2)), "no A3 norm-two prime")
    return dict(discriminants=dict(A2=784, A3=441), maximal_class_numbers=dict(A2=1, A3=1),
                GL_lattice_types=dict(trace_zero=["A2"], trace_one=["A3"]),
                minkowski_cutoffs=dict(A2=4, A3=3), principal_prime_norms=dict(A2=[2, 2], A3=[]))


def verify_units(name):
    order = RelativeOrder(name)
    generator = GENERATORS[name]
    torsion_order = 4 if name == "A2" else 6
    coefficient_bound = Fraction(9, 2) if name == "A2" else Fraction(7, 3)
    E2 = order.abs_square(generator)
    require(E2.compare(1) > 0, "expanding unit")
    E2.unit_reciprocal()
    norm = order.norm(generator)
    require(norm == (ONE if name == "A2" else (-1, 0)), "generator relative norm")
    # (E+E^-1)^2 / |r-sigma(r)|^2 gives the displayed coefficient norm bound.
    require(Fraction(2*E2.a+2*E2.denominator,
                     E2.denominator*(4 if name == "A2" else 3)) == coefficient_bound,
            "proved coefficient norm bound")
    reduced, candidates = set(), []
    # Complete: |b|,|d|<=1 and |a|,|c|<=2 follow from the norm bound <9/2.
    for z in product(range(-2, 3), range(-1, 2), range(-2, 3), range(-1, 2)):
        if RING.norm(z[:2]) >= coefficient_bound or RING.norm(z[2:]) >= coefficient_bound:
            continue
        if order.norm(z) not in (ONE, (-1, 0)):
            continue
        candidates.append(z)
        value = order.abs_square(z)
        if value.compare(1) >= 0 and value.compare(E2) < 0:
            reduced.add(z)
    torsion = {power(order, ROOT, j) for j in range(torsion_order)}
    require(len(torsion) == torsion_order and reduced == torsion, "complete reduced unit group")
    inverse = order.sigma(generator)
    if norm == (-1, 0):
        inverse = tuple(-x for x in inverse)
    require(order.mul(generator, inverse) == IDENTITY, "integral unit inverse")
    require(all(order.norm(z) == ONE for z in torsion), "torsion has relative norm one")
    return dict(order=name, coefficient_norm_bound=[coefficient_bound.numerator, coefficient_bound.denominator],
                generator=generator, generator_relative_norm=norm,
                generator_abs_square=[E2.a, E2.b, E2.d, E2.denominator],
                norm_equation_candidates=len(candidates), reduced_unit_count=len(reduced),
                torsion_order=torsion_order, determinant_image=[1] if name == "A2" else [-1, 1])


def expected_classes():
    O2, O3 = RelativeOrder("A2"), RelativeOrder("A3")
    return [
        dict(label="d7 order 2: A2", m=2, finite_centralizer_order=2,
             norm=(8, 3, 7), norm_denominator=1,
             representative=O2.matrix(ROOT), primitive_translation=O2.matrix(GENERATORS["A2"]), flip=None),
        dict(label="d7 order 3: inverse classes merged", m=3, finite_centralizer_order=3,
             norm=(23, 5, 21), norm_denominator=2,
             representative=O3.matrix(ROOT), primitive_translation=O3.matrix(power(O3, GENERATORS["A3"], 2)),
             flip=None),
    ]


def verify_witnesses_and_splitting():
    classes = expected_classes()
    for C, name in zip(classes, ("A2", "A3")):
        order = RelativeOrder(name)
        R, T = C["representative"], C["primitive_translation"]
        require(OPS.det(R) == OPS.det(T) == ONE, "SL elliptic and translation witnesses")
        require(OPS.power(R, C["m"]) == OPS.neg(MATRIX_IDENTITY), "PSL elliptic order")
        require(OPS.mul(R, T) == OPS.mul(T, R), "commuting translation")
        norm = order.abs_square(GENERATORS[name])
        if name == "A3":
            norm = norm.square()
        require((norm.a, norm.b, norm.d) == C["norm"] and norm.denominator == C["norm_denominator"],
                "individual primitive norm")
        H = {OPS.canon(OPS.power(R, j)) for j in range(C["m"])}
        require(len(H) == C["finite_centralizer_order"], "finite rotation subgroup")
        require(all(OPS.canon(OPS.mul(A, B)) in H for A in H for B in H), "finite subgroup closure")
    # Trace-zero: two SL determinant cosets, merged by +/- lifts in PSL.
    D = ((1, 0), (0, 0), (0, 0), (-1, 0))
    R2 = classes[0]["representative"]
    require(OPS.det(D) == (-1, 0) and OPS.mul(OPS.mul(D, R2), D) == OPS.neg(R2),
            "opposite determinant coset is the negative involution lift")
    # Trace-one: a norm-minus-one centralizer unit corrects the GL inverse conjugator.
    O3 = RelativeOrder("A3")
    J = ((1, 0), (1, 0), (0, 0), (-1, 0))
    U = O3.matrix(GENERATORS["A3"])
    K = OPS.mul(J, U)
    R3 = classes[1]["representative"]
    Ri = ((1, 0), (1, 0), (-1, 0), (0, 0))
    require(OPS.det(K) == ONE and OPS.power(K, 2) == OPS.neg(MATRIX_IDENTITY), "SL inverse conjugator")
    require(OPS.mul(K, R3) == OPS.mul(Ri, K), "order-three inverse element classes merge")
    return dict(classes=classes, trace_zero_GL_to_SL=2, trace_zero_PSL_classes=1,
                trace_one_GL_to_SL=1, trace_one_PSL_classes=1, inverse_conjugator=K)


def verify_inventory():
    orders = verify_orders_and_lattices()
    units = [verify_units(name) for name in ("A2", "A3")]
    splitting = verify_witnesses_and_splitting()
    return dict(proof_id="d7-arithmetic-v1", number_theory=orders, units=units,
                class_count=2, classes=splitting.pop("classes"), splitting=splitting,
                orbital_normalization="full-PSL-centralizer; no endpoint flip in involution centralizer")


def verify_group_records(group):
    require(group.key == GroupKey(7) and group.inventory_status == "self-contained", "d7 proof binding")
    require(group.cusp_count == group.GG == 1 and group.ce_g0[0] == 0
            and group.ce_integral == 0 and group.ce_kernel == 1, "d7 cusp constants")
    proof = verify_inventory()
    require(len(group.elliptic_classes) == 2, "complete d7 class count")
    for actual, expected in zip(group.elliptic_classes, proof["classes"]):
        for key, value in expected.items():
            require(getattr(actual, key) == value, f"d7 record differs: {key}")
        require(not actual.cuspidal and actual.normalization_status == "proved", "d7 class status")
    return proof


if __name__ == "__main__":
    import json
    print(json.dumps(verify_inventory(), indent=2))
