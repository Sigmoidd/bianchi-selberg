"""Exact replay of the arithmetic proof in docs/D2_INVENTORY_PROOF.md.

No bounded matrix-conjugator search, numerical root classification, or
external class-count table is used. Finite unit boxes have proved bounds.
"""
from itertools import product
from math import gcd

from flint import arb, fmpz_mat
from groups.arithmetic import QuadraticRing
from groups.exact import Radical
from groups.matrix import MatrixOps
from groups.relative_orders import QuadraticOrder, BASIS, IDENTITY, ROOT, power

RING = QuadraticRing(2)
ZERO = (0, 0)
ONE = (1, 0)
S = (0, 1)
I2 = (ONE, ZERO, ZERO, ONE)


class RelativeOrder(QuadraticOrder):
    """d=2 orders with their explicit archimedean modulus formulas."""
    PARAMETERS = {"A2": (ZERO, ONE), "S2": (S, (-1, 0)), "A3": (ONE, ONE)}

    def __init__(self, name):
        self.name = name
        super().__init__(RING, *self.PARAMETERS[name])

    def abs_square(self, z):
        a, b, c, d = z
        base = a*a+2*b*b+c*c+2*d*d
        if self.name == "A2":
            return Radical(base, 2*(b*c-a*d), 2)
        if self.name == "S2":
            return Radical(base+2*(b*c-a*d), a*c+2*b*d, 2)
        return Radical(base+a*c+2*b*d, b*c-a*d, 6)

OPS = MatrixOps(RING)
mat_mul, mat_neg, mat_det, mat_canon = OPS.mul, OPS.neg, OPS.det, OPS.canon


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def verify_orders_and_class_numbers():
    orders = {name: RelativeOrder(name) for name in RelativeOrder.PARAMETERS}
    discriminants = {name: order.discriminant() for name, order in orders.items()}
    require(discriminants == dict(A2=1024, S2=256, A3=576), "order discriminants")
    O8 = orders["S2"]
    require(power(O8, ROOT, 4) == (-1, 0, 0, 0), "zeta8 polynomial")
    columns = [power(O8, ROOT, j) for j in range(4)]
    change = fmpz_mat([[columns[j][i] for j in range(4)] for i in range(4)])
    require(abs(int(change.det())) == 1, "S2 equals the zeta8 power order")
    # (y+1)^4+1 is Eisenstein at 2. Together with discr=2^8 this
    # proves maximality, as explained in the local valuation argument.
    eisenstein = (2, 4, 6, 4, 1)
    require(eisenstein[0] % 4 == 2 and all(c % 2 == 0 for c in eisenstein[:-1]),
            "Eisenstein uniformizer proof")
    pi = (1, 0, -1, 0)
    require(O8.absolute_norm(pi) == 2, "principal norm-two ideal")
    require(arb(24)/arb.pi()**2 < 3, "zeta8 Minkowski bound")
    # A3 is the compositum of maximal quadratic orders of coprime
    # discriminants -8 and -3. Their tensor product is maximal locally.
    require(gcd(8, 3) == 1, "coprime quadratic discriminants")
    require(arb(36)/arb.pi()**2 < 4, "A3 Minkowski bound")
    O3 = orders["A3"]
    beta, beta_bar = (0, 1, 1, 0), (0, -1, 1, 0)  # s+alpha, -s+alpha
    require(O3.absolute_norm(beta) == O3.absolute_norm(beta_bar) == 3,
            "both norm-three ideals have principal generators")
    # At 2, alpha^2-alpha+1 has no F2 root: residue degree >=2,
    # so there is no norm-two ideal. At 3, s=+/-1 and alpha=-1.
    require(all((r*r-r+1) % 2 != 0 for r in (0, 1)), "no norm-two prime in A3")
    require({r for r in range(3) if (r*r+2) % 3 == 0} == {1, 2}, "two primes at 3")
    require({r for r in range(3) if (r*r-r+1) % 3 == 0} == {2}, "alpha residue at 3")
    require((1+2) % 3 == 0 and (-1+2) % 3 != 0, "beta selects one prime")
    require((-(-1)+2) % 3 == 0 and (-(1)+2) % 3 != 0, "beta_bar selects the other")
    return dict(discriminants=discriminants, maximal_class_numbers=dict(S2=1, A3=1),
                principal_prime_norms=dict(S2=[2], A3=[3, 3]))


def verify_lattice_types():
    O8 = RelativeOrder("S2")
    # A2 embeds by i=1+s*z, s*i=s-2*z. Its index in S2 is two.
    A_basis = [IDENTITY, BASIS[1], (1, 0, 0, 1), (0, 1, -2, 0)]
    inclusion = fmpz_mat([[A_basis[j][i] for j in range(4)] for i in range(4)])
    require(abs(int(inclusion.det())) == 2, "A2 index")
    s = BASIS[1]
    f_basis = [O8.mul(s, e) for e in BASIS]
    conductor = fmpz_mat([[f_basis[j][i] for j in range(4)] for i in range(4)])
    require(abs(int(conductor.det())) == 4, "conductor norm")
    reduce = lambda z: (z[0] % 2)+2*(z[2] % 2)
    require({reduce(e) for e in A_basis} == {0, 1}, "A2 image in S2/f")
    # S2/f = F2[z]/((z+1)^2), with integers 0,1,z,1+z.
    def mul(x, y):
        a, b, c, d = x % 2, x//2, y % 2, y//2
        return ((a*c+b*d) % 2)+2*((a*d+b*c) % 2)
    subspaces = [frozenset([0]), frozenset([0, 1]), frozenset([0, 2]),
                 frozenset([0, 3]), frozenset(range(4))]
    normalized = []
    for V in subspaces:
        span = {mul(x, v) for x in range(4) for v in V}
        # Close under addition (XOR).
        span |= {x ^ y for x in span for y in span}
        if len(span) == 4:
            normalized.append(V)
    require(normalized == [frozenset([0, 1]), frozenset([0, 2]), frozenset(range(4))],
            "all normalized conductor subspaces")
    require({mul(2, v) for v in normalized[0]} == set(normalized[1]),
            "the integral unit z swaps the two unit lines")
    require(O8.norm(ROOT) == (-1, 0), "z is an invertible integral multiplier")
    return dict(trace_zero_GL_types=["A2", "S2"], trace_one_GL_types=["A3"],
                conductor_quotient_size=4, normalized_subspaces=3)


def verify_unit_reduction(name):
    order = RelativeOrder(name)
    generator = {"A2": (1, 0, 0, -1), "S2": (1, -1, 2, 0),
                 "A3": (-1, 1, 2, 0)}[name]
    torsion_order = {"A2": 4, "S2": 8, "A3": 6}[name]
    root_difference_squared = {"A2": 4, "S2": 2, "A3": 3}[name]
    generator_norm = order.norm(generator)
    require(generator_norm == ((1, 0) if name == "A3" else (-1, 0)), "unit generator norm")
    E2 = order.abs_square(generator)
    require(E2.compare(1) > 0, "unit generator expands")
    require(E2.a*E2.a-E2.d*E2.b*E2.b == 1, "reciprocal absolute modulus")
    require((2*E2.a+2) % root_difference_squared == 0, "exact coefficient bound")
    coefficient_bound = (2*E2.a+2)//root_difference_squared
    require(coefficient_bound in (2, 4), "proved small unit box")
    reduced = set()
    admissible = []
    for z in product(range(-1, 2), repeat=4):
        a, b, c, d = z
        if a*a+2*b*b >= coefficient_bound or c*c+2*d*d >= coefficient_bound:
            continue
        if order.norm(z) not in (ONE, (-1, 0)):
            continue
        admissible.append(z)
        value = order.abs_square(z)
        if value.compare(1) >= 0 and value.compare(E2) < 0:
            reduced.add(z)
    torsion = {power(order, ROOT, j) for j in range(torsion_order)}
    require(len(torsion) == torsion_order and power(order, ROOT, torsion_order) == IDENTITY,
            "torsion generator and distinct powers")
    require(reduced == torsion, "exhaustive reduced unit list equals the stated roots of unity")
    inverse = order.sigma(generator)
    if generator_norm == (-1, 0):
        inverse = tuple(-c for c in inverse)
    require(order.mul(generator, inverse) == IDENTITY, "generator has integral inverse")
    return dict(order=name, reduced_unit_count=len(reduced),
                norm_equation_candidates=len(admissible), coefficient_norm_bound=coefficient_bound,
                generator=list(generator), generator_relative_norm=list(generator_norm),
                generator_abs_square=(E2.a, E2.b, E2.d),
                torsion_order=torsion_order,
                determinant_image=[1] if name == "A3" else [-1, 1])


def expected_classes():
    A2, S2, A3 = (RelativeOrder(name) for name in ("A2", "S2", "A3"))
    RA = A2.matrix(ROOT)
    RS = S2.matrix((1, 0, 0, 1))  # i=1+s*z
    RC = A3.matrix(ROOT)
    RC_inverse = (ONE, ONE, (-1, 0), ZERO)
    TA = A2.matrix(power(A2, (1, 0, 0, -1), 2))
    TS = S2.matrix(S2.mul(ROOT, (1, -1, 2, 0)))
    TC = A3.matrix((-1, 1, 2, 0))
    return [
        dict(label="d2 order 2: multiplier A2", m=2, finite_centralizer_order=4,
             norm=(17, 12, 2), representative=RA, primitive_translation=TA, flip=RS),
        dict(label="d2 order 2: multiplier S2", m=2, finite_centralizer_order=4,
             norm=(3, 2, 2), representative=RS, primitive_translation=TS, flip=RA),
        dict(label="d2 order 3: alpha", m=3, finite_centralizer_order=3,
             norm=(5, 2, 6), representative=RC, primitive_translation=TC, flip=None),
        dict(label="d2 order 3: alpha inverse", m=3, finite_centralizer_order=3,
             norm=(5, 2, 6), representative=RC_inverse,
             primitive_translation=A3.matrix((1, 1, -2, 0)), flip=None),
    ]


def verify_class_witnesses():
    classes = expected_classes()
    A2, S2, A3 = (RelativeOrder(name) for name in ("A2", "S2", "A3"))
    units = [power(A2, (1, 0, 0, -1), 2), S2.mul(ROOT, (1, -1, 2, 0)),
             (-1, 1, 2, 0), (1, 1, -2, 0)]
    for C, order, unit in zip(classes, (A2, S2, A3, A3), units):
        modulus = order.abs_square(unit)
        # The inverse class uses the contracting embedding; its reciprocal
        # is the expanding eigenvalue modulus squared.
        if modulus.compare(1) < 0:
            require(modulus.a**2-modulus.d*modulus.b**2 == 1, "reciprocal norm")
            modulus = Radical(modulus.a, -modulus.b, modulus.d)
        require((modulus.a, modulus.b, modulus.d) == C["norm"], "primitive norm witness")
        require(order.matrix(unit) == C["primitive_translation"], "primitive unit matrix")
    for C in classes:
        R, T = C["representative"], C["primitive_translation"]
        require(mat_det(R) == mat_det(T) == ONE, "SL representative and primitive translation")
        require(mat_mul(T, R) == mat_mul(R, T), "translation commutes")
        square = mat_mul(R, R)
        if C["m"] == 2:
            require(square == mat_neg(I2), "PSL order two")
            X = C["flip"]
            require(mat_det(X) == ONE and mat_mul(X, X) == mat_neg(I2), "SL flip involution")
            require(mat_mul(X, R) == mat_neg(mat_mul(R, X)), "PSL flip commutation")
            H = {mat_canon(Y) for Y in (I2, R, X, mat_mul(R, X))}
        else:
            require(mat_mul(square, R) == mat_neg(I2), "PSL order three")
            H = {mat_canon(Y) for Y in (I2, R, square)}
        require(len(H) == C["finite_centralizer_order"], "finite subgroup size")
        require(all(mat_canon(mat_mul(A, B)) in H for A in H for B in H), "finite subgroup closure")
    # The two trace-zero GL types cannot merge: scalar/non-scalar mod (s).
    reduced = [tuple(a % 2 for a, b in C["representative"]) for C in classes[:2]]
    require(reduced == [(0, 1, 1, 0), (1, 0, 0, 1)], "distinct involution multiplier types")
    # Relative conjugation on A3 has determinant -1 and sends alpha to alpha^-1.
    J = (ONE, ONE, ZERO, (-1, 0))
    require(mat_det(J) == (-1, 0) and mat_mul(J, J) == I2, "trace-one GL inverse witness")
    require(mat_mul(mat_mul(J, classes[2]["representative"]), J)
            == classes[3]["representative"], "inverse classes differ by determinant coset")
    return classes


def verify_inventory():
    number_theory = verify_orders_and_class_numbers()
    lattices = verify_lattice_types()
    units = [verify_unit_reduction(name) for name in ("A2", "S2", "A3")]
    classes = verify_class_witnesses()
    require([c["m"] for c in classes] == [2, 2, 3, 3], "complete PSL element class multiplicity")
    return dict(proof_id="d2-arithmetic-v1", number_theory=number_theory,
                lattices=lattices, units=units, class_count=4, classes=classes,
                orbital_normalization="full-PSL-centralizer")


def verify_group_records(group):
    from groups.identity import GroupKey
    require(group.key == GroupKey(2) and group.inventory_status == "self-contained", "d2 proof binding")
    require(group.cusp_count == group.GG == 1 and group.ce_g0[0] == 0
            and group.ce_integral == 0 and group.ce_kernel == 1, "d2 cusp inputs")
    proof = verify_inventory()
    require(len(group.elliptic_classes) == len(proof["classes"]), "inventory record count")
    for actual, expected in zip(group.elliptic_classes, proof["classes"]):
        require(actual.norm_denominator == 1, "d2 primitive norm denominator")
        for name in ("label", "m", "finite_centralizer_order", "norm", "representative",
                     "primitive_translation", "flip"):
            require(getattr(actual, name) == expected[name], f"inventory record differs: {name}")
        require(not actual.cuspidal and actual.normalization_status == "proved", "class proof status")
    return proof


if __name__ == "__main__":
    import json
    print(json.dumps(verify_inventory(), indent=2))
