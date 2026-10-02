"""Exact replay for docs/D19_INVENTORY_PROOF.md; no classification tables."""
from fractions import Fraction
from itertools import product
from math import gcd
from flint import arb
from groups.arithmetic import QuadraticRing
from groups.exact import Radical
from groups.identity import GroupKey
from groups.matrix import MatrixOps, IDENTITY as MATRIX_IDENTITY
from groups.relative_orders import QuadraticOrder, ROOT, IDENTITY, power

RING = QuadraticRing(19)
OPS = MatrixOps(RING)
ONE = (1, 0)
GENERATORS = {'A2': (5, 3, 8, -3), 'A3': (10, 0, -3, -4)}


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


class RelativeOrder(QuadraticOrder):
    def __init__(self, name):
        if name not in GENERATORS:
            raise ValueError('unknown d19 relative order')
        self.name = name
        super().__init__(RING, (0, 0) if name == 'A2' else ONE, ONE)

    def abs_square(self, z):
        a, b, c, d = z
        base = RING.norm((a, b)) + RING.norm((c, d))
        if self.name == 'A2':
            return Radical(base, b*c-a*d, 19)
        return Radical(2*base+2*a*c+a*d+b*c+10*b*d, b*c-a*d, 57, 2)


def residue(z, tau, root, p):
    a,b,c,d = z
    return (a+b*tau+(c+d*tau)*root) % p


def verify_orders_and_lattices():
    O2,O3 = RelativeOrder('A2'),RelativeOrder('A3')
    require(arb(2)*arb(19).sqrt()/arb.pi() < 3, 'base Minkowski bound below three')
    roots = lambda p: {x for x in range(p) if (x*x-x+5) % p == 0}
    require(roots(2) == set(), 'base prime two inert; base class number one')
    require(O2.discriminant() == 5776 and O3.discriminant() == 3249, 'relative discriminants')
    require(gcd(19,4) == gcd(19,3) == 1, 'coprime discriminants prove maximality')
    require(arb(114)/arb.pi()**2 < 12, 'A2 Minkowski bound below twelve')
    require(arb(171)/(2*arb.pi()**2) < 9, 'A3 Minkowski bound below nine')
    require(roots(3) == set() and roots(5) == {0,1} and roots(7) == {2,6}
            and roots(11) == {3,9}, 'complete small-prime base factorization')
    # A2 prime norms <=11: 4, four at 5, two at 9. None at 2,3,7,11.
    require(O2.absolute_norm((1,0,1,0)) == 4, 'principal norm-four prime')
    require({x for x in range(5) if (x*x+1)%5 == 0} == {2,3}, 'A2 factors at five')
    five = [(2,0,0,1),(2,0,0,-1),(2,0,1,-1),(2,0,-1,1)]
    pairs = list(product((0,1),(2,3)))
    vanishing = []
    for z in five:
        require(O2.absolute_norm(z) == 5, 'principal A2 norm-five prime')
        zeros = [pair for pair in pairs if residue(z,*pair,5) == 0]
        require(len(zeros) == 1, 'unique norm-five residue pair')
        vanishing.append(zeros[0])
    require(set(vanishing) == set(pairs), 'every A2 norm-five prime principal')
    f = lambda z,r,p: tuple(x%p for x in RING.add(z[:2],RING.mul(z[2:],r)))
    f9roots = [(a,b) for a,b in product(range(3),repeat=2)
               if tuple(x%3 for x in RING.add(RING.mul((a,b),(a,b)),ONE)) == (0,0)]
    require(f9roots == [(1,1),(2,2)], 'two A2 roots over F9')
    nine = [(2,1,3,-1),(2,1,-3,1)]
    require(all(O2.absolute_norm(z) == 9 for z in nine), 'principal norm-nine primes')
    require([[f(z,r,3)==(0,0) for r in f9roots] for z in nine]
            == [[False,True],[True,False]], 'distinct norm-nine prime factors')
    require(all(not any((x*x+1)%p == 0 for x in range(p)) for p in (7,11)),
            'no A2 prime norms seven or eleven')
    # A3 prime norms <=8: two at 4, four at 7. At 3 ramified over F9 -> norm9.
    four = [(1,1,2,-1),(3,0,-2,1)]
    require(all(O3.absolute_norm(z) == 4 for z in four), 'principal A3 norm-four primes')
    require([[f(z,r,2)==(0,0) for r in ((0,1),(1,1))] for z in four]
            == [[True,False],[False,True]], 'distinct A3 prime factors over F4')
    require({x for x in range(3) if (x*x-x+1)%3 == 0} == {2}, 'alpha ramified at three')
    require({x for x in range(5) if (x*x-x+1)%5 == 0} == set(), 'no A3 norm-five prime')
    require({x for x in range(7) if (x*x-x+1)%7 == 0} == {3,5}, 'four A3 primes at seven')
    seven = [(2,0,0,-1),(2,-1,0,1),(2,0,-1,1),(1,1,1,-1)]
    pairs = list(product((2,6),(3,5)))
    zeros = []
    for z in seven:
        require(O3.absolute_norm(z) == 7, 'principal A3 norm-seven prime')
        hit = [pair for pair in pairs if residue(z,*pair,7) == 0]
        require(len(hit) == 1, 'unique norm-seven residue pair')
        zeros.append(hit[0])
    require(set(zeros) == set(pairs), 'every A3 norm-seven prime principal')
    return dict(discriminants=dict(A2=5776,A3=3249), maximal_class_numbers=dict(A2=1,A3=1),
                GL_lattice_types=dict(trace_zero=['A2'],trace_one=['A3']),
                minkowski_cutoffs=dict(A2=11,A3=8),
                principal_prime_norms=dict(A2=[4,5,5,5,5,9,9],A3=[4,4,7,7,7,7]),
                A2_norm_five_generators=five,A2_norm_five_residues=vanishing,
                A2_norm_nine_generators=nine,A3_norm_four_generators=four,
                A3_norm_seven_generators=seven,A3_norm_seven_residues=zeros)


def verify_units(name):
    order=RelativeOrder(name)
    generator=GENERATORS[name]
    torsion_order=4 if name=='A2' else 6
    bound=Fraction(171,2) if name=='A2' else Fraction(304,3)
    E2=order.abs_square(generator)
    expected=Radical(170,39,19) if name=='A2' else Radical(302,40,57,2)
    require(E2 == expected and E2.compare(1)>0, 'expanding generator modulus')
    E2.unit_reciprocal()
    norm=order.norm(generator)
    require(norm == (-1,0), 'generator relative norm')
    require(Fraction(2*E2.a+2*E2.denominator,E2.denominator*(4 if name=='A2' else 3)) == bound,
            'proved archimedean coefficient bound')
    ar,br=(range(-11,12),range(-4,5)) if name=='A2' else (range(-12,13),range(-4,5))
    reduced,candidates=set(),[]
    for z in product(ar,br,ar,br):
        if RING.norm(z[:2])>=bound or RING.norm(z[2:])>=bound:
            continue
        if order.norm(z) not in (ONE,(-1,0)):
            continue
        candidates.append(z)
        value=order.abs_square(z)
        if value.compare(1)>=0 and value.compare(E2)<0:
            reduced.add(z)
    torsion={power(order,ROOT,j) for j in range(torsion_order)}
    require(len(torsion)==torsion_order and reduced==torsion, 'complete reduced unit group')
    require(len(candidates)==(12 if name=='A2' else 18), 'norm-equation candidate count')
    inverse=order.sigma(generator)
    if norm==(-1,0): inverse=tuple(-x for x in inverse)
    require(order.mul(generator,inverse)==IDENTITY, 'integral unit inverse')
    require(all(order.norm(z)==ONE for z in torsion), 'torsion relative norms')
    return dict(order=name,coefficient_norm_bound=[bound.numerator,bound.denominator],
                coordinate_bounds=dict(a=11 if name=='A2' else 12,b=4),
                generator=generator,generator_relative_norm=norm,
                generator_abs_square=[E2.a,E2.b,E2.d,E2.denominator],
                norm_equation_candidates=len(candidates),reduced_unit_count=len(reduced),
                torsion_order=torsion_order,determinant_image=[-1,1])


def expected_classes():
    O2,O3=RelativeOrder('A2'),RelativeOrder('A3')
    J2=((1,0),(0,0),(0,0),(-1,0))
    return [dict(label='d19 order 2: A2 with endpoint flip',m=2,finite_centralizer_order=4,
                 norm=(57799,13260,19),norm_denominator=1,representative=O2.matrix(ROOT),
                 primitive_translation=O2.matrix(power(O2,GENERATORS['A2'],2)),
                 flip=OPS.mul(J2,O2.matrix(GENERATORS['A2']))),
            dict(label='d19 order 3: inverse classes merged',m=3,finite_centralizer_order=3,
                 norm=(45601,6040,57),norm_denominator=1,representative=O3.matrix(ROOT),
                 primitive_translation=O3.matrix(power(O3,GENERATORS['A3'],2)),flip=None)]


def verify_witnesses_and_splitting():
    classes=expected_classes()
    O2,O3=RelativeOrder('A2'),RelativeOrder('A3')
    require(O2.abs_square(GENERATORS['A2']).square()==Radical(57799,13260,19), 'primitive involution norm')
    require(O3.abs_square(GENERATORS['A3']).square().compare(Radical(45601,6040,57))==0,
            'primitive order-three norm')
    for C in classes:
        R,T=C['representative'],C['primitive_translation']
        require(OPS.det(R)==OPS.det(T)==ONE, 'SL elliptic and translation witnesses')
        require(OPS.power(R,C['m'])==OPS.neg(MATRIX_IDENTITY), 'PSL element order')
        require(OPS.mul(R,T)==OPS.mul(T,R), 'commuting translation')
        if C['m']==2:
            X=C['flip']
            require(OPS.det(X)==ONE and OPS.power(X,2)==OPS.neg(MATRIX_IDENTITY), 'SL endpoint flip')
            require(OPS.mul(X,R)==OPS.neg(OPS.mul(R,X)), 'PSL flip commutation')
            H={OPS.canon(M) for M in (MATRIX_IDENTITY,R,X,OPS.mul(R,X))}
        else:
            H={OPS.canon(OPS.power(R,j)) for j in range(3)}
        require(len(H)==C['finite_centralizer_order'], 'finite full-centralizer subgroup size')
        require(all(OPS.canon(OPS.mul(A,B)) in H for A in H for B in H), 'finite subgroup closure')
    J3=((1,0),(1,0),(0,0),(-1,0))
    X3=OPS.mul(J3,O3.matrix(GENERATORS['A3']))
    R3=classes[1]['representative']
    require(OPS.det(X3)==ONE and OPS.power(X3,2)==OPS.neg(MATRIX_IDENTITY), 'SL inverse conjugator')
    require(OPS.mul(OPS.mul(X3,R3),OPS.neg(X3))==OPS.mul(OPS.mul(J3,R3),J3),
            'inverse order-three classes merge in SL')
    return dict(classes=classes,trace_zero_GL_to_SL=1,trace_zero_PSL_classes=1,
                trace_one_GL_to_SL=1,trace_one_PSL_classes=1,inverse_conjugator=X3,
                inverse_conjugator_determinant=1,order_three_inverse_classes_merge=True)


def verify_inventory():
    orders=verify_orders_and_lattices()
    units=[verify_units(name) for name in ('A2','A3')]
    split=verify_witnesses_and_splitting()
    return dict(proof_id='d19-arithmetic-v1',number_theory=orders,units=units,
                class_count=2,classes=split.pop('classes'),splitting=split,
                orbital_normalization='full-PSL-centralizer; involution endpoint flip included')


def verify_group_records(group):
    require(group.key==GroupKey(19) and group.inventory_status=='self-contained', 'd19 proof binding')
    require(group.cusp_count==group.GG==1 and group.ce_g0[0]==0 and group.ce_integral==0
            and group.ce_kernel==1, 'd19 cusp constants')
    proof=verify_inventory()
    require(len(group.elliptic_classes)==2, 'complete d19 class count')
    for actual,expected in zip(group.elliptic_classes,proof['classes']):
        for key,value in expected.items():
            require(getattr(actual,key)==value, f'd19 record differs: {key}')
        require(not actual.cuspidal and actual.normalization_status=='proved', 'd19 class status')
    return proof


if __name__=='__main__':
    import json
    print(json.dumps(verify_inventory(),indent=2))
