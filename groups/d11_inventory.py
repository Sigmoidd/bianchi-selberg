"""Exact replay for docs/D11_INVENTORY_PROOF.md; no classification tables."""
from fractions import Fraction
from itertools import product
from math import gcd
from flint import arb
from groups.arithmetic import QuadraticRing
from groups.exact import Radical
from groups.identity import GroupKey
from groups.matrix import MatrixOps, IDENTITY as MATRIX_IDENTITY
from groups.relative_orders import QuadraticOrder, ROOT, IDENTITY, power

RING = QuadraticRing(11)
OPS = MatrixOps(RING)
ONE = (1, 0)
GENERATORS = {'A2': (1, 1, 2, -1), 'A3': (-3, 2, 4, 0)}


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


class RelativeOrder(QuadraticOrder):
    def __init__(self, name):
        if name not in GENERATORS:
            raise ValueError('unknown d11 relative order')
        self.name = name
        super().__init__(RING, (0, 0) if name == 'A2' else ONE, ONE)

    def abs_square(self, z):
        a, b, c, d = z
        base = RING.norm((a, b)) + RING.norm((c, d))
        if self.name == 'A2':
            return Radical(base, b*c-a*d, 11)
        return Radical(2*base+2*a*c+a*d+b*c+6*b*d, b*c-a*d, 33, 2)


def residue(z, tau, root, p):
    a,b,c,d = z
    return (a+b*tau+(c+d*tau)*root) % p


def verify_orders_and_lattices():
    O2,O3 = RelativeOrder('A2'),RelativeOrder('A3')
    require(O2.discriminant() == 1936 and O3.discriminant() == 1089, 'relative discriminants')
    require(gcd(11,4) == gcd(11,3) == 1, 'coprime discriminants prove maximality')
    require(arb(66)/arb.pi()**2 < 7, 'A2 Minkowski bound below seven')
    require(arb(99)/(2*arb.pi()**2) < 6, 'A3 Minkowski bound below six')
    roots = lambda p: {x for x in range(p) if (x*x-x+3) % p == 0}
    require(roots(2) == set(), 'base prime two is inert')
    require(roots(3) == {0,1}, 'base prime three splits')
    require(roots(5) == {2,4}, 'base prime five splits')
    # A2: one ramified prime at 2, residue F4; no prime of norm 2 or 3.
    require(O2.absolute_norm((1,0,1,0)) == 4, 'principal norm-four prime')
    require({x for x in range(3) if (x*x+1) % 3 == 0} == set(), 'no A2 norm-three prime')
    require({x for x in range(5) if (x*x+1) % 5 == 0} == {2,3}, 'four A2 primes at five')
    generators = [(0,1,1,0),(0,1,-1,0),(1,-1,1,0),(1,-1,-1,0)]
    pairs = list(product((2,4),(2,3)))
    vanishing = []
    for z in generators:
        require(O2.absolute_norm(z) == 5, 'principal A2 norm-five prime')
        zeros = [pair for pair in pairs if residue(z,*pair,5) == 0]
        require(len(zeros) == 1, 'unique norm-five residue pair')
        vanishing.append(zeros[0])
    require(set(vanishing) == set(pairs), 'every A2 norm-five prime principal')
    # A3: two unramified norm-four primes at 2; two ramified norm-three at 3.
    # In F4 tau^2=tau+1: alpha=tau or tau+1. These elements vanish separately.
    four = [(0,1,-1,0),(-1,1,1,0)]
    require(all(O3.absolute_norm(z) == 4 for z in four), 'principal A3 norm-four primes')
    f4 = lambda z,root: tuple(x % 2 for x in RING.add(z[:2],RING.mul(z[2:],root)))
    require([[f4(z,root) for root in ((0,1),(1,1))] for z in four]
            == [[(0,0),(1,0)],[(1,0),(0,0)]], 'distinct A3 norm-four residue factors')
    require({x for x in range(3) if (x*x-x+1) % 3 == 0} == {2}, 'alpha ramified at three')
    three = [(-2,0,0,1),(-2,0,1,-1)]
    zeros = []
    for z in three:
        require(O3.absolute_norm(z) == 3, 'principal A3 norm-three prime')
        zeros.append([t for t in (0,1) if residue(z,t,2,3) == 0])
    require(zeros == [[1],[0]], 'distinct A3 norm-three primes')
    require({x for x in range(5) if (x*x-x+1) % 5 == 0} == set(), 'no A3 norm-five prime')
    return dict(discriminants=dict(A2=1936,A3=1089), maximal_class_numbers=dict(A2=1,A3=1),
                GL_lattice_types=dict(trace_zero=['A2'],trace_one=['A3']),
                minkowski_cutoffs=dict(A2=6,A3=5),
                principal_prime_norms=dict(A2=[4,5,5,5,5],A3=[3,3,4,4]),
                A2_norm_five_generators=generators,A2_norm_five_residues=vanishing,
                A3_norm_three_generators=three,A3_norm_four_generators=four)


def verify_units(name):
    order=RelativeOrder(name)
    generator=GENERATORS[name]
    torsion_order=4 if name=='A2' else 6
    bound=Fraction(11,2) if name=='A2' else Fraction(16)
    E2=order.abs_square(generator)
    expected=Radical(10,3,11) if name=='A2' else Radical(46,8,33,2)
    require(E2 == expected and E2.compare(1)>0, 'expanding generator modulus')
    E2.unit_reciprocal()
    norm=order.norm(generator)
    require(norm == ((-1,0) if name=='A2' else ONE), 'generator relative norm')
    require(Fraction(2*E2.a+2*E2.denominator,E2.denominator*(4 if name=='A2' else 3)) == bound,
            'proved archimedean coefficient bound')
    ar,br=(range(-2,3),range(-1,2)) if name=='A2' else (range(-4,5),range(-2,3))
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
    require(len(candidates)==(12 if name=='A2' else 10), 'norm-equation candidate count')
    inverse=order.sigma(generator)
    if norm==(-1,0): inverse=tuple(-x for x in inverse)
    require(order.mul(generator,inverse)==IDENTITY, 'integral unit inverse')
    require(all(order.norm(z)==ONE for z in torsion), 'torsion relative norms')
    return dict(order=name,coefficient_norm_bound=[bound.numerator,bound.denominator],
                coordinate_bounds=dict(a=2 if name=='A2' else 4,b=1 if name=='A2' else 2),
                generator=generator,generator_relative_norm=norm,
                generator_abs_square=[E2.a,E2.b,E2.d,E2.denominator],
                norm_equation_candidates=len(candidates),reduced_unit_count=len(reduced),
                torsion_order=torsion_order,determinant_image=[-1,1] if name=='A2' else [1])


def expected_classes():
    O2,O3=RelativeOrder('A2'),RelativeOrder('A3')
    J2=((1,0),(0,0),(0,0),(-1,0))
    J3=((1,0),(1,0),(0,0),(-1,0))
    R2,R3=O2.matrix(ROOT),O3.matrix(ROOT)
    U2=O2.matrix(GENERATORS['A2'])
    T2=O2.matrix(power(O2,GENERATORS['A2'],2))
    T3=O3.matrix(GENERATORS['A3'])
    return [dict(label='d11 order 2: A2 with endpoint flip',m=2,finite_centralizer_order=4,
                 norm=(199,60,11),norm_denominator=1,representative=R2,
                 primitive_translation=T2,flip=OPS.mul(J2,U2)),
            dict(label='d11 order 3: alpha',m=3,finite_centralizer_order=3,
                 norm=(23,4,33),norm_denominator=1,representative=R3,primitive_translation=T3,flip=None),
            dict(label='d11 order 3: alpha inverse',m=3,finite_centralizer_order=3,
                 norm=(23,4,33),norm_denominator=1,representative=OPS.mul(OPS.mul(J3,R3),J3),
                 primitive_translation=OPS.mul(OPS.mul(J3,T3),J3),flip=None)]


def verify_witnesses_and_splitting():
    classes=expected_classes()
    O2,O3=RelativeOrder('A2'),RelativeOrder('A3')
    n2=O2.abs_square(GENERATORS['A2']).square()
    n3=O3.abs_square(GENERATORS['A3'])
    require(n2==Radical(199,60,11), 'primitive involution norm')
    require(n3.compare(Radical(23,4,33))==0, 'primitive order-three norm')
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
    require(OPS.det(J3)==(-1,0) and OPS.power(J3,2)==MATRIX_IDENTITY, 'GL inverse intertwiner')
    require(OPS.mul(OPS.mul(J3,classes[1]['representative']),J3)==classes[2]['representative'], 'inverse GL conjugacy')
    # No determinant-minus-one A3 unit exists, by the complete unit proof.
    return dict(classes=classes,trace_zero_GL_to_SL=1,trace_zero_PSL_classes=1,
                trace_one_GL_to_SL=2,trace_one_PSL_classes=2,inverse_conjugator=J3,
                inverse_conjugator_determinant=-1,order_three_inverse_classes_merge=False)


def verify_inventory():
    orders=verify_orders_and_lattices()
    units=[verify_units(name) for name in ('A2','A3')]
    split=verify_witnesses_and_splitting()
    return dict(proof_id='d11-arithmetic-v1',number_theory=orders,units=units,
                class_count=3,classes=split.pop('classes'),splitting=split,
                orbital_normalization='full-PSL-centralizer; involution endpoint flip included')


def verify_group_records(group):
    require(group.key==GroupKey(11) and group.inventory_status=='self-contained', 'd11 proof binding')
    require(group.cusp_count==group.GG==1 and group.ce_g0[0]==0 and group.ce_integral==0
            and group.ce_kernel==1, 'd11 cusp constants')
    proof=verify_inventory()
    require(len(group.elliptic_classes)==3, 'complete d11 class count')
    for actual,expected in zip(group.elliptic_classes,proof['classes']):
        for key,value in expected.items():
            require(getattr(actual,key)==value, f'd11 record differs: {key}')
        require(not actual.cuspidal and actual.normalization_status=='proved', 'd11 class status')
    return proof


if __name__=='__main__':
    import json
    print(json.dumps(verify_inventory(),indent=2))
