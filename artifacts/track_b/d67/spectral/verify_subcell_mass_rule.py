"""Exact partition and CR Gram identities, with a floating-probe cross-check.

This verifies the integration rule, not the full coefficient/matrix ledger.
"""
from fractions import Fraction as F
import json
from adaptive_error_bounds import subdivide
from coefficient_bounds import volume
from adaptive_matrix_probe import integration_rule
from reference_mesh import ROOT


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def barycentric(point):
    return (1-sum(point), *point)


def gram(tet):
    values=[tuple(1-3*x for x in barycentric(p)) for p in tet]
    total=[sum(row[j] for row in values) for j in range(4)]
    return tuple(tuple((total[i]*total[j]+sum(row[i]*row[j] for row in values))/20
                       for j in range(4)) for i in range(4))


def verify(depth):
    base=((F(0),F(0),F(0)),(F(1),F(0),F(0)),
          (F(0),F(1),F(0)),(F(0),F(0),F(1)))
    cells=[base]
    for _ in range(depth):
        cells=[child for parent in cells for child in subdivide(parent)]
    weights, floating_grams=integration_rule(depth)
    require(len(cells)==8**depth, 'subcell count')
    accumulated=[[F(0) for _ in range(4)] for _ in range(4)]
    maximum_error=0.0
    for k,cell in enumerate(cells):
        require(volume(cell)==volume(base)/8**depth, 'subcell volume')
        exact=gram(cell)
        for v,point in enumerate(cell):
            require(tuple(F(float(x)) for x in weights[k,v])==barycentric(point),
                    'floating barycentric rule differs from exact partition')
        for i in range(4):
            for j in range(4):
                accumulated[i][j]+=exact[i][j]/8**depth
                maximum_error=max(maximum_error,abs(float(exact[i][j])-floating_grams[k,i,j]))
    require(tuple(map(tuple, accumulated))==gram(base), 'constant-density Gram conservation')
    require(sum(map(volume,cells),F(0))==volume(base), 'volume conservation')
    require(maximum_error<=2e-16, 'floating Gram cross-check')
    return dict(depth=depth,subcells=len(cells),exact_volume_conservation=True,
                exact_constant_density_gram_conservation=True,
                floating_gram_maximum_difference=maximum_error)


if __name__=='__main__':
    result=dict(schema='d67-subcell-mass-rule-check/v1',d=67,
                checks=[verify(depth) for depth in range(4)],
                full_finite_mass_enclosure_verified=False,
                matrix_positivity_verified=False,spectral_exclusion_certified=False)
    (ROOT/'subcell_mass_rule_check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
