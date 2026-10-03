"""Exact locally balanced enclosure and its comparison scalar budget."""
from fractions import Fraction as F
from math import isqrt
import hashlib,json
from reference_mesh import ROOT


def balance(amin,bmin,radius_squared):
    t=amin*radius_squared/bmin
    if not t:return amin,bmin,F()
    scaled=t*2**40
    n=isqrt(scaled.numerator//scaled.denominator)
    if n*n*scaled.denominator<scaled.numerator:n+=1
    p=min(F(1,2),F(n,2**20))
    reserve=p*bmin
    return amin*reserve/(reserve+amin*radius_squared),(1-p)*bmin,p


def budget(eta=F(1,10)):
    if eta<=0:raise ValueError("eta must be positive")
    old=ROOT/'adaptive_error_bounds.json'
    data=json.loads(old.read_text())
    if data['d']!=67 or data['spectral_exclusion_certified'] is not False:
        raise ArithmeticError('source scalar identity')
    # Centroid improvement times an exact uniform inverse-metric comparison.
    ratio=F(4983,8000)*F(9,8)
    gamma=ratio*F(data['gamma_squared_upper'])
    c0=1-(1+1/eta)*gamma
    return dict(schema='d67-balanced-envelope-scalar-budget/v1',d=67,
        source_bound_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),
        kappa_squared_upper='1661/15000',inverse_metric_comparison='9/8',
        combined_ratio=str(ratio),gamma_squared_upper=str(gamma),gamma_squared_float=float(gamma),
        eta=str(eta),c0=str(c0),c0_float=float(c0),scalar_budget_arithmetic_pass=c0>=0,
        threshold_matrix_positivity_verified=False,spectral_exclusion_certified=False)


if __name__=='__main__':print(json.dumps(budget(),indent=2))
