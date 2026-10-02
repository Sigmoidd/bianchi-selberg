"""Fail-closed status for the threshold route; not a matrix verifier.

This command deliberately cannot issue an exclusion certificate before an
independent adaptive mesh/moment-map and interval matrix verifier exists.
"""
from fractions import Fraction as F
import hashlib,json
from reference_mesh import ROOT
from verify_reference_mesh import verify_mesh


def status():
    mesh=verify_mesh()
    report=ROOT/'adaptive_error_bounds.json'
    if not report.exists():report=ROOT/'adaptive_error_bounds_200k.json'
    local=json.loads(report.read_text())
    if local['schema']!='d67-adaptive-error-bound-diagnostic/v1' or local['d']!=67:
        raise ArithmeticError('adaptive bound identity mismatch')
    if local['spectral_exclusion_certified'] is not False:
        raise ArithmeticError('adaptive local bounds cannot claim exclusion')
    gamma=F(local['gamma_squared_upper']);sigma=F(local['sigma_squared_upper'])
    eta=F(1,5);theta=F(9,10);rho=F(5)
    c_e=1-(1+1/eta)*gamma-rho*(1/theta-1)*sigma
    failures=[]
    if c_e<0:failures.append('SCALAR_ERROR_BUDGET_FAIL')
    failures.extend(['ADAPTIVE_FACE_MOMENT_MAP_NOT_VERIFIED',
                     'ENCLOSED_THRESHOLD_MATRICES_AND_MEAN_VECTOR_MISSING',
                     'VERIFIED_THRESHOLD_MATRIX_POSITIVITY_MISSING',
                     'INDEPENDENT_FULL_SPECTRAL_REPLAY_MISSING'])
    return dict(schema='d67-threshold-route-status/v1',d=67,subgroup='full',level_hnf=[1,0,1],target_interval='(0,1)',
        exact_reference_mesh=mesh,local_bound_report=report.name,local_bound_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),
        eta=str(eta),theta=str(theta),rho=str(rho),finite_mean_penalty=str(rho*(1-theta)),
        c_e=str(c_e),c_e_float=float(c_e),local_scalar_budget_arithmetic_pass=c_e>=0,
        failed_conditions=failures,spectral_exclusion_certified=False)

if __name__=='__main__':
    result=status();(ROOT/'status.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    raise SystemExit(2)
