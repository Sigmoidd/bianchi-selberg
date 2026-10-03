"""Replay direct negative-index scalar arithmetic; never issue exclusion."""
from fractions import Fraction as F
import hashlib
import json
from reference_mesh import ROOT


def budget():
    path=ROOT/'tuned_error_bounds.json'
    data=json.loads(path.read_text())
    if (data['schema']!='d67-tuned-exact-scalar-ledger/v1' or data['d']!=67
            or data['vertical_energy_fraction']!='3/4'
            or data['spectral_exclusion_certified'] is not False):
        raise ArithmeticError('tuned scalar identity')
    for name,key in [('adaptive_refinement_plan.json.gz','plan_sha256'),
                     ('tuned_error_bounds.py','scalar_producer_sha256')]:
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=data[key]:
            raise ArithmeticError('tuned scalar binding: '+name)
    gamma=F(data['gamma_squared_upper']);eta=F(1,10)
    c0=1-(1+1/eta)*gamma
    return dict(schema='d67-direct-index-scalar-budget/v1',d=67,
        ledger_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        eta=str(eta),gamma_squared_upper=str(gamma),c0=str(c0),c0_float=float(c0),
        scalar_budget_arithmetic_pass=c0>=0,
        finite_matrix='Q_h + b_h/2 - 3*t_h*t_h^T/4 - 11*M_h/10',
        sufficient_check='finite_matrix + alpha*z*z^T is PSD for any rational z and alpha > 0',
        weighted_mean_enclosure_required=False,
        threshold_matrix_positivity_verified=False,spectral_exclusion_certified=False)


if __name__=='__main__':print(json.dumps(budget(),indent=2))
