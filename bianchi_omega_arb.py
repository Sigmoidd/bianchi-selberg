"""Compatibility entry point for the two historical level-one Arb runs.

The shared implementation now lives in core/, fields/, and groups/. The
Eisenstein elliptic normalization remains frozen, with the issue documented.
"""
from core.assemble import evaluate, integ, setp, PREC, RTOL
from core.bspline import gbspline, g0, gpp0
from core.testfunctions import sinc2k, h_isig, ce_integral
from fields.quadratic import Lval1, Lprime1, zetaK2, zetaK_logder, FpF, eta, _Lseries
from groups import get_group

CHI = {"i": (4, {1: 1, 3: -1}), "omega": (3, {1: 1, 2: -1})}

def field(kind):
    return get_group(kind).analytic_data()

def compute_B(kind, k=2, frac=0.999, R=40, verbose=True):
    return evaluate(kind, k, frac, R, verbose).bound

if __name__=="__main__":
    print("="*62); print("CERTIFIED (Arb) -- validate Q(i) [expect B ~ 0.31]:"); print("="*62)
    Bi=compute_B("i")
    print(f"  B(Q(i)) enclosure: [{float(Bi.lower()):.6f}, {float(Bi.upper()):.6f}]")
    print()
    print("="*62); print("Q(omega) (Eisenstein-Picard), level 1:"); print("="*62)
    print("Historical elliptic normalization is frozen; see docs/NORMALIZATION_ISSUE.md.")
    Bw=compute_B("omega")
    lo,hi=float(Bw.lower()),float(Bw.upper())
    print(f"  B(Q(omega)) enclosure: [{lo:.6f}, {hi:.6f}]")
    print()
    if hi<1:
        print(f"  *** CERTIFIED B < 1 (upper bound {hi:.6f}) => lambda_1(PSL(2,Z[omega])) >= 1 ***")
    else:
        print(f"  upper bound {hi:.6f} not < 1; tighten R/prec.")
