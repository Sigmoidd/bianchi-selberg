"""Level-one assembly with explicit group inputs and inventory gates.

Mechanical screens omit elliptic terms; they are not spectral certificates.
"""
from dataclasses import dataclass
import math
from flint import arb, acb, ctx
from groups import get_group
from core.bspline import g0, gpp0, validate_k
from core.testfunctions import sinc2k, h_isig, ce_integral
from core.terms import prime_term
from groups.systoles import verify_systole

PREC = 100
RTOL = 1e-5

def setp():
    ctx.prec = PREC

@dataclass
class Evaluation:
    group: object
    k: int
    delta: float
    R: float
    include_elliptic: bool
    precision_bits: int
    bound: object
    terms: dict

def integ(f,a,b,panel=4.0):
    a=arb(a); b=arb(b); n=max(1,int(float((b-a)/panel)+0.5)); step=(b-a)/n
    tot=acb(0)
    for i in range(n):
        tot += acb.integral(lambda z,_: f(z), a+i*step, a+(i+1)*step, rel_tol=RTOL)
    return tot

def evaluate(kind, k=2, frac=0.999, R=40, verbose=True, include_elliptic=True):
    validate_k(k)
    if not (arb(frac) > 0 and arb(frac) <= 1):
        raise ValueError("support fraction must lie in (0, 1]")
    if not math.isfinite(float(R)) or R < 7:
        raise ValueError("shifted-line tail bound requires finite R >= 7")
    setp()
    group = get_group(kind)
    if group.cusp_count != 1:
        raise ValueError("this assembly requires exactly one cusp")
    # Ensure an arbitrary field label never substitutes an unproved systole.
    verify_systole(group)
    F = group.analytic_data(require_inventory=include_elliptic)
    # Float delta is an exact binary input. Rounding is followed by an Arb check.
    d = math.nextafter(float(F["systole"].lower())*float(frac)/(2*k), -math.inf)
    supp_ball = 2*k*arb(d)
    if not (arb(d) > 0 and supp_ball.upper() <= F["systole"].lower()):
        raise ValueError("support 2k*delta is not provably within the systole")
    supp = 2*k*d
    G0=g0(k,d)
    # I
    I = F["vol"]/(2*arb.pi())*(-gpp0(k,d))
    # NCE
    NCE = G0*F["C_ell"] if include_elliptic else arb(0)
    # CE integral (holomorphic B-spline polynomial per knot-interval)
    CEi = ce_integral(k,d,supp,F["ckern"]) if include_elliptic else arb(0)
    CE = F["CEg0"]*G0 + F["CEint"]*CEi if include_elliptic else arb(0)
    # SCATT + PAR h(0)
    Ch0 = arb(1)/4 + (arb(1)/F["GG"])*(arb(1)/4)
    # PAR g0
    PARg0 = (arb(1)/F["GG"])*G0*(F["eta"]/2 - arb.const_euler())
    # PSI = -(1/GG)(1/2pi) 2 int_0^inf h Re psi(1+ir) dr
    psi_main = 2*integ(lambda z: sinc2k(z,k,d)*(1+acb(0,1)*z).digamma(), 0, R).real
    d_=arb(d); Ra=arb(R)
    # psi_main is the integral over [-R,R], so this is twice the positive
    # half-line majorant on [R,infinity).
    psi_tail = 2*d_**(-2*k)*(((1+Ra).log()+1+arb.pi()/2)*Ra**(1-2*k)/(2*k-1)+Ra**(1-2*k)/(2*k-1)**2)
    PSI = -(arb(1)/F["GG"])*(1/(2*arb.pi()))*(psi_main + arb(0).union(psi_tail).union(-psi_tail))
    # PHIINT via the CONTOUR-SHIFT form (elementary + digamma; derived in
    # RIGOR_GAPS.md and numerically cross-checked against the direct zeta form).
    # phi_K(s)=(2pi/sqrt|D|)zeta_K(s-1)/((s-1)zeta_K(s)):
    #   PHIINT = (1/2pi) int h/(1+r^2) dr  +  prime  -  (1/2pi) int h(t-i) BR(t) dt,
    #   BR(t) = 1/(2+it)+1/(1+it)+C_K+psi(2+it),  C_K = (1/2)log|D| - log(2pi).
    # prime term: sum mult*log(Np)/Np^m*g(m*log(Np)); for Q(i) the
    # norm-2 prime (1+i) survives (log2<supp), for Q(omega) none (min norm 3 > supp).
    CK = arb(F["D"]).log()/2 - (2*arb.pi()).log()
    # part_elem
    pe_main = 2*integ(lambda z: sinc2k(z,k,d)/(1+acb(z)**2), 0, R).real
    # pe_main is over [-R,R]; h/(1+r^2) is nonnegative on R.
    pe_tail = 2*d_**(-2*k)*Ra**(-1-2*k)/(2*k+1)
    part_elem = (1/(2*arb.pi()))*(pe_main + arb(0).union(pe_tail))
    # part_prime
    part_prime = prime_term(group.field, k, d)
    # part_shift
    def shift_integrand(z):
        it=acb(0,1)*z
        return sinc2k(z-acb(0,1),k,d)*(1/(2+it)+1/(1+it)+CK+(2+it).digamma())
    ps_main = 2*integ(shift_integrand, 0, R, panel=2.0).real
    ch=d_.cosh()
    ck_bound = abs(CK).union((2*arb.pi()).log()).upper()
    # ps_main is over [-R,R], hence the leading factor 2 here as well.
    ps_tail = 2*(ch/d_)**(2*k)*(((2+Ra).log()+3+ck_bound)*Ra**(1-2*k)/(2*k-1)
                              + Ra**(1-2*k)/(2*k-1)**2)
    part_shift = -(1/(2*arb.pi()))*(ps_main + arb(0).union(ps_tail).union(-ps_tail))
    PHIINT = part_elem + part_prime + part_shift
    hi=h_isig(k,d)
    B = I+NCE+CE+Ch0+PARg0+PSI+PHIINT-hi
    if verbose:
        for nm,v in [("I",I),("NCE",NCE),("CE",CE),("SCATT+PARh0",Ch0),("PARg0",PARg0),
                     ("PSI",PSI),("PHIINT",PHIINT),("h(i)",hi),("B",B)]:
            print(f"    {nm:12s} = {v}")
    terms = dict(I=I, NCE=NCE, CE=CE, Ch0=Ch0, PARg0=PARg0,
                 PSI=PSI, PHIINT=PHIINT, hi=hi, prime=part_prime)
    return Evaluation(group, k, d, R, include_elliptic, ctx.prec, B, terms)
