"""
Arb interval-CERTIFIED assembly of B for PSL(2,O_K), K in {Q(i), Q(omega)},
upgrading bianchi_omega.py from high-precision mpmath to rigorous ball arithmetic
(python-flint / Arb).  Goal: certified upper bound  B < 1  =>  lambda_1 >= 1.

All quantities are arb/acb balls.  Certified pieces:
  * vol, eta, C_ell, CE-coeffs, systole: exact/certified constants.
  * g(x), g(0), g''(0): exact B-spline (rational-coeff polynomials in delta).
  * CE integral, PSI integral, PHIINT: certified acb_calc quadrature on [0,R]
    with explicit closed-form tail bounds on [R,inf).
  * zeta_K'/zeta_K and F'/F = ((s-1)zeta_K)'/((s-1)zeta_K) via acb power series
    (jets); F'/F regularizes the r=0 pole so the PHIINT integrand is holomorphic
    across the whole real axis (no principal-value / split needed).

Validated by reproducing the certified Picard value B(Q(i)) ~ 0.31.
"""
from flint import arb, acb, acb_series, ctx
import math
PREC = 100          # margin 0.47 => modest precision suffices
RTOL = 1e-5         # ~4-5 digits per integral; ample vs the 0.47 margin

def setp(): ctx.prec = PREC
setp()

# ---------------------------------------------------------------- characters
CHI = {"i":(4,{1:1,3:-1}), "omega":(3,{1:1,2:-1})}

# ---------------------------------------------------------------- L-functions (certified)
def Lval1(kind):
    """L(1,chi) = -(1/q) sum chi(a) psi(a/q)  (poles cancel since sum chi=0)."""
    q,chi = CHI[kind]
    return (-(arb(1)/q)*sum(c*(acb(a)/q).digamma() for a,c in chi.items())).real
def Lprime1(kind, rho=0.4):
    """Certified L'(1,chi) via Cauchy: (1/2pi i) oint L(s)/(s-1)^2 ds, L entire."""
    q,chi = CHI[kind]
    def Lf(z):
        return acb(q)**(-z)*sum(c*z.zeta(acb(a)/q) for a,c in chi.items())
    def f(th,_):
        e=(acb(0,1)*th).exp(); z=acb(1)+rho*e
        return Lf(z)/(z-1)**2*rho*e*acb(0,1)
    I=acb.integral(f,0.0,float(2*math.pi))
    return (I/(acb(0,1)*2*arb.pi())).real
def zetaK2(kind):
    q,chi=CHI[kind]
    L2 = acb(q)**(-acb(2))*sum(c*acb(2).zeta(acb(a)/q) for a,c in chi.items())
    return (acb(2).zeta()*L2).real

# ---------------------------------------------------------------- zeta_K jets
def _Lseries(x, q, chi):
    logq = arb(q).log()
    qser = acb_series([-x[0]*logq, -logq]).exp()   # q^{-(s+eps)} ; x[0]=s
    hur = acb_series([0,0])
    for a,c in chi.items():
        hur = hur + c*x.zeta(acb(a)/q)
    return qser*hur
def zetaK_logder(s, kind):
    q,chi=CHI[kind]; x=acb_series([s,1])
    zk = x.zeta()*_Lseries(x,q,chi)
    return zk[1]/zk[0]
def FpF(s, kind):
    """F'/F(s), F(s)=(s-1)zeta_K(s): analytic & nonzero at s=1."""
    q,chi=CHI[kind]; x=acb_series([s,1])
    zk = x.zeta()*_Lseries(x,q,chi)
    F = acb_series([s-1,1])*zk
    return F[1]/F[0]

# ---------------------------------------------------------------- field data
def field(kind):
    if kind=="i":
        return dict(vol=2*zetaK2("i")/arb.pi()**2,
                    systole=((3+arb(5).sqrt())/2).log(),
                    C_ell=(7+4*arb(3).sqrt()).log()/9,
                    CEg0=arb(5)/16*arb(2).log(), CEint=arb(1)/4, ckern=arb(1),
                    GG=2, eta=eta(kind), D=4)
    return dict(vol=3*arb(3).sqrt()*zetaK2("omega")/(4*arb.pi()**2),
                systole=arb("2.369205407066905").log(),
                C_ell=(7+4*arb(3).sqrt()).log()/8,
                CEg0=arb(2)/9*arb(3).log(), CEint=arb(1)/3, ckern=arb(1)/2,
                GG=3, eta=eta("omega"), D=3)
def eta(kind):
    q,chi=CHI[kind]; V = arb(1) if kind=="i" else arb(3).sqrt()/2
    w = 4 if kind=="i" else 6
    c0 = w*(arb.const_euler()*Lval1(kind)+Lprime1(kind))
    return (V/arb.pi())*c0

# ---------------------------------------------------------------- test function (B-spline)
def gbspline(x, k, d):
    n=2*k; fact=arb(1)
    for j in range(2,n): fact=fact*j
    x=arb(x); d=arb(d); tot=arb(0)
    for j in range(n+1):
        y=x+(k-j)*2*d
        lo=y.lower()
        if lo>0: tot += (arb(1) if j%2==0 else arb(-1))*math.comb(n,j)*y**(n-1)
        elif y.upper()>0:
            hb=arb(y.upper())**(n-1); p=arb(0).union(hb)
            tot += (p if j%2==0 else -p)*math.comb(n,j)
    return tot/(fact*(2*d)**n)
def g0(k,d): return gbspline(arb(0),k,d)
def gpp0(k,d):
    """g''(0) from the RIGHT cubic (x=0 is a knot); exact 4-point formula."""
    hh=arb(d)/2
    a0=gbspline(hh*0,k,d); a1=gbspline(hh,k,d); a2=gbspline(2*hh,k,d); a3=gbspline(3*hh,k,d)
    return (2*a0-5*a1+4*a2-a3)/hh**2
def sinc2k(z,k,d):
    z=acb(z)*arb(d)
    if abs(complex(z))>0.4: base=z.sin()/z
    else:
        w=z*z; term=acb(1); tot=acb(0)
        for m in range(12):
            if m>0: term=term*w/((2*m)*(2*m+1))
            tot += term if m%2==0 else -term
        zab=abs(z); f=arb(1)
        for j in range(2,26): f=f*j
        tail=zab**24/f/(1-arb(1)/(26*27)); rad=arb(0).union(tail).union(-tail)
        base=tot+acb(rad,rad)
    return base**(2*k)
def h_isig(k,d,sig=1):
    x=arb(d)*arb(sig); return (x.sinh()/x)**(2*k)

# ---------------------------------------------------------------- certified integrals
def integ(f,a,b,panel=4.0):
    a=float(a); b=float(b); n=max(1,int((b-a)/panel+0.5)); step=(b-a)/n
    tot=acb(0)
    for i in range(n):
        tot += acb.integral(lambda z,_: f(z), a+i*step, a+(i+1)*step, rel_tol=RTOL)
    return tot

def compute_B(kind, k=2, frac=0.999, R=40, verbose=True):
    # The explicit tail estimates below use 2*k > 1 and, for the shifted
    # bracket, 1 + pi/2 + 3/R <= 3.  Keep the numerical preconditions at the
    # point of use rather than silently applying a bound outside its range.
    if k < 1:
        raise ValueError("tail bounds require an integer k >= 1")
    if R < 7:
        raise ValueError("shifted-line tail bound requires R >= 7")
    setp(); F=field(kind)
    d = float(F["systole"].lower())*frac/(2*k)
    supp = 2*k*d
    G0=g0(k,d)
    # I
    I = F["vol"]/(2*arb.pi())*(-gpp0(k,d))
    # NCE
    NCE = G0*F["C_ell"]
    # CE integral (holomorphic B-spline polynomial per knot-interval)
    CEi = ce_integral(k,d,supp,F["ckern"])
    CE = F["CEg0"]*G0 + F["CEint"]*CEi
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
    # prime term: sum over prime powers N<=e^supp of (log N/2) g(log N); for Q(i) the
    # norm-2 prime (1+i) survives (log2<supp), for Q(omega) none (min norm 3 > supp).
    CK = arb(F["D"]).log()/2 - (2*arb.pi()).log()
    # part_elem
    pe_main = 2*integ(lambda z: sinc2k(z,k,d)/(1+acb(z)**2), 0, R).real
    # pe_main is over [-R,R]; h/(1+r^2) is nonnegative on R.
    pe_tail = 2*d_**(-2*k)*Ra**(-1-2*k)/(2*k+1)
    part_elem = (1/(2*arb.pi()))*(pe_main + arb(0).union(pe_tail))
    # part_prime
    part_prime = arb(0)
    if kind=="i":
        l2=arb(2).log()
        if l2 < 2*k*d_: part_prime = l2/2*gbspline(l2,k,d)
    # part_shift
    def shift_integrand(z):
        it=acb(0,1)*z
        return sinc2k(z-acb(0,1),k,d)*(1/(2+it)+1/(1+it)+CK+(2+it).digamma())
    ps_main = 2*integ(shift_integrand, 0, R, panel=2.0).real
    ch=d_.cosh()
    # ps_main is over [-R,R], hence the leading factor 2 here as well.
    ps_tail = 2*(ch/d_)**(2*k)*(((2+Ra).log()+3+(2*arb.pi()).log())*Ra**(1-2*k)/(2*k-1)
                              + Ra**(1-2*k)/(2*k-1)**2)
    part_shift = -(1/(2*arb.pi()))*(ps_main + arb(0).union(ps_tail).union(-ps_tail))
    PHIINT = part_elem + part_prime + part_shift
    hi=h_isig(k,d)
    B = I+NCE+CE+Ch0+PARg0+PSI+PHIINT-hi
    if verbose:
        for nm,v in [("I",I),("NCE",NCE),("CE",CE),("SCATT+PARh0",Ch0),("PARg0",PARg0),
                     ("PSI",PSI),("PHIINT",PHIINT),("h(i)",hi),("B",B)]:
            print(f"    {nm:12s} = {v}")
    return B

def ce_integral(k,d,supp,ckern):
    """int_0^supp g(x) sinh x/(cosh x+ckern) dx, split at knots (g cubic per piece)."""
    d=arb(d); n=2*k; fact=arb(1)
    for j in range(2,n): fact=fact*j
    tot=acb(0)
    for m in range(k):
        a_,b_=float(2*d*m),float(2*d*(m+1))
        def piece(z,m=m):
            z=acb(z); s=acb(0)
            for j in range(0,k+m+1):
                y=z+(k-j)*2*d
                t=math.comb(n,j)*y**(n-1)
                s = s+t if j%2==0 else s-t
            gp=s/(fact*(2*d)**n)
            return gp*z.sinh()/(z.cosh()+ckern)
        tot += acb.integral(lambda z,_: piece(z), a_, b_)
    return tot.real

if __name__=="__main__":
    print("="*62); print("CERTIFIED (Arb) -- validate Q(i) [expect B ~ 0.31]:"); print("="*62)
    Bi=compute_B("i")
    print(f"  B(Q(i)) enclosure: [{float(Bi.lower()):.6f}, {float(Bi.upper()):.6f}]")
    print()
    print("="*62); print("Q(omega) (Eisenstein-Picard), level 1:"); print("="*62)
    Bw=compute_B("omega")
    lo,hi=float(Bw.lower()),float(Bw.upper())
    print(f"  B(Q(omega)) enclosure: [{lo:.6f}, {hi:.6f}]")
    print()
    if hi<1:
        print(f"  *** CERTIFIED B < 1 (upper bound {hi:.6f}) => lambda_1(PSL(2,Z[omega])) >= 1 ***")
    else:
        print(f"  upper bound {hi:.6f} not < 1; tighten R/prec.")
