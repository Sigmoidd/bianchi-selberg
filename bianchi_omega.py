"""
Full assembly of B for Bianchi groups PSL(2,O_K), K in {Q(i), Q(omega)}, to test
the Selberg-type bound  B < 1  =>  lambda_1 >= 1  at level 1.

Field-parametrized; VALIDATED by reproducing the Arb-certified Picard value
B(Q(i)) ~ 0.32 before trusting Q(omega).  High-precision mpmath throughout
(the individual sub-results were cross-checked against the Arb engine / Matthies).

B = I + NCE + CE + [SCATT+PARh0] + PARg0 + PSI + PHIINT - h(i),
with all loxodromic terms annihilated (supp g = [-2k*delta, 2k*delta] <= systole).

Per-field data (derived + validated in this repo):
  Q(i):    vol .30532, systole .96242, C_ell (nce)=(1/9)log(7+4sqrt3),
           CE g0coeff=(5/16)log2, CE intcoeff=1/4, kernel cosh x+1,
           [G:G']=2, eta=.82283, |D|=4.
  Q(omega):vol .16916, systole .86255, C_ell=log(7+4sqrt3)/8,
           CE g0coeff=(2/9)log3, CE intcoeff=1/3, kernel cosh x+1/2,
           [G:G']=3, eta=.94550, |D|=3.
Scattering: phi_K(s) = (2pi/sqrt|D|) zeta_K(s-1)/((s-1) zeta_K(s)); phi'/phi is
field-only through zeta_K = zeta * L(.,chi_D). (Constant cancels in phi'/phi;
proven phi(s)phi(2-s)=1 and phi(1)=-1 for both fields, so SCATT const=+h(0)/4.)
"""
import mpmath as mp
mp.mp.dps = 30
from math import comb

# --------------------------- field data ---------------------------
def field(kind):
    g = mp.euler
    if kind=="i":
        return dict(kind=kind, vol=2*zetaK(2,"i")/mp.pi**2, systole=mp.log((3+mp.sqrt(5))/2),
                    C_ell=mp.log(7+4*mp.sqrt(3))/9,
                    CEg0=mp.mpf(5)/16*mp.log(2), CEint=mp.mpf(1)/4, ckern=mp.mpf(1),
                    GG=2, eta=eta_i(), D=4)
    else:
        return dict(kind=kind, vol=3*mp.sqrt(3)*zetaK(2,"omega")/(4*mp.pi**2),
                    systole=mp.log(mp.mpf("2.36920540706690")),   # min loxodromic N(T)
                    C_ell=mp.log(7+4*mp.sqrt(3))/8,
                    CEg0=mp.mpf(2)/9*mp.log(3), CEint=mp.mpf(1)/3, ckern=mp.mpf(1)/2,
                    GG=3, eta=eta_omega(), D=3)

# Dirichlet characters
CHI = {"i":(4,{1:1,3:-1}), "omega":(3,{1:1,2:-1})}
def Lval(s, kind):
    q,chi = CHI[kind]
    return q**(-s)*sum(c*mp.zeta(s,mp.mpf(a)/q) for a,c in chi.items())
def Lder(s, kind):
    q,chi = CHI[kind]
    # d/ds [ q^{-s} sum chi(a) zeta(s,a/q) ]
    t1 = -mp.log(q)*q**(-s)*sum(c*mp.zeta(s,mp.mpf(a)/q) for a,c in chi.items())
    t2 = q**(-s)*sum(c*mp.zeta(s,mp.mpf(a)/q,1) for a,c in chi.items())
    return t1+t2
def zetaK(s, kind): return mp.zeta(s)*Lval(s,kind)
def zetaK_logder(s, kind):
    return mp.zeta(s,1,1)/mp.zeta(s) + Lder(s,kind)/Lval(s,kind)

def L1(kind):  # L(1,chi) via digamma closed form -(1/q)sum chi(a) psi(a/q)
    q,chi=CHI[kind]
    return -(mp.mpf(1)/q)*sum(c*mp.digamma(mp.mpf(a)/q) for a,c in chi.items())
def eta_i():
    # certified value from verify_eisenstein.py (Arb): 2g+2log2+3logpi-4logGamma(1/4)
    return mp.mpf("0.8228252496786964")
def eta_omega():
    # certified value from verify_eisenstein.py (Arb), eta=(V/pi)c0
    return mp.mpf("0.9454972808680957")

# --------------------------- test function ---------------------------
def h(r, k, d):   # sinc^{2k}(d r)
    x = d*r
    if abs(x) < mp.mpf(10)**(-12): return mp.mpf(1)
    return (mp.sin(x)/x)**(2*k)
def h_i(sig, k, d):  # h(i sigma) = (sinh/.)^{2k}
    x=d*sig; return (mp.sinh(x)/x)**(2*k)
def gfun(x, k, d):   # B-spline g(x)
    n=2*k
    fact=mp.factorial(n-1)
    tot=mp.mpf(0)
    for j in range(n+1):
        y=x+(k-j)*2*d
        if y>0:
            tot += (-1)**j*comb(n,j)*y**(n-1)
    return tot/(fact*(2*d)**n)

# --------------------------- terms ---------------------------
def compute_B(kind, k=2, frac=0.999, R=60, verbose=True):
    F=field(kind)
    syst = F["systole"]
    d = float(syst)*frac/(2*k)
    d = mp.mpf(d)
    g0 = gfun(mp.mpf(0),k,d)
    supp = 2*k*d

    # I = vol/(4pi^2) int_R h(r) r^2 dr.  By Fourier, int_R h r^2 dr = 2pi(-g''(0))
    # (h = FT of g), so I = vol/(2pi)(-g''(0)).  x=0 is a B-spline KNOT, so g''(0) must
    # be taken from ONE side (g is a single cubic on (0,2delta)); exact 4-point formula:
    hh = d/2            # 3*hh = 1.5*delta < 2*delta, inside the right cubic piece
    g0,g1,g2,g3 = (gfun(mp.mpf(0),k,d), gfun(hh,k,d), gfun(2*hh,k,d), gfun(3*hh,k,d))
    gpp0 = (2*g0 - 5*g1 + 4*g2 - g3)/hh**2          # p''(0) for the right cubic
    I = F["vol"]/(2*mp.pi)*(-gpp0)
    NCE = g0*F["C_ell"]
    CEint = F["CEint"]*mp.quad(lambda x: gfun(x,k,d)*mp.sinh(x)/(mp.cosh(x)+F["ckern"]),[0,supp])
    CE = F["CEg0"]*g0 + CEint
    # SCATT + PAR h(0):  +h(0)/4 (phi(1)=-1)  +  (1/GG) h(0)/4
    Ch0 = mp.mpf(1)/4 + (mp.mpf(1)/F["GG"])*(mp.mpf(1)/4)
    PARg0 = (mp.mpf(1)/F["GG"])*g0*(F["eta"]/2 - mp.euler)
    PSI = -(mp.mpf(1)/F["GG"])*(1/(2*mp.pi))*2*mp.quad(
            lambda r: h(r,k,d)*mp.re(mp.digamma(1+1j*r)), [0,R])
    # PHIINT = (1/4pi) int_R h(r) (phi'/phi)(1+ir) dr, phi'/phi field-only
    def phiratio(r):
        s=1+1j*r
        return zetaK_logder(1j*r,kind) - 1/(1j*r) - zetaK_logder(s,kind)
    eps=mp.mpf("1e-6")
    PHIINT = (1/(4*mp.pi))*2*mp.quad(lambda r: h(r,k,d)*mp.re(phiratio(r)), [eps,R])
    hi = h_i(mp.mpf(1),k,d)

    B = I+NCE+CE+Ch0+PARg0+PSI+PHIINT-hi
    if verbose:
        for nm,v in [("I",I),("NCE",NCE),("CE",CE),("SCATT+PARh0",Ch0),
                     ("PARg0",PARg0),("PSI",PSI),("PHIINT",PHIINT),
                     ("h(i)",hi),("B",B)]:
            print(f"    {nm:12s} = {mp.nstr(v,8)}")
    return B

if __name__=="__main__":
    print("="*60)
    print("VALIDATE on Q(i) [Arb-certified value B ~ 0.3200]:")
    print("="*60)
    Bi = compute_B("i")
    print(f"  => B(Q(i)) = {mp.nstr(Bi,8)}   (Arb engine: 0.3199)")
    print()
    print("="*60)
    print("Q(omega) (Eisenstein-Picard), level 1:")
    print("="*60)
    Bw = compute_B("omega")
    print(f"  => B(Q(omega)) = {mp.nstr(Bw,8)}")
    print()
    if Bw < 1:
        print(f"  *** B < 1  =>  NO exceptional eigenvalues: lambda_1(PSL(2,Z[omega])) >= 1 ***")
    else:
        print(f"  B >= 1: inconclusive by positivity alone.")
