"""
Direct computation of the CUSPIDAL ELLIPTIC (CE) data for Bianchi groups,
validated on Q(i) then applied to Q(omega).

CE term (Friedman arXiv:math/0612807 Lemma 4.3.2, after the log A divergence
cancels against the parabolic term):
  CE = sum_i [ 1/(|C(g_i)| |1-eps_i^2|^2) ]
             * [ 2 g(0) log|c_i|  +  Integral_i ],
  Integral_i = int_0^inf g(x) sinh x / (cosh x - 1 + |1-eps_i^2|^2/2) dx.

Cuspidal elliptic elements fix the cusp at infinity (up to Gamma):
  g = [[u, b],[0, 1/u]],  u a root of unity in O_K^*, u^2 != 1,
acting as z -> u^2 z + u b, with the OTHER (finite) fixed point
  z* = u b / (1 - u^2),  a cusp of Gamma.
 - The Gamma-conjugacy class is labelled by z* mod O_K (translations).
 - |c_i| = |denominator of z* in lowest terms|  (c_i = lower-left of gamma_i
   sending infinity -> z*).
 - |1-eps_i^2|^2:  order 2 (Q(i)) -> |1-(-1)|^2 = 4;  order 3 (Q(omega)) -> 3.
 - |C(g_i)| = #{ p in O_K^* : b (p - 1/p)/(u - 1/u) in O_K }  (centralizer in SL).

Validation target (Q(i), Friedman/EGM):
  4 classes, |C|=4, |c_i| = {1, 2, 2, sqrt2}  =>
  g(0)-coeff = sum 2 log|c_i|/(4*|C_i|) = (5/16) log2,
  integral-coeff = sum 1/(4*|C_i|) = 1/4  (kernel sinh/(cosh+1)=tanh(x/2)).
"""
import math, cmath, itertools
from collections import Counter
from fractions import Fraction

# ---------------- O_K arithmetic with Euclidean gcd ----------------
class OK:
    def __init__(self, kind):
        self.kind = kind
        self.g = 1j if kind=="i" else cmath.exp(2j*math.pi/3)
    def mul(self,x,y):
        a,b=x;c,d=y
        if self.kind=="i": return (a*c-b*d, a*d+b*c)
        return (a*c-b*d, a*d+b*c-b*d)          # omega^2=-omega-1
    def add(self,x,y): return (x[0]+y[0],x[1]+y[1])
    def sub(self,x,y): return (x[0]-y[0],x[1]-y[1])
    def neg(self,x): return (-x[0],-x[1])
    def emb(self,x): return x[0]+x[1]*self.g
    def norm(self,x):
        a,b=x; return a*a+b*b if self.kind=="i" else a*a-a*b+b*b
    def conjugate(self,x):
        """Exact nontrivial Galois conjugation of an O_K element."""
        a,b=x
        return (a,-b) if self.kind=="i" else (a-b,-b)
    def divides(self,x,y):
        """Exact test x in y O_K, avoiding floating-point Euclidean division."""
        n=self.norm(y)
        if n == 0: raise ZeroDivisionError("division by zero in O_K")
        a,b=self.mul(x,self.conjugate(y))
        return a % n == 0 and b % n == 0
    def units(self):
        if self.kind=="i": return [(1,0),(-1,0),(0,1),(0,-1)]
        return [(1,0),(-1,0),(0,1),(0,-1),(-1,-1),(1,1)]  # +-1,+-omega,+-omega^2
        # omega=(0,1); omega^2 = -omega-1 = (-1,-1); -omega^2=(1,1)
    def round(self, zc):
        """nearest O_K element to complex zc."""
        if self.kind=="i":
            return (round(zc.real), round(zc.imag))
        n = round(2*zc.imag/math.sqrt(3))
        m = round(zc.real + n/2.0)
        return (m,n)
    def divmod(self,a,b):
        """a = q*b + r, N(r)<N(b). Euclidean."""
        q = self.round(self.emb(a)/self.emb(b))
        r = self.sub(a, self.mul(q,b))
        return q,r
    def gcd(self,a,b):
        while b!=(0,0):
            _,r=self.divmod(a,b); a,b=b,r
        return a
    def is_unit(self,x): return self.norm(x)==1
    def exact_div(self,a,b):
        q,r=self.divmod(a,b)
        assert r==(0,0), f"{a} not divisible by {b}"
        return q

def zstar_key_and_c(R, num, den):
    """Class label = z* reduced to the O_K fundamental domain (complex, so that
       i/2 ~ 3i/2 etc. merge correctly). |c| = |denominator in lowest terms|."""
    g = R.gcd(num, den)
    numr = R.exact_div(num, g); denr = R.exact_div(den, g)   # lowest terms
    c = abs(R.emb(denr))
    z = R.emb(numr) / R.emb(denr)
    z = z - R.emb(R.round(z))          # reduce mod O_K into fundamental domain
    # fold boundary: among O_K-translates, pick the canonical one (min modulus,
    # then max real, then max imag) so opposite boundary edges (i/2 ~ -i/2) merge.
    best=None
    for m in range(-1,2):
        for n in range(-1,2):
            w = z + R.emb((m,n))
            cand = (round(abs(w),6), -round(w.real,6), -round(w.imag,6))
            if best is None or cand < best: best=cand
    key = (best[0], best[1], best[2])
    return key, c

def centralizer_order(R, u, b):
    """|C(g)| = #{p in O_K^*: b(p-1/p)/(u-1/u) in O_K}. u a unit."""
    # u - 1/u : 1/u = conjugate-unit; for units, inverse is a unit.
    uinv = unit_inverse(R,u)
    umm = R.sub(u,uinv)          # u - 1/u
    cnt=0
    for p in R.units():
        pinv=unit_inverse(R,p)
        pmm=R.sub(p,pinv)
        num=R.mul(b,pmm)
        # need num/umm in O_K
        if R.divides(num,umm): cnt+=1
    return cnt

def unit_inverse(R,u):
    for v in R.units():
        if R.mul(u,v)==(1,0): return v
    raise ValueError("not a unit")

def residue_representatives(R, q, B=1):
    """Exact representatives in a small box for O_K/(q).

    In the only call below N(q)=3, so the assertion that this yields exactly
    three classes is independently checked against the algebraic ideal norm.
    """
    representatives=[]
    for a in range(-B,B+1):
        for b in range(-B,B+1):
            x=(a,b)
            if not any(R.divides(R.sub(x,y),q) for y in representatives):
                representatives.append(x)
    assert len(representatives)==R.norm(q), (q, representatives, R.norm(q))
    return representatives

def exact_eisenstein_CE_check(verbose=True):
    """Exact finite algebra behind the six Eisenstein CE candidates.

    Every CE class has a representative g_{u,b} with u in {omega,omega^2}
    and b modulo (u-u^{-1}); the latter ideal has norm 3.  The global
    exhaustiveness conclusion then follows from Friedman's divergence identity
    (documented in RIGOR_GAPS.md), while all arithmetic checked here is integer
    arithmetic in Z[omega].
    """
    R=OK("omega")
    us=[(0,1),(-1,-1)]             # omega, omega^2 modulo +/- I
    weights=Fraction(0)
    c2=Counter()
    candidate_count=0
    for u in us:
        uinv=unit_inverse(R,u)
        q=R.sub(u,uinv)
        assert R.norm(q)==3
        reps=residue_representatives(R,q)
        assert len(reps)==3
        den=R.sub((1,0),R.mul(u,u))
        assert R.norm(den)==3       # den is an associate of q
        for b in reps:
            C=centralizer_order(R,u,b)
            assert C==6, (u,b,C)
            # The norm-3 denominator is prime: it cancels exactly for the
            # zero residue, and otherwise the reduced denominator has norm 3.
            c2[1 if R.divides(R.mul(u,b),den) else 3] += 1
            weights += Fraction(1,3*C)
            candidate_count += 1
    assert candidate_count==6
    assert c2==Counter({1:2,3:4}), c2
    assert weights==Fraction(1,3), weights
    assert 2*weights+Fraction(1,3)==1  # Friedman's identity, exactly
    if verbose:
        print("EXACT EISENSTEIN CE CHECK")
        print("  2 rotation multipliers x 3 residue classes = 6 candidates")
        print("  all centralizer orders are 6; |c|^2 multiset = {1: 2, 3: 4}")
        print("  2*sum 1/(3|C|) + 1/3 = 1 exactly")
    return dict(candidates=candidate_count, c2=dict(c2), weight=weights)

def run_CE(kind, Bmax=6, verbose=True):
    R=OK(kind)
    if kind=="i":
        us=[(0,1)]                 # u=i (u=-i redundant); order 2
        eps2=4; m_order=2
    else:
        us=[(0,1),(-1,-1)]         # u=omega, u=omega^2 ; order 3
        eps2=3; m_order=3
    classes={}   # canonical z* -> dict(u,b,|c|,|C|)
    for u in us:
        u2=R.mul(u,u)
        one=(1,0)
        den=R.sub(one,u2)          # 1-u^2
        if den==(0,0): continue
        for br in range(-Bmax,Bmax+1):
            for bi in range(-Bmax,Bmax+1):
                b=(br,bi)
                num=R.mul(u,b)     # u*b ; z* = num/den
                zkey,cc = zstar_key_and_c(R,num,den)   # handles num=0 -> |c|=1, z*=0
                # class label also carries the rotation multiplier u^2 (g vs g^{-1}
                # are DISTINCT conjugacy classes for order 3): key=(z*, u^2).
                key = (zkey, u2)
                if key not in classes:
                    classes[key]=dict(u=u,b=b,c=cc,C=centralizer_order(R,u,b))
    # assemble
    g0coeff=0.0; intcoeff=0.0
    reps=[]
    for key,d in classes.items():
        w=1.0/(d["C"]*eps2)
        g0coeff += 2*math.log(d["c"])*w
        intcoeff += w
        reps.append((d["c"],d["C"]))
    if verbose:
        print(f"K=Q({kind}) cuspidal elliptic (order {m_order}, |1-eps^2|^2={eps2}):")
        print(f"  #classes = {len(classes)}")
        from collections import Counter
        print(f"  (|c|,|C|) multiset: {Counter([(round(c,4),C) for c,C in reps])}")
        print(f"  sum 1/(|C||1-eps^2|^2) = {intcoeff:.6f}  (group-id check: expect (k-l/[G:G'])... )")
        print(f"  g(0)-coefficient  sum 2log|c|/(|C||1-eps^2|^2) = {g0coeff:.8f}")
        # group relation: 2*sum 1/(|C||1-eps^2|^2) + l/[G:G'] = k  (=1 for trivial rep)
        GG = 2 if kind=="i" else 3
        lhs = 2*intcoeff + 1.0/GG
        print(f"  group identity 2*sum1/(|C||1-eps^2|^2)+1/[G:G']={lhs:.6f} (expect 1.0)")
    return dict(g0coeff=g0coeff, intcoeff=intcoeff, eps2=eps2, nclasses=len(classes))

if __name__=="__main__":
    print("="*64); print("VALIDATE on Q(i): expect g(0)-coeff=(5/16)log2=",
          f"{5/16*math.log(2):.8f}, integral-coeff=1/4"); print("="*64)
    ri=run_CE("i")
    ok = abs(ri["g0coeff"]-5/16*math.log(2))<1e-6 and abs(ri["intcoeff"]-0.25)<1e-9
    print(f"  >> {'VALIDATED' if ok else 'MISMATCH'}")
    assert ok, "Picard validation must pass before applying the Eisenstein check"
    print()
    print("="*64); print("APPLY to Q(omega):"); print("="*64)
    rw=run_CE("omega")
    print()
    exact_eisenstein_CE_check()
