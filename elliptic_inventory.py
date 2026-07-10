"""
Direct computation of the NON-CUSPIDAL ELLIPTIC term of the Selberg trace
formula for Bianchi groups PSL(2, O_K), K in {Q(i), Q(omega)}.

The term (Friedman arXiv:math/0612807 Thm 4.1.1; Balkanova et al 1712.00880):
   NCE = g(0) * C_ell,
   C_ell = sum_{ {R} nce }  logN(T0) / ( 4 |E(R)| sin^2(pi k / m(R)) ),
the sum over Gamma-conjugacy classes of non-cuspidal elliptic elements.
In the counting-function limit g(0)=k/pi, so a3-elliptic = (1/pi) C_ell.

Validation target (Matthies, via Aurich-Steiner-Then):
   Q(i):  C_ell = (1/9) log(7+4 sqrt3) = (2/9) log(2+sqrt3).

Method (exact O_K arithmetic, no floats in the group theory):
 1. Enumerate primitive non-cuspidal elliptic elements M in PSL(2,O_K)
    (trace t: t=+-1 order 3 for Q(i); t=0 order 2 for Q(omega)), bounded entries.
 2. Cluster into Gamma-conjugacy classes by conjugating with a bounded generating
    set (union-find). [Validated: must give ONE class for Q(i) order 3.]
 3. Per class rep M:
      - m(R) = order (2 or 3), sin^2(pi/m).
      - |E(R)| = # elliptic elements in the centralizer C(M) (mod +-I): count
        (u,v) in O_K^2 with u^2 + t u v + v^2 = 1 and |u + v e^{i theta}| = 1.
      - N(T0) = min |u + v e^{i theta}|^2 > 1 over the SAME (u,v) solutions
        (loxodromic units of the order O_K[M]); the fundamental unit.
 4. C_ell = sum over classes of logN(T0)/(4|E(R)| sin^2(pi/m)).
"""
import math, cmath, itertools
from fractions import Fraction

# ---------------------------------------------------------------- O_K ring
class Ring:
    """O_K = Z[g].  Q(i): g=i, g^2=-1.  Q(omega): g=omega, g^2=-g-1."""
    def __init__(self, kind):
        self.kind = kind
        self.g = 1j if kind == "i" else cmath.exp(2j*math.pi/3)
    def mul(self, x, y):
        a,b = x; c,d = y
        if self.kind == "i":
            return (a*c - b*d, a*d + b*c)
        else:  # omega^2 = -omega - 1
            return (a*c - b*d, a*d + b*c - b*d)
    def add(self, x, y): return (x[0]+y[0], x[1]+y[1])
    def sub(self, x, y): return (x[0]-y[0], x[1]-y[1])
    def neg(self, x): return (-x[0], -x[1])
    def emb(self, x): return x[0] + x[1]*self.g          # complex embedding
    def norm(self, x):
        a,b = x
        return a*a + b*b if self.kind == "i" else a*a - a*b + b*b
ZERO=(0,0); ONE=(1,0)

# ---------------------------------------------------------------- 2x2 matrices over O_K
class M2:
    def __init__(self, R): self.R = R
    def mul(self, X, Y):
        R=self.R; (A,B,C,D)=X; (E,F,G,H)=Y
        return (R.add(R.mul(A,E),R.mul(B,G)), R.add(R.mul(A,F),R.mul(B,H)),
                R.add(R.mul(C,E),R.mul(D,G)), R.add(R.mul(C,F),R.mul(D,H)))
    def det(self, X):
        R=self.R; (A,B,C,D)=X; return R.sub(R.mul(A,D),R.mul(B,C))
    def trace(self, X): return self.R.add(X[0], X[3])
    def neg(self, X): return tuple(self.R.neg(e) for e in X)
    def eq(self, X, Y): return X==Y
    def is_pmI(self, X):
        I=(ONE,ZERO,ZERO,ONE)
        return X==I or X==self.neg(I)
    def canon(self, X):
        """PSL representative: min(X, -X) lexicographically."""
        nX=self.neg(X); return min(X,nX)

# ---------------------------------------------------------------- enumeration
def okelts(B):
    return [(a,b) for a in range(-B,B+1) for b in range(-B,B+1)]

def enumerate_elliptic(R, t, B=2):
    """Primitive non-cuspidal elliptic PSL elements of trace t, entries coords in [-B,B]."""
    MM=M2(R); elts=set()
    elset=okelts(B)
    for A in elset:
        D=(t - A[0], -A[1])            # trace = t (real) => A+D=(t,0)
        for Bm in elset:
            for C in elset:
                X=(A,Bm,C,D)
                if MM.det(X)!=ONE: continue
                if MM.is_pmI(X): continue
                elts.add(MM.canon(X))
    return MM, elts

def conjugators(R, B=1):
    MM=M2(R); gs=[]
    elset=okelts(B)
    for A in elset:
        for Bm in elset:
            for C in elset:
                for D in elset:
                    X=(A,Bm,C,D)
                    if MM.det(X)==ONE:
                        gs.append(X)
    return gs

def inv_SL(R, X):
    A,B,C,D=X; return (D, R.neg(B), R.neg(C), A)   # det 1 inverse

def cluster(MM, elts, gs):
    """union-find conjugacy clustering."""
    R=MM.R; elts=list(elts); idx={e:i for i,e in enumerate(elts)}
    parent=list(range(len(elts)))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    def union(a,b):
        ra,rb=find(a),find(b)
        if ra!=rb: parent[ra]=rb
    eset=set(elts)
    for M in elts:
        for g in gs:
            gi=inv_SL(R,g)
            conj=MM.mul(MM.mul(g,M),gi)
            c=MM.canon(conj)
            if c in idx:
                union(idx[M], idx[c])
    classes={}
    for e in elts:
        classes.setdefault(find(idx[e]), []).append(e)
    return list(classes.values())

def centralizer_data(R, t, theta, B2=4):
    """|E(R)| (mod +-I) and N(T0) via (u,v) in O_K^2, u^2+t uv+v^2 = 1."""
    e_ith = cmath.exp(1j*theta)
    tw=(t,0)
    fin=set(); Nmin=None
    for u in okelts(B2):
        for v in okelts(B2):
            # det = u^2 + t u v + v^2 in O_K, must equal 1
            u2=R.mul(u,u); v2=R.mul(v,v); uv=R.mul(u,v)
            det=R.add(R.add(u2,v2), R.mul(tw,uv))
            if det!=ONE: continue
            lam = R.emb(u) + R.emb(v)*e_ith
            mod2 = abs(lam)**2
            if abs(mod2-1) < 1e-9:
                # elliptic unit; record mod +-I : (u,v) ~ (-u,-v)
                key=min((u,v),(R.neg(u),R.neg(v)))
                fin.add(key)
            elif mod2 > 1+1e-9:
                if Nmin is None or mod2 < Nmin-1e-9:
                    Nmin=mod2
    return len(fin), Nmin

def run(kind, verbose=True):
    R=Ring(kind); MM=M2(R)
    if kind=="i":
        t, m, theta = 1, 3, math.pi/3      # order-3 non-cuspidal
    else:
        t, m, theta = 0, 2, math.pi/2      # order-2 non-cuspidal
    _, elts = enumerate_elliptic(R, t, B=2)
    # also include trace -t elements? In PSL, -M has trace -t; canon merges M,-M,
    # so enumerating t and using canon covers both signs.
    gs = conjugators(R, B=1)
    classes = cluster(MM, elts, gs)
    nE, Nmin = centralizer_data(R, t, theta, B2=4)
    sin2 = math.sin(math.pi/m)**2
    nclasses=len(classes)
    logN = math.log(Nmin) if Nmin else float('nan')
    C_ell = nclasses * logN / (4*nE*sin2)
    if verbose:
        print(f"K=Q({kind}): non-cuspidal elliptic order {m}")
        print(f"  #elliptic elements enumerated (mod +-I): {len(elts)}")
        print(f"  #conjugacy classes (clustered): {nclasses}")
        print(f"  |E(R)| (elliptic centralizer, mod +-I): {nE}")
        print(f"  N(T0) = {Nmin:.6f}   (log = {logN:.6f})")
        print(f"  sin^2(pi/{m}) = {sin2:.6f}")
        print(f"  C_ell = {nclasses} * log N(T0) / (4*{nE}*{sin2:.4f}) = {C_ell:.8f}")
    return dict(m=m, nclasses=nclasses, nE=nE, Nmin=Nmin, C_ell=C_ell)

if __name__ == "__main__":
    print("="*66)
    print("VALIDATION on Q(i): target C_ell = (1/9) log(7+4 sqrt3) = "
          f"{math.log(7+4*math.sqrt(3))/9:.8f}")
    print("="*66)
    ri = run("i")
    target = math.log(7+4*math.sqrt(3))/9
    ok = abs(ri["C_ell"] - target) < 1e-6
    print(f"  >> {'VALIDATED' if ok else 'MISMATCH'} (got {ri['C_ell']:.8f}, target {target:.8f})")
    print()
    print("="*66)
    print("APPLY to Q(omega) (Eisenstein-Picard):")
    print("="*66)
    rw = run("omega")
    if ok:
        print(f"\n  => C_ell(Z[omega]) = {rw['C_ell']:.8f}  [trusting validated method]")
    else:
        print("\n  method NOT validated on Q(i); Z[omega] number withheld.")
