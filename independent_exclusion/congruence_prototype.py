"""M-p0 prototype: exclusion criterion for Gamma_0(p) (CONGRUENCE.md).
Float, conforming diagnostic (M0-analog at level p).

Reference-cell architecture: one 24-split mesh of the level-1 cell K,
element data computed once; N(p)+1 copies scattered by exact coset
combinatorics over R = Z[i]/p; interfaces between distinct copies
conforming (CR face dofs merged), self-identifications relaxed.
Two-cusp pencil:

  mu(lam) = min { Q - lam*M - beta_inf*t_inf^2 - beta_0*t_0^2 : L_s = 0 }

Criterion needs mu(lam) > 0 for all lam in (0,1).

Supported levels (set via set_level); see P1_ACTION.md:
  (2+i)  -> F_5,     norm=5,  index=6,   i |-> 3
  (3)    -> F_9,     norm=9,  index=10,  pairs (a,b), i^2=-1
  (3+2i) -> F_13,    norm=13, index=14,  i |-> 8
  (5)    -> F_5xF_5, norm=25, index=36,  CRT (2+i)(2-i), i |-> (3,2)
"""

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu, lobpcg, LinearOperator

Y = 1.25
NP = 5                       # N(n) = norm; updated by set_level
INDEX = 6                    # |P1(R)| = [Gamma:Gamma_0(n)]; set_level
LEVEL = "(2+i)"
_RING = None                 # residue ring; lazy via _ring()


def yf(x1, x2):
    return np.sqrt(1.0 - x1 ** 2 - x2 ** 2)


# ----------------------------------------------------------------------
# reference cell: 24-split mesh
# ----------------------------------------------------------------------

def build_reference(N1, N2, N3):
    """24-split mesh of K. Returns nodes X, tets, boundary triangle lists
    per pairing face, top-face triangles.  N1 must be even."""
    assert N1 % 2 == 0
    x1g = np.linspace(-0.5, 0.5, N1 + 1)
    x2g = np.linspace(0.0, 0.5, N2 + 1)
    dz2 = (x1g[1] - x1g[0]) ** 2 + (x2g[1] - x2g[0]) ** 2
    base = yf(x1g[:, None], x2g[None, :])
    cl = np.empty((N1, N2))
    for i in range(N1):
        for j in range(N2):
            yfm = min(base[i, j], base[i + 1, j], base[i, j + 1],
                      base[i + 1, j + 1])
            cl[i, j] = 0.125 * dz2 / yfm ** 3 * (1 + 1e-9)
    lift = np.zeros((N1 + 1, N2 + 1))
    for i in range(N1 + 1):
        for j in range(N2 + 1):
            lift[i, j] = cl[max(i - 1, 0):min(i + 1, N1),
                            max(j - 1, 0):min(j + 1, N2)].max()

    nodes = {}
    X = []

    def nid(key, xyz):
        if key not in nodes:
            nodes[key] = len(X)
            X.append(xyz)
        return nodes[key]

    def grid_node(i, j, k):
        tt = k / N3
        yb = base[i, j] + lift[i, j]
        return nid(("g", i, j, k),
                   (x1g[i], x2g[j], yb * (1 - tt) + Y * tt))

    def face_center(corners, key, floor_cell=None):
        P = np.array([X[c] for c in corners])
        c = P.mean(0)
        if floor_cell is not None:           # floor: place on lifted sphere
            ci, cj = floor_cell
            c = np.array([c[0], c[1], yf(c[0], c[1]) + cl[ci, cj]])
        return nid(key, tuple(c))

    tets = []
    # boundary triangles per pairing face, with reference-coordinate keys
    btri = {"x1m": [], "x1p": [], "x2m": [], "x2p": [], "floor": [],
            "top": []}

    for i in range(N1):
        for j in range(N2):
            for k in range(N3):
                c = {}
                for a in (0, 1):
                    for b in (0, 1):
                        for d in (0, 1):
                            c[(a, b, d)] = grid_node(i + a, j + b, k + d)
                ctr = nid(("c", i, j, k),
                          tuple(np.mean([X[v] for v in c.values()], 0)))
                # six faces: (corner quad in cyclic order, tag or None)
                faces = [
                    ([(0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)],
                     "x1m" if i == 0 else None),
                    ([(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)],
                     "x1p" if i == N1 - 1 else None),
                    ([(0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)],
                     "x2m" if j == 0 else None),
                    ([(0, 1, 0), (1, 1, 0), (1, 1, 1), (0, 1, 1)],
                     "x2p" if j == N2 - 1 else None),
                    ([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)],
                     "floor" if k == 0 else None),
                    ([(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)],
                     "top" if k == N3 - 1 else None),
                ]
                for fi, (quad, tag) in enumerate(faces):
                    corners = [c[q] for q in quad]
                    # canonical key: shared with the neighboring hex
                    fc = face_center(
                        corners, ("f",) + tuple(sorted(corners)),
                        floor_cell=(i, j) if tag == "floor" else None)
                    for e in range(4):
                        tri = (corners[e], corners[(e + 1) % 4], fc)
                        tets.append((tri[0], tri[1], tri[2], ctr))
                        if tag is not None:
                            btri[tag].append(tri)
    return np.array(X), np.array(tets, dtype=int), btri


def ref_geometry(X, tets):
    nt = len(tets)
    vol = np.empty(nt)
    grads = np.empty((nt, 4, 3))
    for e, tet in enumerate(tets):
        P = X[list(tet)]
        A = np.hstack([np.ones((4, 1)), P])
        det = np.linalg.det(A)
        if det < 0:
            tet = (tet[1], tet[0], tet[2], tet[3])
            tets[e] = tet
            P = X[list(tet)]
            A = np.hstack([np.ones((4, 1)), P])
            det = np.linalg.det(A)
        assert det > 1e-16
        vol[e] = det / 6.0
        grads[e] = np.linalg.inv(A)[1:4, :].T
    return vol, grads


def ref_elements(X, tets, vol, grads):
    """Per-tet CR element matrices with exact-weight quadrature (deg 2)."""
    QA, QB = (5 + 3 * np.sqrt(5)) / 20, (5 - np.sqrt(5)) / 20
    qp = np.array([[QA if i == j else QB for j in range(3)]
                   for i in range(3)] + [[QB] * 3])
    nt = len(tets)
    Se = np.empty((nt, 4, 4))
    Me = np.empty((nt, 4, 4))
    ae = np.empty((nt, 4))
    Mloc = np.full((4, 4), -1 / 20.0) + np.eye(4) * (9 / 20.0)
    for e in range(nt):
        P = X[list(tets[e])]
        gphi = -3.0 * grads[e]
        wq = wm = 0.0
        av = np.zeros(4)
        for q in qp:
            lam = np.array([1 - q.sum(), *q])
            y = lam @ P[:, 2]
            wq += 0.25 / y
            wm += 0.25 / y ** 3
            av += 0.25 / y ** 3 * (1 - 3 * lam)
        Se[e] = wq * vol[e] * (gphi @ gphi.T)
        Me[e] = wm * vol[e] * Mloc
        ae[e] = vol[e] * av
    return Se, Me, ae


# ----------------------------------------------------------------------
# residue ring R = Z[i]/n and P^1(R) coset combinatorics
# (field levels + square-free CRT products; see P1_ACTION.md)
# ----------------------------------------------------------------------

class FieldResidueRing:
    """R = Z[i]/p for a single prime p (field, or F_9 for p=(3)).

    Elements:
      prime field (norm in {5,13}): int in 0..norm-1
      F_9 (norm=9): tuple (a, b) for a + b*i with a,b in Z/3
    """

    def __init__(self, kind, norm, mod=None, i_img=None, label=None):
        self.kind = kind          # 'prime' | 'f9'
        self.norm = int(norm)
        self.NP = self.norm       # alias used historically
        self.index = self.norm + 1
        self.level = label or f"N={norm}"
        self.is_product = False
        if kind == "prime":
            self.mod, self.i_img = int(mod), int(i_img)
        elif kind == "f9":
            self.mod = 3
        else:
            raise ValueError(kind)
        self.zero = self.embed(0)
        self.one = self.embed(1)
        self.i = self.embed_i()
        assert self.mul(self.i, self.i) == self.embed(-1), "i^2 != -1"

    def embed(self, n):
        n = int(n)
        if self.kind == "prime":
            return n % self.mod
        return (n % 3, 0)

    def embed_i(self):
        if self.kind == "prime":
            return self.i_img
        return (0, 1)

    def elements(self):
        if self.kind == "prime":
            return list(range(self.mod))
        return [(a, b) for a in range(3) for b in range(3)]

    def is_zero(self, a):
        return a == self.zero

    def is_unit(self, a):
        return not self.is_zero(a)

    def add(self, a, b):
        if self.kind == "prime":
            return (a + b) % self.mod
        return ((a[0] + b[0]) % 3, (a[1] + b[1]) % 3)

    def neg(self, a):
        if self.kind == "prime":
            return (-a) % self.mod
        return ((-a[0]) % 3, (-a[1]) % 3)

    def mul(self, a, b):
        if self.kind == "prime":
            return (a * b) % self.mod
        return ((a[0] * b[0] - a[1] * b[1]) % 3,
                (a[0] * b[1] + a[1] * b[0]) % 3)

    def inv(self, a):
        assert self.is_unit(a), "division by non-unit in residue ring"
        if self.kind == "prime":
            return pow(int(a), self.mod - 2, self.mod)
        nrm = (a[0] * a[0] + a[1] * a[1]) % 3
        assert nrm != 0
        ninv = 1 if nrm == 1 else 2
        conj = (a[0], (-a[1]) % 3)
        return ((conj[0] * ninv) % 3, (conj[1] * ninv) % 3)

    def mat(self, entries):
        out = []
        for row in entries:
            r = []
            for e in row:
                if e == "i":
                    r.append(self.i)
                elif e == "-i":
                    r.append(self.neg(self.i))
                else:
                    r.append(self.embed(e))
            out.append(r)
        return out

    def mat_mul(self, A, B):
        C = [[self.zero, self.zero], [self.zero, self.zero]]
        for i in range(2):
            for j in range(2):
                C[i][j] = self.add(self.mul(A[i][0], B[0][j]),
                                   self.mul(A[i][1], B[1][j]))
        return C

    def mat_eq(self, A, B):
        return all(A[i][j] == B[i][j] for i in range(2) for j in range(2))

    def mat_pm_I(self, A):
        I = self.mat([[1, 0], [0, 1]])
        mI = self.mat([[-1, 0], [0, -1]])
        return self.mat_eq(A, I) or self.mat_eq(A, mI)

    def p1_normalize(self, c, d):
        # field: unimodular iff (c,d) != (0,0)
        assert not (self.is_zero(c) and self.is_zero(d))
        if self.is_unit(c):
            invc = self.inv(c)
            return (self.one, self.mul(d, invc))
        # c = 0 ⇒ d unit
        return (self.zero, self.one)

    def p1_points(self):
        return ([(self.zero, self.one)]
                + [(self.one, d) for d in self.elements()])

    def act(self, pt, M):
        c, d = pt
        return self.p1_normalize(
            self.add(self.mul(c, M[0][0]), self.mul(d, M[1][0])),
            self.add(self.mul(c, M[0][1]), self.mul(d, M[1][1])))

    def local_cusp_is_inf(self, pt):
        return self.is_zero(pt[0])

    def generator_matrices(self):
        T1 = self.mat([[1, 1], [0, 1]])
        R = self.mat([["i", 0], [0, "-i"]])
        Ti = self.mat([[1, "i"], [0, 1]])
        TiR = self.mat_mul(Ti, R)
        S = self.mat([[0, -1], [1, 0]])
        return {"T1": T1, "R": R, "TiR": TiR, "S": S}

    def inverses(self, mats):
        inv = {"T1": self.mat([[1, -1], [0, 1]]),
               "R": mats["R"], "TiR": mats["TiR"], "S": mats["S"]}
        for name in ("R", "TiR", "S"):
            assert self.mat_pm_I(self.mat_mul(mats[name], mats[name])), \
                f"{name}^2 not ±I in PSL over {self.level}"
        assert self.mat_eq(self.mat_mul(mats["T1"], inv["T1"]),
                           self.mat([[1, 0], [0, 1]]))
        return inv


class ProductResidueRing:
    """R ≅ ∏ R_j via CRT for square-free n = p1...pk.

    P¹(R) ≅ ∏ P¹(R_j): points are tuples of local field points.
    Matrices have entries that are tuples of local ring elements.
    Right action is componentwise (see P1_ACTION.md).
    """

    def __init__(self, factors, level):
        assert len(factors) >= 2
        self.factors = list(factors)
        self.level = level
        self.is_product = True
        self.norm = 1
        self.index = 1
        for F in self.factors:
            self.norm *= F.norm
            self.index *= F.index
        self.NP = self.norm
        # ring 0/1/i as tuples (for mat())
        self.zero = tuple(F.zero for F in self.factors)
        self.one = tuple(F.one for F in self.factors)
        self.i = tuple(F.i for F in self.factors)
        for F in self.factors:
            assert F.mul(F.i, F.i) == F.embed(-1)

    def embed(self, n):
        return tuple(F.embed(n) for F in self.factors)

    def add(self, a, b):
        return tuple(F.add(a[j], b[j]) for j, F in enumerate(self.factors))

    def neg(self, a):
        return tuple(F.neg(a[j]) for j, F in enumerate(self.factors))

    def mul(self, a, b):
        return tuple(F.mul(a[j], b[j]) for j, F in enumerate(self.factors))

    def mat(self, entries):
        # each local ring builds the same integer/'i' matrix, then zip
        local = [F.mat(entries) for F in self.factors]
        out = [[None, None], [None, None]]
        for i in range(2):
            for j in range(2):
                out[i][j] = tuple(local[k][i][j] for k in range(len(local)))
        return out

    def mat_mul(self, A, B):
        C = [[self.zero, self.zero], [self.zero, self.zero]]
        for i in range(2):
            for j in range(2):
                C[i][j] = self.add(self.mul(A[i][0], B[0][j]),
                                   self.mul(A[i][1], B[1][j]))
        return C

    def mat_eq(self, A, B):
        return all(A[i][j] == B[i][j] for i in range(2) for j in range(2))

    def mat_pm_I(self, A):
        I = self.mat([[1, 0], [0, 1]])
        mI = self.mat([[-1, 0], [0, -1]])
        return self.mat_eq(A, I) or self.mat_eq(A, mI)

    def _local_mat(self, M, j):
        return [[M[a][b][j] for b in range(2)] for a in range(2)]

    def p1_points(self):
        # cartesian product of local P¹
        from itertools import product
        locals_pts = [F.p1_points() for F in self.factors]
        return [tuple(p) for p in product(*locals_pts)]

    def act(self, pt, M):
        # pt = (pt_0, ..., pt_{k-1}); componentwise right action
        return tuple(
            F.act(pt[j], self._local_mat(M, j))
            for j, F in enumerate(self.factors))

    def local_cusp_pattern(self, pt):
        # tuple of 0/1: local infinity?
        return tuple(1 if F.local_cusp_is_inf(pt[j]) else 0
                     for j, F in enumerate(self.factors))

    def generator_matrices(self):
        T1 = self.mat([[1, 1], [0, 1]])
        R = self.mat([["i", 0], [0, "-i"]])
        Ti = self.mat([[1, "i"], [0, 1]])
        TiR = self.mat_mul(Ti, R)
        S = self.mat([[0, -1], [1, 0]])
        return {"T1": T1, "R": R, "TiR": TiR, "S": S}

    def inverses(self, mats):
        inv = {"T1": self.mat([[1, -1], [0, 1]]),
               "R": mats["R"], "TiR": mats["TiR"], "S": mats["S"]}
        for name in ("R", "TiR", "S"):
            assert self.mat_pm_I(self.mat_mul(mats[name], mats[name])), \
                f"{name}^2 not ±I in PSL over {self.level}"
        assert self.mat_eq(self.mat_mul(mats["T1"], inv["T1"]),
                           self.mat([[1, 0], [0, 1]]))
        return inv


def make_ring(level):
    """Build residue ring for a supported level string / norm."""
    if level in (5, "(2+i)", "2+i"):
        return FieldResidueRing("prime", 5, mod=5, i_img=3, label="(2+i)")
    if level in (9, "(3)", "3", 3):
        return FieldResidueRing("f9", 9, label="(3)")
    if level in (13, "(3+2i)", "3+2i"):
        return FieldResidueRing("prime", 13, mod=13, i_img=8, label="(3+2i)")
    if level in (25, "(5)", "5", 5):
        # n = (5) = (2+i)(2-i); R ≅ F5 × F5, i ↦ 3 and i ↦ 2
        f1 = FieldResidueRing("prime", 5, mod=5, i_img=3, label="(2+i)")
        f2 = FieldResidueRing("prime", 5, mod=5, i_img=2, label="(2-i)")
        return ProductResidueRing([f1, f2], level="(5)")
    raise ValueError(
        f"unsupported level {level!r}; "
        f"use '(2+i)', '(3)', '(3+2i)', or '(5)'")


# backward-compatible name
GaussianResidueRing = make_ring


def set_level(level):
    """Select level n; updates NP (=norm), INDEX, LEVEL, residue ring."""
    global NP, LEVEL, _RING, INDEX
    _RING = make_ring(level)
    NP = _RING.norm
    INDEX = _RING.index
    LEVEL = _RING.level
    return _RING


def _ring():
    global _RING
    if _RING is None:
        _RING = make_ring(LEVEL)
    return _RING


def build_gluing(level=None):
    """Copy gluing: for pairing delta with source face A, target face A',
    copy p's A-face glues to copy act(p, delta^{-1})'s A'-face.

    Right action on P¹(R); see P1_ACTION.md.  If level is given, calls
    set_level first.
    """
    if level is not None:
        R = set_level(level)
    else:
        R = _ring()
    mats = R.generator_matrices()
    inv = R.inverses(mats)
    pts = R.p1_points()
    NC = R.index
    assert len(pts) == NC, f"|P1|={len(pts)} != index={NC}"
    idx = {p: n for n, p in enumerate(pts)}
    glue = {}
    for name in mats:
        perm = [idx[R.act(p, inv[name])] for p in pts]
        assert sorted(perm) == list(range(NC)), f"{name} not a permutation"
        for n, p in enumerate(pts):
            assert R.act(pts[perm[n]], mats[name]) == p, \
                f"{name} round-trip failed at copy {n}"
        glue[name] = perm

    # Cusp labels: field → {0,1}; product → local-infinity pattern ids
    if not getattr(R, "is_product", False):
        cusp_class = [0 if R.local_cusp_is_inf(p) else 1 for p in pts]
        assert cusp_class.count(0) == 1 and cusp_class.count(1) == R.norm, \
            f"prime cusp classes failed: {cusp_class.count(0)}, " \
            f"{cusp_class.count(1)} (want 1, {R.norm})"
    else:
        patterns = [R.local_cusp_pattern(p) for p in pts]
        pat_ids = {}
        cusp_class = []
        for pat in patterns:
            if pat not in pat_ids:
                pat_ids[pat] = len(pat_ids)
            cusp_class.append(pat_ids[pat])
        assert len(cusp_class) == NC
        # (5): expect patterns (1,1),(1,0),(0,1),(0,0) with sizes 1,5,5,25
        if R.level == "(5)":
            from collections import Counter
            ctr = Counter(patterns)
            assert ctr[(1, 1)] == 1 and ctr[(1, 0)] == 5 \
                and ctr[(0, 1)] == 5 and ctr[(0, 0)] == 25, ctr
    return pts, glue, cusp_class


def verify_all_levels():
    """Assert suite for primes + composite (5) (restores prior level)."""
    prev = LEVEL
    try:
        for lev in ("(2+i)", "(3)", "(3+2i)", "(5)"):
            pts, glue, cc = build_gluing(lev)
            R = _ring()
            from collections import Counter
            print(f"  gluing OK  level={lev}  norm={R.norm}  "
                  f"index={R.index}  |P1|={len(pts)}  "
                  f"cusp_hist={dict(Counter(cc))}")
            for name, perm in glue.items():
                assert len(set(perm)) == R.index
    finally:
        set_level(prev)
    print("  all levels: residue/gluing asserts passed")


# ----------------------------------------------------------------------
# assembly of N(p)+1 copies
# ----------------------------------------------------------------------

def tri_key(X, tri, mapping=None):
    pts = []
    for n in tri:
        p = X[n]
        if mapping is not None:
            p = mapping(p)
        pts.append((round(p[0], 12), round(p[1], 12), round(p[2], 12)))
    return tuple(sorted(pts))


PAIRINGS = [
    # (delta, source tag, target tag, reference map on coordinates)
    ("T1", "x1m", "x1p", lambda p: (p[0] + 1.0, p[1], p[2])),
    ("R", "x2m", "x2m", lambda p: (-p[0], p[1], p[2])),
    ("TiR", "x2p", "x2p", lambda p: (-p[0], p[1], p[2])),
    ("S", "floor", "floor", lambda p: (-p[0], p[1], p[2])),
]


def assemble_level_p(N1, N2, N3, verbose=True):
    X, tets, btri = build_reference(N1, N2, N3)
    vol, grads = ref_geometry(X, tets)
    Se, Me, ae = ref_elements(X, tets, vol, grads)
    nt = len(tets)

    # reference CR face table
    fid = {}
    tf = np.empty((nt, 4), int)
    for e in range(nt):
        for a in range(4):
            key = tuple(sorted(np.delete(tets[e], a)))
            if key not in fid:
                fid[key] = len(fid)
            tf[e, a] = fid[key]
    nfr = len(fid)

    # boundary triangle -> reference CR dof
    def dof_of(tri):
        return fid[tuple(sorted(tri))]

    # per-pairing bijection on reference boundary dofs
    pair_maps = {}
    for name, src, dst, mp in PAIRINGS:
        dst_lookup = {tri_key(X, t): dof_of(t) for t in btri[dst]}
        m = {}
        for t in btri[src]:
            key = tri_key(X, t, mapping=mp)
            m[dof_of(t)] = dst_lookup[key]
        pair_maps[name] = m

    # top-face vector and areas
    top_pairs = []
    for t in btri["top"]:
        P = X[list(t)]
        area = 0.5 * np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))
        top_pairs.append((dof_of(t), area))

    # index many copies + union-find on (copy, ref dof)
    pts, glue, cusp_class = build_gluing()
    NC = INDEX
    parent = list(range(NC * nfr))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        parent[find(x)] = find(y)

    for name, src, dst, mp in PAIRINGS:
        perm = glue[name]
        for c in range(NC):
            j = perm[c]
            if j == c:
                continue                     # self-identification: relax
            for d_src, d_dst in pair_maps[name].items():
                union(c * nfr + d_src, j * nfr + d_dst)

    gid = {}
    def g(x):
        r = find(x)
        if r not in gid:
            gid[r] = len(gid)
        return gid[r]

    gmap = np.array([g(x) for x in range(NC * nfr)])
    ng = len(gid)

    # sparse assembly (COO triplets)
    rows, cols, qv, mv = [], [], [], []
    avec = np.zeros(ng)
    t_inf = np.zeros(ng)
    t_0 = np.zeros(ng)
    for c in range(NC):
        gd = gmap[c * nfr:(c + 1) * nfr]
        ix = gd[tf]                                # (nt, 4)
        rows.append(np.repeat(ix, 4, axis=1).ravel())
        cols.append(np.tile(ix, (1, 4)).ravel())
        qv.append(Se.ravel())
        mv.append(Me.ravel())
        np.add.at(avec, ix.ravel(), ae.ravel())
        tv = t_inf if cusp_class[c] == 0 else t_0
        for d, area in top_pairs:
            tv[gd[d]] += area
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    Q = coo_matrix((np.concatenate(qv), (rows, cols)), (ng, ng)).tocsr()
    M = coo_matrix((np.concatenate(mv), (rows, cols)), (ng, ng)).tocsr()

    if verbose:
        one = np.ones(ng)
        volw = one @ (M @ one)
        # rough per-copy weight volume at 6x3x3 ~0.145; mesh-dependent
        print(f"  level={LEVEL} N(p)={NP} copies={NC} tets/copy={nt} "
              f"global CR dofs={ng} (merged {NC * nfr - ng})")
        print(f"  check 1'M1 = {volw:.6f}  ({NC}*vol_w(K_ref))")
        print(f"  check |Q@1| = {np.abs(Q @ one).max():.2e}   "
              f"t_inf(1) = {t_inf @ one:.6f}   "
              f"t_0(1) = {t_0 @ one:.6f}")
        # two-cusp prime: t_∞=1/2, t_0=N/2 (assembly still binary-buckets
        # non-∞ classes into t_0 — multi-cusp widths not modeled yet)
        if not getattr(_ring(), "is_product", False):
            assert abs(t_inf @ one - 0.5) < 1e-9
            assert abs(t_0 @ one - NP / 2) < 1e-9
        from scipy.sparse.csgraph import connected_components
        ncomp, _ = connected_components(M != 0, directed=False)
        print(f"  check connected components = {ncomp} (want 1)")
        assert ncomp == 1
    return Q, M, avec, t_inf, t_0


def constrained_min_eig_sparse(Ash, M, lowrank, ell, k=2, tol=1e-8):
    """Smallest k eigs of (Ash + sum c_i v_i v_i', M) s.t. ell'x = 0.
    Ash sparse symmetric; lowrank = [(c_i, v_i)]. LOBPCG with the
    constraint imposed via the M-orthogonal complement of M^{-1} ell."""
    ng = Ash.shape[0]

    def matvec(x):
        x = np.asarray(x).ravel()
        y = Ash @ x
        for c, v in lowrank:
            y = y + (c * (v @ x)) * v
        return y

    A = LinearOperator((ng, ng), matvec=matvec, dtype=float)
    prec = splu((Ash + 2.0 * M).tocsc())     # SPD-ish sparse preconditioner
    Mprec = LinearOperator((ng, ng), matvec=prec.solve, dtype=float)
    Yc = splu(M.tocsc()).solve(ell).reshape(-1, 1)
    rng = np.random.default_rng(0)
    X = rng.standard_normal((ng, max(k, 4)))
    # tol governs the worst vector in the block; the k we report sit
    # 2-3 orders below it (verified against dense on the coarse mesh)
    vals, _ = lobpcg(A, X, B=M, M=Mprec, Y=Yc, largest=False,
                     tol=1e-6, maxiter=300)
    return np.sort(vals)[:k]


def constrained_min_eig(A, M, ell, k=2):
    n = A.shape[0]
    v = ell.copy()
    v[0] += np.sign(ell[0] if ell[0] != 0 else 1.0) * np.linalg.norm(ell)
    beta = 2.0 / (v @ v)

    def hah(B):
        w = B @ v
        return (B - beta * np.outer(v, w) - beta * np.outer(w, v)
                + beta ** 2 * (v @ w) * np.outer(v, v))[1:, 1:]

    return eigh(hah(A), hah(M), eigvals_only=True,
                subset_by_index=[0, k - 1])


def run(N1, N2, N3, lams=(0.05, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 0.999),
        force_sparse=False, level=None):
    if level is not None:
        set_level(level)
    NC = INDEX
    print(f"--- Gamma_0{LEVEL}: {NC} copies of reference {N1}x{N2}x{N3} "
          f"(24-split), Y={Y}, N(n)={NP}, index={INDEX} ---")
    Q, M, avec, t_inf, t_0 = assemble_level_p(N1, N2, N3)
    ng = Q.shape[0]
    dense = ng <= 6000 and not force_sparse
    if dense:
        Qd, Md = Q.toarray(), M.toarray()
    print(f"  solver: {'dense eigh' if dense else 'sparse LOBPCG'}")
    print(f"  {'lambda':>8} {'s':>7} {'mu1':>10} {'mu2':>10}")
    mus = []
    for lam in lams:
        s = np.sqrt(1 - lam)
        b_inf = 2 * (1 - s) / Y ** 2
        b_0 = 2 * (1 - s) / (NP * Y ** 2)
        ell = avec + (t_inf + t_0) / ((1 + s) * Y ** 2)
        if dense:
            A = (Qd - lam * Md - b_inf * np.outer(t_inf, t_inf)
                 - b_0 * np.outer(t_0, t_0))
            m = constrained_min_eig(A, Md, ell)
        else:
            m = constrained_min_eig_sparse(
                (Q - lam * M).tocsr(), M,
                [(-b_inf, t_inf), (-b_0, t_0)], ell)
        mus.append(m[0])
        print(f"  {lam:8.4f} {s:7.4f} {m[0]:10.5f} {m[1]:10.5f}")
    print(f"  margin: min mu = {min(mus):.5f}")
    return min(mus)


if __name__ == "__main__":
    import sys
    # usage:
    #   python congruence_prototype.py [mode] [level]
    #   mode: all | validate | medium | fine | gluing | margin
    #   level: (2+i) | (3) | (3+2i)   default (2+i)
    args = sys.argv[1:]
    mode = args[0] if args else "all"
    level = args[1] if len(args) > 1 else "(2+i)"
    if mode == "gluing":
        verify_all_levels()
    else:
        set_level(level)
        if mode in ("all", "validate"):
            run(4, 2, 2, lams=(0.05, 0.999))
            run(4, 2, 2, lams=(0.05, 0.999), force_sparse=True)
            print()
        if mode in ("all", "medium", "margin"):
            run(6, 3, 3)
        if mode == "fine":
            run(8, 4, 4, lams=(0.05, 0.5, 0.999))
