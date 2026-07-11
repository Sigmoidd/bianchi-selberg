"""Pencil probe on 8x3x2 for N(p)=13."""
import warnings
import numpy as np
from flint import arb
from scipy.sparse.linalg import splu, lobpcg, LinearOperator

import congruence_prototype as cp
import m3p_certify as m

cp.set_level("(3+2i)")
import sys
N1 = int(sys.argv[1]) if len(sys.argv) > 1 else 8
N2 = int(sys.argv[2]) if len(sys.argv) > 2 else 3
N3 = int(sys.argv[3]) if len(sys.argv) > 3 else 2
NU = float(sys.argv[4]) if len(sys.argv) > 4 else 1.001
NWIN = int(sys.argv[5]) if len(sys.argv) > 5 else 16
m.NU_STAR = NU
m.NWIN = NWIN

print(f"mesh {N1}x{N2}x{N3} NU={NU} NWIN={NWIN}", flush=True)
ref = m.reference_arb(N1, N2, N3)
print(
    f"gamma={m.upper(ref['gamma']):.4f} S_Q={m.upper(ref['S_Q']):.5f} "
    f"S_M={m.upper(ref['S_M']):.4f}",
    flush=True,
)
glob = m.assemble_global(ref)
ng = glob["ng"]
print(f"dofs={ng} A_GB~{ng*ng*8/1e9:.2f}", flush=True)

volw = float(np.ones(ng) @ (glob["Mm"] @ np.ones(ng)))
t1, t2 = glob["tim"].sum(), glob["t0m"].sum()
a1 = glob["am"].sum()
sgrid = np.linspace(0.0, 1.0, NWIN + 1)


def slack_of(th, th2, al, th4, rt):
    slack = np.inf
    for w in range(NWIN):
        co = m.window_coeffs(
            ref, sgrid[w], sgrid[w + 1], th, th2, al, arb(rt), arb(NU), th4,
        )
        ce = float(co["c_e"].mid() - co["c_e"].rad())
        de = m.upper(co["d_e"])
        if float(co["c_S"].mid()) <= 0 or ce <= 0:
            return -1.0
        slack = min(slack, ce / de - 1)
        k0 = float(co["kap0"].mid())
        L1 = a1 + k0 * (t1 + t2)
        D1 = (
            float(co["lam_t"].mid()) * (1 + al) * volw
            + float(co["bt_inf"].mid()) * t1 ** 2
            + float(co["bt_0"].mid()) * t2 ** 2
        )
        slack = min(slack, rt * (1 - th) * L1 ** 2 / D1 - 1)
    return slack


def pencil(th, th2, al, th4, rt, w):
    co = m.window_coeffs(
        ref, sgrid[w], sgrid[w + 1], th, th2, al, arb(rt), arb(NU), th4,
    )
    cQ = float(co["c_Q"].mid())
    k0 = float(co["kap0"].mid())
    ell = glob["am"] + k0 * (glob["tim"] + glob["t0m"])
    rt1 = rt * (1 - th)
    lt = float(co["lam_t"].mid()) * (1 + al)
    bi = float(co["bt_inf"].mid())
    b0 = float(co["bt_0"].mid())
    Nsp = (cQ * (glob["Qm"] + (1 - th4) * glob["Rm"])).tocsr()
    Dsp = (lt * glob["Mm"]).tocsr()

    def nmv(x):
        x = np.asarray(x).ravel()
        return Nsp @ x + rt1 * (ell @ x) * ell

    def dmv(x):
        x = np.asarray(x).ravel()
        return (
            Dsp @ x
            + bi * (glob["tim"] @ x) * glob["tim"]
            + b0 * (glob["t0m"] @ x) * glob["t0m"]
        )

    prec = splu((Nsp + Dsp).tocsc())
    X = np.random.default_rng(1).standard_normal((ng, 6))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        vals, _ = lobpcg(
            LinearOperator((ng, ng), matvec=nmv),
            X,
            B=LinearOperator((ng, ng), matvec=dmv),
            M=LinearOperator((ng, ng), matvec=prec.solve),
            largest=False,
            tol=1e-6,
            maxiter=400,
        )
    return float(np.sort(vals)[0]), float(co["c_Q"].mid())


cands = []
for th in (0.3, 0.35, 0.4, 0.45, 0.5):
    for th2 in (0.9, 0.95, 0.98):
        for al in (0.1, 0.12, 0.15, 0.2):
            for th4 in (0.85, 0.9):
                for rt in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2):
                    s = slack_of(th, th2, al, th4, rt)
                    if s > 0:
                        cands.append((s, th, th2, al, th4, rt))
print(f"scalar-feasible {len(cands)}", flush=True)

# trial: low rt first + highest slack
trial = sorted(cands, key=lambda c: (c[5], -c[0]))[:12]
trial += sorted(cands, reverse=True, key=lambda c: c[0])[:4]
seen = set()
best = None
for c in trial:
    key = c[1:]
    if key in seen:
        continue
    seen.add(key)
    evs = []
    cQ0 = None
    for w in (0, NWIN // 2, NWIN - 1):
        ev, cQ = pencil(*c[1:], w)
        evs.append(ev)
        if w == 0:
            cQ0 = cQ
    ev = min(evs)
    score = min(c[0], ev - 1)
    print(
        f"params={c[1:]} slack={c[0]:.3f} ev={ev:.4f} "
        f"c_Q0={cQ0:.4f} score={score:.4f}",
        flush=True,
    )
    if best is None or score > best[0]:
        best = (score, c[1:], ev, c[0])
print("BEST", best, flush=True)
