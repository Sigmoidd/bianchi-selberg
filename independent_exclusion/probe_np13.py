"""Parameter probe for Gamma_0(3+2i), N(p)=13. Worktree np9 only."""
import warnings
import numpy as np
from flint import arb
from scipy.sparse.linalg import splu, lobpcg, LinearOperator

import congruence_prototype as cp
import m3p_certify as m

cp.set_level("(3+2i)")

import sys
from collections import defaultdict

NU = float(sys.argv[1]) if len(sys.argv) > 1 else 1.01
NWIN = int(sys.argv[2]) if len(sys.argv) > 2 else 12
N1 = int(sys.argv[3]) if len(sys.argv) > 3 else 6
N2 = int(sys.argv[4]) if len(sys.argv) > 4 else 3
N3 = int(sys.argv[5]) if len(sys.argv) > 5 else 2

m.NU_STAR = NU
m.NWIN = NWIN

print(f"building ref {N1}x{N2}x{N3} nu*={NU} NWIN={NWIN}...", flush=True)
ref = m.reference_arb(N1, N2, N3)
print(
    f"  gamma={m.upper(ref['gamma']):.4f} "
    f"alpha={m.upper(ref['alpha_ref']):.4f} "
    f"S_M={m.upper(ref['S_M']):.3f} "
    f"S_Q={m.upper(ref['S_Q']):.2e}",
    flush=True,
)
glob = m.assemble_global(ref)
volw = float(np.ones(glob["ng"]) @ (glob["Mm"] @ np.ones(glob["ng"])))
t1, t2 = glob["tim"].sum(), glob["t0m"].sum()
a1 = glob["am"].sum()
print(f"  volw={volw:.4f} t_inf={t1:.2f} t0={t2:.2f} a1={a1:.4f}", flush=True)

sgrid = np.linspace(0.0, 1.0, NWIN + 1)
cands = []
for th in (0.4, 0.45, 0.5, 0.55, 0.6, 0.65):
    for th2 in (0.7, 0.8, 0.9, 0.95):
        for al in (0.12, 0.15, 0.2, 0.25, 0.35):
            for th4 in (0.8, 0.85, 0.9):
                for rt in (0.8, 1.0, 1.15, 1.25, 1.5, 2.0, 2.5):
                    slack = np.inf
                    ok = True
                    for w in range(NWIN):
                        co = m.window_coeffs(
                            ref, sgrid[w], sgrid[w + 1],
                            th, th2, al, arb(rt), arb(NU), th4,
                        )
                        ce = float(co["c_e"].mid() - co["c_e"].rad())
                        de = m.upper(co["d_e"])
                        if float(co["c_S"].mid()) <= 0 or ce <= 0:
                            ok = False
                            break
                        slack = min(slack, ce / de - 1)
                        k0 = float(co["kap0"].mid())
                        L1 = a1 + k0 * (t1 + t2)
                        D1 = (
                            float(co["lam_t"].mid()) * (1 + al) * volw
                            + float(co["bt_inf"].mid()) * t1 ** 2
                            + float(co["bt_0"].mid()) * t2 ** 2
                        )
                        slack = min(slack, rt * (1 - th) * L1 ** 2 / D1 - 1)
                    if ok and slack > 0:
                        cands.append((slack, th, th2, al, th4, rt))

print(f"scalar-feasible {len(cands)}", flush=True)
by = defaultdict(list)
for c in cands:
    by[c[5]].append(c)
for rt in sorted(by):
    best = max(by[rt], key=lambda c: c[0])
    print(f"  best slack rt={rt}: {best[0]:.3f} {best[1:]}", flush=True)


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
    ng = glob["ng"]

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
    return float(np.sort(vals)[0])


if not by:
    print("NO SCALAR FEASIBLE")
else:
    best_score = None
    # low rt first + a few mid; also try top-slack overall
    trial = []
    for rt in sorted(by):
        trial.append(max(by[rt], key=lambda c: c[0]))
    # also top 5 by slack regardless of rt
    for c in sorted(cands, reverse=True, key=lambda c: c[0])[:5]:
        if c not in trial:
            trial.append(c)
    for c in trial:
        # check endpoints and a middle window
        ws = (0, NWIN // 2, NWIN - 1)
        ev = min(pencil(*c[1:], w) for w in ws)
        score = min(c[0], ev - 1)
        print(
            f"PENCIL params={c[1:]} slack={c[0]:.3f} "
            f"ev={ev:.4f} score={score:.4f}",
            flush=True,
        )
        if best_score is None or score > best_score[0]:
            best_score = (score, c[1:], ev, c[0])
    print("BEST", best_score, flush=True)
