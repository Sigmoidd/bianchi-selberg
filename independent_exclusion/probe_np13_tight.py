"""Near-tight certification params for N(p)=13."""
import warnings
import numpy as np
from flint import arb
from scipy.sparse.linalg import splu, lobpcg, LinearOperator

import congruence_prototype as cp
import m3p_certify as m

cp.set_level("(3+2i)")
ref = m.reference_arb(6, 3, 2)
glob = m.assemble_global(ref)
ng = glob["ng"]
Y = 1.25

lam = 0.999
s = np.sqrt(1 - lam)
b_inf = 2 * (1 - s) / Y ** 2
b_0 = b_inf / 13
ell = glob["am"] + (glob["tim"] + glob["t0m"]) / ((1 + s) * Y ** 2)
Ash = (glob["Qm"] - lam * glob["Mm"]).tocsr()


def matvec(x):
    x = np.asarray(x).ravel()
    y = Ash @ x
    y = y - b_inf * (glob["tim"] @ x) * glob["tim"]
    y = y - b_0 * (glob["t0m"] @ x) * glob["t0m"]
    return y


prec = splu((Ash + 2 * glob["Mm"]).tocsc())
Yc = splu(glob["Mm"].tocsc()).solve(ell).reshape(-1, 1)
X = np.random.default_rng(0).standard_normal((ng, 6))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    vals, _ = lobpcg(
        LinearOperator((ng, ng), matvec=matvec),
        X,
        B=glob["Mm"],
        M=LinearOperator((ng, ng), matvec=prec.solve),
        Y=Yc,
        largest=False,
        tol=1e-7,
        maxiter=400,
    )
print("raw constrained mu", np.sort(vals)[:3], flush=True)

volw = float(np.ones(ng) @ (glob["Mm"] @ np.ones(ng)))
t1, t2 = glob["tim"].sum(), glob["t0m"].sum()
a1 = glob["am"].sum()


def eval_params(NU, NWIN, th, th2, al, th4, rt):
    sgrid = np.linspace(0.0, 1.0, NWIN + 1)
    slack = np.inf
    for w in range(NWIN):
        co = m.window_coeffs(
            ref, sgrid[w], sgrid[w + 1], th, th2, al, arb(rt), arb(NU), th4,
        )
        ce = float(co["c_e"].mid() - co["c_e"].rad())
        de = m.upper(co["d_e"])
        if float(co["c_S"].mid()) <= 0 or ce <= 0:
            return -1, None
        slack = min(slack, ce / de - 1)
        k0 = float(co["kap0"].mid())
        L1 = a1 + k0 * (t1 + t2)
        D1 = (
            float(co["lam_t"].mid()) * (1 + al) * volw
            + float(co["bt_inf"].mid()) * t1 ** 2
            + float(co["bt_0"].mid()) * t2 ** 2
        )
        slack = min(slack, rt * (1 - th) * L1 ** 2 / D1 - 1)

    def pencil(w):
        co = m.window_coeffs(
            ref, sgrid[w], sgrid[w + 1], th, th2, al, arb(rt), arb(NU), th4,
        )
        cQ = float(co["c_Q"].mid())
        k0 = float(co["kap0"].mid())
        ellv = glob["am"] + k0 * (glob["tim"] + glob["t0m"])
        rt1 = rt * (1 - th)
        lt = float(co["lam_t"].mid()) * (1 + al)
        bi = float(co["bt_inf"].mid())
        b0v = float(co["bt_0"].mid())
        Nsp = (cQ * (glob["Qm"] + (1 - th4) * glob["Rm"])).tocsr()
        Dsp = (lt * glob["Mm"]).tocsr()

        def nmv(x):
            x = np.asarray(x).ravel()
            return Nsp @ x + rt1 * (ellv @ x) * ellv

        def dmv(x):
            x = np.asarray(x).ravel()
            return (
                Dsp @ x
                + bi * (glob["tim"] @ x) * glob["tim"]
                + b0v * (glob["t0m"] @ x) * glob["t0m"]
            )

        pr = splu((Nsp + Dsp).tocsc())
        X0 = np.random.default_rng(1).standard_normal((ng, 8))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            vv, _ = lobpcg(
                LinearOperator((ng, ng), matvec=nmv),
                X0,
                B=LinearOperator((ng, ng), matvec=dmv),
                M=LinearOperator((ng, ng), matvec=pr.solve),
                largest=False,
                tol=1e-7,
                maxiter=500,
            )
        return float(np.sort(vv)[0]), float(co["c_Q"].mid()), float(
            co["c_e"].mid()
        ), float(co["d_e"].mid())

    evs = []
    for w in (0, NWIN // 2, NWIN - 1):
        ev, cQ, ce, de = pencil(w)
        evs.append(ev)
        print(
            f"  w={w} ev={ev:.4f} c_Q={cQ:.4f} c_e={ce:.3f} d_e={de:.3f}",
            flush=True,
        )
    return slack, min(evs)


candidates = [
    (1.001, 24, 0.30, 0.98, 0.08, 0.90, 0.45),
    (1.001, 24, 0.32, 0.97, 0.10, 0.90, 0.50),
    (1.001, 20, 0.35, 0.95, 0.10, 0.85, 0.50),
    (1.002, 24, 0.30, 0.98, 0.08, 0.88, 0.40),
    (1.001, 32, 0.28, 0.98, 0.08, 0.90, 0.42),
    (1.0005, 24, 0.30, 0.98, 0.10, 0.90, 0.48),
    (1.001, 16, 0.35, 0.95, 0.10, 0.85, 0.50),
]

best = None
for cand in candidates:
    NU, NWIN, th, th2, al, th4, rt = cand
    print(f"\nTRY NU={NU} NWIN={NWIN} params={(th, th2, al, th4, rt)}", flush=True)
    slack, ev = eval_params(NU, NWIN, th, th2, al, th4, rt)
    if slack is None or slack < 0:
        print(f"  scalar FAIL slack={slack}", flush=True)
        continue
    score = min(slack, (ev if ev is not None else 0) - 1)
    print(f"  slack={slack:.3f} min_ev={ev:.4f} score={score:.4f}", flush=True)
    if best is None or score > best[0]:
        best = (score, cand, slack, ev)

print("\nBEST", best, flush=True)
