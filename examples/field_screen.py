"""Arb mechanical screen; omitting elliptics never certifies a spectral gap."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core.assemble import evaluate
from core.bspline import g0
from groups.systoles import verify_systole


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("d", type=int, nargs="*", default=[2, 7, 11, 19])
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--R", type=float, default=40)
    args = parser.parse_args()
    print("MECHANICAL SCREEN ONLY: inventories are incomplete; no new spectral certificate.")
    for d in args.d:
        proof = verify_systole(d)
        E = evaluate(d, k=args.k, R=args.R, verbose=False, include_elliptic=False)
        budget = 1-E.bound
        print(f"d={d}: systole={proof['systole']}; B_mech={E.bound}; "
              f"NCE+CE budget={budget}; C_ell budget if CE=0: {budget/g0(E.k, E.delta)}")


if __name__ == "__main__":
    main()
