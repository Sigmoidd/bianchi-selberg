"""Exact Klein-four witness in the PSL centralizer of the Eisenstein involution.

Works from any working directory. A bounded census is diagnostic and does not
prove that the finite subgroup is maximal or settle the orbital normalization.
"""
import argparse
import importlib.util
from pathlib import Path
import sys

REPO = next((p for p in Path(__file__).resolve().parents
             if (p/"elliptic_inventory.py").is_file()), None)


def load_inventory(repo=REPO):
    if repo is None:
        raise FileNotFoundError("Cannot discover a checkout; use --repo CHECKOUT")
    path = Path(repo).resolve()/"elliptic_inventory.py"
    if not path.is_file():
        raise FileNotFoundError(f"No elliptic_inventory.py at {path}; use --repo CHECKOUT")
    spec = importlib.util.spec_from_file_location("bianchi_inventory_evidence", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.path.insert(0, str(path.parent))
    from groups.arithmetic import QuadraticRing
    module.ExactRing = QuadraticRing
    return module


def centralizer_census(E, kind, M, B=1):
    if isinstance(B, bool) or not isinstance(B, int) or not 1 <= B <= 3:
        raise ValueError("census bound must be an integer in [1, 3]")
    ring = E.ExactRing(kind)
    MM = E.M2(ring)
    exact, flips = set(), set()
    for X in E.conjugators(ring, B=B):
        XM, MX = MM.mul(X, M), MM.mul(M, X)
        if XM == MX:
            exact.add(MM.canon(X))
        elif XM == MM.neg(MX):
            flips.add(MM.canon(X))
    return exact, flips


def klein_four_witness(E):
    ring = E.ExactRing("omega")
    MM = E.M2(ring)
    I = (E.ONE, E.ZERO, E.ZERO, E.ONE)
    R = (E.ZERO, (-1, 0), E.ONE, E.ZERO)
    w, w2 = (0, 1), (-1, -1)
    X = (w, w2, w2, ring.neg(w))
    conditions = (
        MM.det(R) == E.ONE,
        MM.det(X) == E.ONE,
        MM.mul(R, R) == MM.neg(I),
        MM.mul(X, X) == MM.neg(I),
        MM.mul(X, R) == MM.neg(MM.mul(R, X)),
    )
    if not all(conditions):
        raise ArithmeticError("the exact determinant/order/PSL commutation checks failed")
    H = {MM.canon(Y) for Y in (I, R, X, MM.mul(R, X))}
    if len(H) != 4 or any(MM.canon(MM.mul(A, B)) not in H for A in H for B in H):
        raise ArithmeticError("the proposed Klein-four subgroup is not closed of order four")
    if any(not MM.is_pmI(MM.mul(A, A)) for A in H):
        raise ArithmeticError("a proposed involution does not square to identity in PSL")
    return R, X, H


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=REPO)
    parser.add_argument("--bound", type=int, default=1, choices=(1, 2, 3))
    args = parser.parse_args(argv)
    E = load_inventory(args.repo)
    R, X, H = klein_four_witness(E)
    exact, flips = centralizer_census(E, "omega", R, args.bound)
    if E.M2(E.ExactRing("omega")).canon(X) not in flips:
        raise ArithmeticError("explicit flip is missing from the bounded census")
    print("PASS: det(R)=det(X)=1; R^2=X^2=-I; XR=-RX in SL2.")
    print(f"PASS: {{I,R,X,RX}} is a closed Klein-four subgroup of order {len(H)} in PSL2.")
    print(f"Bound {args.bound}: {len(exact)} exact-commuting and {len(flips)} sign-flipping PSL elements.")
    print("PROVED: the full PSL centralizer contains a finite subgroup of order 4.")
    print("NOT PROVED: maximality, primitive loxodromic norm, class completeness, or orbital normalization.")
    # A trace-one lift cannot satisfy XMX^-1=-M: conjugation preserves its
    # trace, whereas -M has trace -1. The census illustrates that distinction.
    _, elements = E.enumerate_elliptic(E.ExactRing("i"), 1, B=1)
    M3 = sorted(elements)[0]
    _, flips3 = centralizer_census(E, "i", M3, args.bound)
    if flips3:
        raise ArithmeticError("trace-one representative has an impossible sign flip")
    print("PASS: Picard trace-one order-3 contrast has no sign flips.")


if __name__ == "__main__":
    main()
