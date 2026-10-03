"""Replay the finite d=11 sinc-family comparison and export its selected certificate.

This optimizes the certified upper bound on an explicit finite grid. It does
not prove a global optimum over k, continuous delta or all test functions.
"""
import argparse
import json
import math
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core.assemble import evaluate
from core.certificate import certificate_payload

SCREEN_FRACTIONS = (.8, .9, .95, .975, .99, .999, 1.)
REFINEMENT_FRACTIONS = (.98, .99, .995, .999, 1.)


def run(output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    def sample(k, frac, R):
        E = evaluate(11, k=k, frac=frac, R=R, verbose=False)
        rows.append(dict(k=k, frac=frac, R=R, delta_hex=E.delta.hex(),
                         lower=math.nextafter(float(E.bound.lower()), -math.inf),
                         upper=math.nextafter(float(E.bound.upper()), math.inf),
                         bound_ball=E.bound.str(40)))
        print(k, frac, R, E.bound)
        return E
    for k in range(2, 9):
        for frac in SCREEN_FRACTIONS:
            sample(k, frac, 40)
    winner = min(rows, key=lambda row: row["upper"])
    # The committed refinement is deliberately tied to the observed winner.
    if (winner["k"], winner["frac"]) != (2, 1.):
        raise ArithmeticError("screen winner changed; review the refinement before exporting")
    for frac in REFINEMENT_FRACTIONS:
        sample(2, frac, 128)
    E = sample(2, 1., 256)
    if min(rows, key=lambda row: row["upper"]) is not rows[-1]:
        raise ArithmeticError("final parameters no longer minimize the sampled upper bounds")
    payload = dict(schema_version=1, group_d=11, family="sinc^(2k)",
                   objective="minimum certified upper bound among evaluated candidates",
                   global_optimality_proved=False,
                   screen=dict(k=list(range(2, 9)), fractions=list(SCREEN_FRACTIONS), R=40),
                   refinement=dict(k=2, fractions=list(REFINEMENT_FRACTIONS), R=128),
                   evaluations=rows, selected=dict(k=2, frac=1., R=256,
                       certificate="d11-k2.json", delta_hex=E.delta.hex()),
                   note="Finite parameter comparison, not a proof of continuous or global optimality.")
    for name, value in (("d11-k2.json", certificate_payload(E)), ("d11-search.json", payload)):
        (output_dir/name).write_text(json.dumps(value, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    run(parser.parse_args().output_dir)
