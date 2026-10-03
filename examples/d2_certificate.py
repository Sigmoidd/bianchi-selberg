"""Replay the d=2 inventory proof and export the full Arb certificate."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core.assemble import evaluate
from core.certificate import certificate_payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--R", type=float, default=40)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    E = evaluate(2, k=args.k, R=args.R, verbose=False)
    payload = certificate_payload(E)
    text = json.dumps(payload, indent=2)+"\n"
    if args.output:
        args.output.write_text(text)
        print(f"Wrote {args.output}; B = {E.bound}, upper endpoint < 1.")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
