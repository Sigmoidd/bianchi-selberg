"""Serialize a regression or mechanical report, with its proof status."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core.assemble import evaluate
from core.certificate import report_payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fields", nargs="+", help="i, omega, or a supported d")
    parser.add_argument("--mechanical", action="store_true")
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--R", type=float, default=40)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    reports = []
    for kind in args.fields:
        kind = int(kind) if kind.isdigit() else kind
        E = evaluate(kind, k=args.k, R=args.R, verbose=False,
                     include_elliptic=not args.mechanical)
        reports.append(report_payload(E))
    text = json.dumps(dict(schema_version=1, reports=reports), indent=2)+"\n"
    if args.output:
        args.output.write_text(text)
        print(f"Wrote {args.output}; report only, not a new spectral certificate.")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
