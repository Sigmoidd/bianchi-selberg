"""Export a finite generated image and its regular adjacency permutations."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from groups.identity import LevelIdeal
from quotients import generated_action

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("d", type=int)
    parser.add_argument("--level-generator", type=int, nargs=2, default=(3, 0))
    parser.add_argument("--max-vertices", type=int, default=10000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        ideal = LevelIdeal.principal(args.d, tuple(args.level_generator))
        action = generated_action(ideal, max_vertices=args.max_vertices)
    except ValueError as exc:
        parser.error(str(exc))
    args.output.write_text(json.dumps(action.payload(), indent=2)+"\n")
    print(f"Wrote {args.output}: {action.order} vertices, degree {action.degree}.")
