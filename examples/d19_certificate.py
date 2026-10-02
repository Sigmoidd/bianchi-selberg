"""Replay the selected d=19 sinc^4 certificate."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from examples.group_certificate import main

if __name__ == "__main__":
    main(default_field=19, default_frac=1.0, default_R=256)
