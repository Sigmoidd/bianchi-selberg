from __future__ import annotations

import unittest

from flint import arb, ctx

from track_b_floor_adaptive import _child_specs, _choose, _summary


def record(index: list[int], remainder: str) -> dict[str, object]:
    return {
        "cell_index": index,
        "adaptive_id": "base:" + ",".join(str(q) for q in index),
        "polynomial_l2_squared_upper": "1",
        "remainder_l2_upper": remainder,
        "combined_cell_l2_upper": "2",
        "_adaptive_bounds": [["0", "1"], ["0", "1"], ["0", "1"]],
    }


class TrackBFloorAdaptiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        ctx.prec = 160

    def test_contribution_selection_is_minimal_prefix(self) -> None:
        rows = [record([k, 0, 0], str(q)) for k, q in enumerate((4, 3, 1))]
        selected = _choose(rows, 0.0, 0.90, None)
        self.assertEqual(len(selected), 2)  # (16+9)/(16+9+1) > .90

    def test_child_partition_has_exact_cartesian_cardinality(self) -> None:
        parent = record([0, 0, 0], "1")
        children = _child_specs(parent, (0, 2))
        self.assertEqual(len(children), 4)
        self.assertEqual(len({q["adaptive_id"] for q in children}), 4)
        volumes = []
        for child in children:
            volume = arb(1)
            for lo, hi in child["bounds"]:
                volume *= hi - lo
            volumes.append(volume)
        self.assertTrue((sum(volumes, arb(0)) - 1).contains(0))

    def test_summary_uses_disjoint_l2_squares(self) -> None:
        rows = [record([0, 0, 0], "3"), record([1, 0, 0], "4")]
        result = _summary(rows, arb(10))
        self.assertTrue((arb(result["remainder_l2_upper"]) - 5).contains(0))
        self.assertTrue((arb(result["polynomial_l2_upper"]) ** 2 - 2).contains(0))


if __name__ == "__main__":
    unittest.main()
