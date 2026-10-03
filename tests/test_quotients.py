from fractions import Fraction
import subprocess
import sys
import unittest
from groups.identity import LevelIdeal
from quotients import ResidueRing, ProjectiveMatrices, generated_action, standard_generators


class QuotientTests(unittest.TestCase):
    def test_import_has_no_analytic_dependencies(self):
        subprocess.run([sys.executable, '-c',
            "import quotients,sys; assert not any(x == 'flint' or x.startswith('core') for x in sys.modules)"], check=True)

    def test_split_primes_and_inert_prime_images(self):
        P, Q = (LevelIdeal.principal(7, x) for x in ((0, 1), (1, -1)))
        self.assertNotEqual(P, Q)
        self.assertEqual(ResidueRing(P).reduce((0, 1)), (0, 0))
        self.assertEqual(ResidueRing(Q).reduce((0, 1)), (1, 0))
        for I, size in ((P, 6), (Q, 6), (LevelIdeal.rational(7, 3), 360)):
            action = generated_action(I)
            self.assertEqual(action.order, size)
            self.assertTrue(action.symmetric)
            self.assertEqual(action.degree, 6)
            self.assertEqual(action.adjacency_product([1]*size), (6,)*size)
            self.assertEqual(action.transition_product([Fraction(1)]*size), (1,)*size)
            self.assertTrue(all(sorted(p) == list(range(size)) for p in action.permutations))
            # Symmetry includes multiplicities, giving a symmetric operator.
            counts = {}
            for p in action.permutations:
                for i, j in enumerate(p):
                    counts[i,j] = counts.get((i,j), 0)+1
            self.assertTrue(all(count == counts.get((j,i), 0) for (i,j), count in counts.items()))

    def test_ring_cosets_composite_level_and_signs(self):
        I = LevelIdeal.principal(7, (0, 1))
        ring = ResidueRing(I)
        for x in ((-17, 9), (12, -5)):
            for v in I.basis:
                self.assertEqual(ring.reduce(x), ring.reduce((x[0]+v[0], x[1]+v[1])))
        I = LevelIdeal.rational(7, 5)
        matrices = ProjectiveMatrices(ResidueRing(I))
        for g in standard_generators():
            G = matrices.reduce(g)
            self.assertEqual(G, matrices.reduce(tuple((-a,-b) for a,b in g)))
            self.assertEqual(matrices.mul(G, matrices.inverse(G)), matrices.identity)
        self.assertEqual(generated_action(LevelIdeal.rational(7, 2)).order, 36)

    def test_size_cap_and_operator_input_fail_closed(self):
        with self.assertRaisesRegex(ValueError, 'size cap'):
            generated_action(LevelIdeal.rational(7, 3), max_vertices=100)
        action = generated_action(LevelIdeal.principal(7, (0,1)))
        with self.assertRaises(ValueError):
            action.adjacency_product([1])
        with self.assertRaises(ValueError):
            generated_action(action.ideal, generators=[])
        with self.assertRaises(ValueError):
            ProjectiveMatrices(ResidueRing(action.ideal)).reduce(((1,0),(0,0),(0,0),(0,0)))
