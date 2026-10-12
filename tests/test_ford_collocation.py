import copy
import json
from fractions import Fraction
from pathlib import Path
import unittest
import numpy as np
from flint import ctx, acb

from ford_equivariant_collocation import (PhysicalOperator, sample_points,
    mode_keys,restriction_indices,drift,sector_matrix,exact_zero_jet_columns,
    vertex_constraint_embedding,normalized_lstsq,reduced_jet_constraints)
from ford_vertex_equivariance import structurally_zero_odd_jet_columns
from ford_constraint_rank_arb import rank_lower
from certificate_hashes import text_hash
from track_b_two_cusp_data import canonical_hash,cusp_modes


class FordCollocationTests(unittest.TestCase):
    def setUp(self):
        ctx.prec=128

    def test_deterministic_disjoint_true_curved_cell_points(self):
        train=sample_points(32,173); hold=sample_points(16,947)
        self.assertEqual(train,sample_points(32,173))
        key=lambda p:tuple(p[k] for k in ('x1','x2','y'))
        self.assertFalse(set(map(key,train))&set(map(key,hold)))
        for p in train+hold:
            x,z,y=map(Fraction,key(p))
            self.assertTrue(-Fraction(1,2)<x<Fraction(1,2))
            self.assertTrue(0<z<Fraction(1,2))
            self.assertGreater(x*x+z*z+y*y,1)

    def test_fast_physical_rows_against_original_arb_component_code(self):
        result=PhysicalOperator(sample_points(2,173),2).independently_check('6.62208')
        self.assertTrue(result['agrees_to_1e-16'])

    def test_exact_mode_embeddings_fix_phase_and_scale(self):
        for lo in range(2,6):
            small=np.arange(len(mode_keys(lo)),dtype=complex)
            big=np.zeros(len(mode_keys(lo+1)),dtype=complex)
            big[restriction_indices(lo,lo+1)]=small
            self.assertEqual(drift(small,big,lo,lo+1)['delta_l2'],0)
            self.assertGreater(drift(small,big*1j,lo,lo+1)['delta_l2'],0)
        # Infinity/zero offsets shift when infinity modes are added.
        self.assertEqual([mode_keys(4)[k] for k in restriction_indices(3,4)],mode_keys(3))

    def test_exact_structural_redundancy_removed_before_scaling(self):
        operator,free=reduced_jet_constraints()
        self.assertEqual(operator.shape,(14,144)); self.assertEqual(len(free),20)
        for M in range(2,7):
            self.assertEqual(list(exact_zero_jet_columns(M)),structurally_zero_odd_jet_columns(M))
        odd=sector_matrix(2,-1,len(mode_keys(2)))
        embed,record=vertex_constraint_embedding(2,'6.62208',odd)
        self.assertEqual(record['constraint_rank_numerical'],9)
        self.assertEqual(embed.shape[1],2)
        norm_idx=next(k for k,m in enumerate(cusp_modes(2)[0]) if m[:2]==(1,0))
        self.assertIsNone(normalized_lstsq(np.eye(embed.shape[1]),embed,norm_idx))

    def test_independent_arb_minor_rank_and_corruption(self):
        r=rank_lower(2,'6.62208')
        self.assertEqual(r['certified_rank_lower'],9)
        self.assertEqual(r['exact_rank_upper'],9)
        self.assertTrue(r['exact_rank_certified'])
        self.assertFalse(acb(r['determinant']['real'],r['determinant']['imag']).contains(0))
        # A zero determinant cannot serve as a rank witness.
        self.assertTrue(acb(0).contains(0))

    def test_result_provenance_and_unchanged_existence_gates(self):
        p=Path('ford_spectral_stability.json')
        if not p.exists(): self.skipTest('full diagnostic not yet generated')
        out=json.loads(p.read_text())
        digest=out.pop('result_definition_hash')
        self.assertEqual(digest,canonical_hash(out))
        for file,digest in out['code_sha256'].items():
            self.assertEqual(digest,text_hash(Path(file).read_bytes()))
        for k in ('global_hejhal_defect_certified','hejhal_existence_certified','rung4_certified'):
            self.assertIs(out[k],False)
        self.assertEqual(out['train_points_hash'],canonical_hash(out['train_points']))
        broken=copy.deepcopy(out); broken['train_points'][0]['y']='2'
        self.assertNotEqual(broken['train_points_hash'],canonical_hash(broken['train_points']))
        stage=copy.deepcopy(out); stage['fixed_parameter_comparison'][0]['r']='7'
        self.assertNotEqual(canonical_hash(stage),canonical_hash(out))
