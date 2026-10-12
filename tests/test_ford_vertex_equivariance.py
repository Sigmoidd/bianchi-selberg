"""Exact geometry, transport, jet, and fail-closed regression tests."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

import ford_vertex_equivariance as v
from global_cusp_gaussian import QMAT, S, T1, TI, R, minv, mmul


class FordVertexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = v.build_certificate()

    def test_serialized_certificate_replays_exactly(self):
        self.assertEqual(self.cert, json.loads(Path(v.__file__).with_suffix('.json').read_text()))
        self.assertTrue(self.cert['all_exact_checks_pass'])

    def test_stabilizers_orbits_and_jet_dimensions(self):
        expected = [(4,[2,2,2],3,3,1,1),(6,[3,3],2,2,2,2),
                    (6,[3,3],2,2,2,2),(6,[3,3],2,2,2,2),
                    (12,[6],1,1,0,0),(12,[6],1,1,0,0)]
        for name, (order,orbits,dv,dg,ov,og) in zip(v.NAMES,expected):
            row=self.cert['vertices'][name]
            with self.subTest(name=name):
                self.assertEqual(row['stabilizer_order'],order)
                self.assertEqual(sorted(map(len,row['fiber_orbits'])),orbits)
                self.assertEqual(int(row['stabilizer_value_dimension']),dv)
                self.assertEqual(int(row['stabilizer_gradient_dimension']),dg)
                self.assertEqual(row['allowed_odd_value_dimension'],ov)
                self.assertEqual(row['allowed_odd_gradient_dimension'],og)
                self.assertTrue(all(row['allowed_odd_projector_checks'].values()))

    def test_action_and_jacobian_chain_rule(self):
        # Independent group-composition identity at off-stratum rational points.
        for p in [(F(2,7),F(-3,11),F(5,4)),(F(-1,3),F(1,6),F(7,9))]:
            for a in (S,T1,TI,R,QMAT):
                for b in (S,T1,QMAT):
                    bp,jb=v.action_t_jacobian(b,p)
                    abp,ja=v.action_t_jacobian(a,bp)
                    cp,jc=v.action_t_jacobian(mmul(a,b),p)
                    self.assertEqual(cp,abp)
                    self.assertEqual(jc,v.matmul(ja,jb))

    def test_all_q_images_and_reductions(self):
        for name,p in v.VERTICES.items():
            target,d=v.Q_REDUCTIONS[name]
            qp,_=v.action_t_jacobian(QMAT,p)
            self.assertEqual(qp,(-p[1],p[0],p[2]))
            reduced,_=v.action_t_jacobian(d,qp)
            self.assertEqual(reduced,v.VERTICES[target])
            self.assertEqual(v.mdet(d),v.ONE)
            self.assertEqual(v.evaluate_group_word(v.REDUCING_WORDS[name]),v.geo.canon(d))
            self.assertTrue(self.cert['vertices'][name]['reducing_transformation_in_H_pi'])
        for name in v.NAMES:
            for matrix,word in v.stabilizer_words()[name].items():
                self.assertEqual(v.evaluate_group_word(word),matrix)

    def test_reynolds_and_odd_projectors_annihilate_constraints(self):
        for comp in self.cert['simultaneous_components']:
            self.assertTrue(comp['constraint_times_projector_zero'])
            self.assertTrue(all(comp['projector_checks'].values()))
        self.assertEqual([c['simultaneous_odd_jet_dimension'] for c in self.cert['simultaneous_components']],[2,4,0])
        for name in v.NAMES:
            p=v.vertex_projector(name)
            for row in v.stabilizer_data(name):
                self.assertEqual(v.matmul(row['jet'],p),p)

    def test_scalar_nodes_do_not_annihilate_every_six_fiber_jet(self):
        for name in ('v_00','v_mhh','v_phh'):
            p=v.allowed_odd_vertex_projector(name)
            self.assertTrue(all(x==0 for i in (0,6,7,8) for x in p[i]))
        self.assertEqual(v.trace(v.allowed_odd_vertex_projector('v_00')),2)
        self.assertEqual(v.trace(v.allowed_odd_vertex_projector('v_mhh')),0)

    def test_integer_odd_fourier_nodal_jets_mode_by_mode(self):
        # At the three scalar nodes all C4 phases coincide (a sign); the
        # signed constant and horizontal frequency sums vanish exactly.
        for x,y in [(F(0),F(0)),(F(-1,2),F(1,2)),(F(1,2),F(1,2))]:
            for a in range(-4,5):
                for b in range(-4,5):
                    orbit=[(a,b),(-b,a),(-a,-b),(b,-a)]
                    phases=[(u*x+w*y)%1 for u,w in orbit]
                    self.assertEqual(len(set(phases)),1)
                    self.assertEqual(sum((-1)**k for k in range(4)),0)
                    self.assertEqual([sum((-1)**k*orbit[k][j] for k in range(4)) for j in range(2)],[0,0])

    def test_atlas_isometry_conjugacy_and_fail_closed(self):
        ledger=self.cert['atlas']
        self.assertTrue(ledger['finite_incidence_complete'])
        self.assertTrue(all(r['Q_conjugate_intertwining_exact'] for r in ledger['transition_words']))
        with self.assertRaises(ValueError): v.require_complete_transition_coverage(ledger)
        for key in ('complete_continuum_transition_cover','all_transition_defects_covered'):
            altered=copy.deepcopy(ledger); altered[key]=True
            with self.assertRaises(ValueError): v.require_complete_transition_coverage(altered)
        for flag in ('spectral_certificate','global_hejhal_defect_certified','rung4_certified','hejhal_existence_certified'):
            self.assertIs(self.cert[flag],False)

    def test_universal_added_vertex_rank_bound(self):
        data=self.cert['redundant_fourier_jet_certificate']
        self.assertTrue(data['exact_checks_pass'])
        self.assertEqual(data['parabolic_odd_ambient_jet_dimension'],20)
        self.assertEqual(data['universal_added_coefficient_row_rank_upper'],14)
        self.assertEqual([r['additional_rank'] for r in data['components']],[2,8,4])
        structural=v.constraint_rows(parabolic_only=True)
        reduced=v.reduced_vertex_constraint_rows()
        self.assertEqual(len(reduced),14)
        self.assertEqual(v.rref_sparse(structural+reduced),v.rref_sparse(v.constraint_rows()))

    def test_exact_zero_fourier_jet_columns_and_provenance(self):
        self.assertEqual(v.structurally_zero_odd_jet_columns(2),[1,9])
        self.assertEqual(v.structurally_zero_odd_jet_columns(3),[1,9])
        self.assertEqual([r['exact_vertex_constraint_rank_upper'] for r in v.cutoff_structural_rank_bounds()],[9,11,14,14,14])
        self.assertEqual(self.cert['provenance']['source_commit'],v.SOURCE_COMMIT)
        self.assertFalse(self.cert['provenance']['historical_digests_rewritten'])
        for evidence in self.cert['provenance']['artifacts'].values():
            self.assertEqual(evidence['canonical_hash_policy'],'sha256-utf8-lf/v1')
            self.assertEqual(len(evidence['raw_sha256']),64)

if __name__=='__main__': unittest.main()
