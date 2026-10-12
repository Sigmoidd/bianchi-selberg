import copy
import json
from pathlib import Path
import unittest
import numpy as np
from scipy.linalg import svd
from ford_candidate_audit import assess
from ford_cusp_mass_verify import identity,verify
from ford_stability_verify import identity_checks


class FordCandidateTests(unittest.TestCase):
    def test_fixed_metric_spectrum_is_invariant_under_kernel_rotations(self):
        rng=np.random.default_rng(811)
        a=rng.normal(size=(19,7))
        orthogonal=np.linalg.qr(rng.normal(size=(7,7)))[0]
        self.assertTrue(np.allclose(svd(a,compute_uv=False),svd(a@orthogonal,compute_uv=False),rtol=1e-12))
        equilibrate=lambda m:m/np.linalg.norm(m,axis=0)
        self.assertFalse(np.allclose(svd(equilibrate(a),compute_uv=False),svd(equilibrate(a@orthogonal),compute_uv=False),rtol=1e-3))

    def test_numeric_gate_is_distinct_from_global_certificate_gate(self):
        p=json.loads(Path('ford_spectral_stability.json').read_text())
        e=json.loads(Path('ford_spectral_extension.json').read_text())
        out=assess(p,e)
        self.assertTrue(out['systems']['odd']['numerical_gate_passed'])
        self.assertTrue(out['systems']['odd_vertex_jets']['numerical_gate_passed'])
        self.assertFalse(out['global_certificate_gate_passed'])
        missing=copy.deepcopy(e)
        missing['spectral_scans']=[q for q in missing['spectral_scans'] if q['M']!=7]
        self.assertFalse(assess(p,missing)['systems']['odd']['numerical_gate_passed'])
        bad=copy.deepcopy(p);bad['hejhal_existence_certified']=True
        self.assertFalse(all(identity_checks(bad).values()))

    def test_conditional_cusp_witness_replays_and_rejects_mismatched_trials(self):
        c=json.loads(Path('ford_experimental_candidate.json').read_text())
        w=json.loads(Path('ford_cusp_mass_witness.json').read_text())
        out=verify(c,w)
        self.assertTrue(out['verified'])
        self.assertFalse(out['global_projected_mass_certified'])
        for key in ('r_decimal','coefficient_vector_hash','candidate_id'):
            bad=copy.deepcopy(w);bad[key]='different'
            self.assertFalse(all(identity(c,bad).values()))
        bad=copy.deepcopy(w);bad['global_projected_mass_certified']=True
        self.assertFalse(all(identity(c,bad).values()))
        for key in ('r','lambda_exact'):
            bad=copy.deepcopy(c);bad[key]='1'
            self.assertFalse(all(identity(bad,w).values()))
