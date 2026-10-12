#!/usr/bin/env python3
"""Bounded M7/M8 follow-up justified by the Picard-region M5/M6 plateau."""
import json
from pathlib import Path

import ford_equivariant_collocation as diagnostic
from certificate_hashes import text_hash
from track_b_two_cusp_data import canonical_hash


def run_extension():
    primary=json.loads(Path('ford_spectral_stability.json').read_text())
    diagnostic.REGIONS={'picard_region':diagnostic.REGIONS['picard_region']}
    out=diagnostic.run([7,8],32,16,17)
    out['extension_justification']='M5/M6 Picard odd parameter agreement 2.56e-10 and coefficient delta0.000603; test two new cutoffs'
    out['primary_result_definition_hash']=primary['result_definition_hash']
    out['boundary_drift']={}
    for r in diagnostic.FIXED_R:
        for name in ('baseline','odd','odd_vertex_jets'):
            lo=next(q for q in primary['fixed_parameter_comparison'] if q['M']==6 and q['r']==r)
            hi=next(q for q in out['fixed_parameter_comparison'] if q['M']==7 and q['r']==r)
            out['boundary_drift'][f'fixed:{r}:{name}']=diagnostic.drift_records([lo,hi],lambda q:q['systems'][name])
    for name in ('baseline','odd','odd_vertex_jets'):
        lo=next(q for q in primary['spectral_scans'] if q['M']==6 and q['region']=='picard_region' and q['system']==name)
        hi=next(q for q in out['spectral_scans'] if q['M']==7 and q['region']=='picard_region' and q['system']==name)
        out['boundary_drift'][f'scan:picard_region:{name}']=diagnostic.drift_records([lo,hi],lambda q:q)
    out['code_sha256_policy']='sha256-utf8-lf/v1'
    out['code_sha256']={p:text_hash(Path(p).read_bytes()) for p in
        ('ford_equivariant_collocation.py','ford_vertex_equivariance.py','global_quarter_turn.py',
         'track_b_two_cusp_data.py','track_b_two_cusp_hejhal.py','ford_spectral_extension.py')}
    out['result_definition_hash']=canonical_hash(out)
    Path('ford_spectral_extension.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__=='__main__': run_extension()
