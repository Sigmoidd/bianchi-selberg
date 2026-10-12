#!/usr/bin/env python3
"""Freeze one EXPERIMENTAL candidate and certify a conditional cusp mass witness.

This is not track-b-frozen-trial/v1 and does not enable a global consumer.
The complete atlas, Dom(Delta) construction and residual are still blocked.
"""
from fractions import Fraction
import json
from pathlib import Path

from flint import acb,arb,ctx
import numpy as np

from certificate_hashes import text_hash
from track_b_two_cusp_data import canonical_hash,cusp_modes


def assess(primary,extension):
    records=primary['spectral_scans']+extension['spectral_scans']
    out={'schema':'ford-numerical-stability-assessment/v1',
         'criterion':{'last_two_coefficient_deltas_below':0.1,'last_two_r_changes_below':1e-4,
                      'heldout_rms_below':1e-6,'interior_minima_required':True,'full_rank_required':True},
         'systems':{},'global_certificate_gate_passed':False}
    def hashed(data):
        payload=dict(data); expected=payload.pop('result_definition_hash',None)
        return expected==canonical_hash(payload)
    out['provenance_checks']={
        'primary_hash':hashed(primary),'extension_hash':hashed(extension),
        'primary_extension_link':extension.get('primary_result_definition_hash')==primary.get('result_definition_hash'),
        'identical_physical_points':primary.get('train_points_hash')==extension.get('train_points_hash')
             and primary.get('heldout_points_hash')==extension.get('heldout_points_hash')}
    for name in ('baseline','odd','odd_vertex_jets'):
        pts=sorted([q for q in records if q['region']=='picard_region' and q['system']==name and q['M']>=6],key=lambda q:q['M'])
        ds=extension['boundary_drift']['scan:picard_region:'+name]+extension['cutoff_drift']['scan:picard_region:'+name]
        tests={'exact_cutoff_records':[q['M'] for q in pts]==[6,7,8],
               'exact_transition_records':[(d['from_M'],d['to_M']) for d in ds]==[(6,7),(7,8)],
               'provenance':all(out['provenance_checks'].values()),
               'two_small_deltas':all(d['delta_l2']<0.1 for d in ds),
               'two_small_r_changes':all(d['r_change']<1e-4 for d in ds),
               'interior_minima':all(not q['boundary_minimum'] for q in pts),
               'observable_numerical_rank':all(q['numerical_rank']==q['unknowns'] for q in pts),
               'heldout_margin':all(q['heldout_physical_residual_rms']<1e-6 for q in pts)}
        out['systems'][name]={'numerical_gate_passed':all(tests.values()),'checks':tests,'last_two_drift':ds,
            'interpretation':'bounded numerical stabilization; not certified cutoff convergence',
            'includes_new_infinity_shell_at_M8':True}
    return out


def freeze_candidate(primary,extension):
    chosen=next(q for q in extension['spectral_scans'] if q['M']==8 and q['system']=='odd')
    r=Fraction(chosen['r_M'])
    coefficients=[[str(Fraction(str(re))),str(Fraction(str(im)))] for re,im in chosen['coefficient_vector']]
    mi,mz=cusp_modes(8); keys=[['infinity',a,b,n] for a,b,n in mi]+[['zero',a,b,n] for a,b,n in mz]
    idx={tuple(k[:3]):j for j,k in enumerate(keys)}
    normal=idx[('infinity',1,0)]
    assert coefficients[normal]==['1','0']
    for k,p in idx.items():
        cusp,a,b=k; q=idx[(cusp,-b,a)]
        assert all(Fraction(coefficients[p][t])==-Fraction(coefficients[q][t]) for t in (0,1))
    atlas=json.loads(Path('ford_vertex_equivariance.json').read_text())
    result={'schema':'ford-experimental-candidate/v1','tier':'EXPERIMENTAL frozen exact rational coefficients from a floating singular vector',
            'audited_bundle_base_commit':'b5f29bd1b314e4d385e5af01805729ba27e3f3e2',
            'remote_tree_identical_base_commit':'b337589612006e9e8752ab074e73914899ecd645',
            'r':str(r),'r_decimal':chosen['r_M'],'lambda_exact':str(1+r*r),
            'cutoff_M':8,'mode_order':keys,'coefficients_exact_rational_real_imag':coefficients,
            'coefficient_vector_hash':canonical_hash(coefficients),'normalization':'a_infinity,(1,0)=1 exactly',
            'odd_coefficient_relations_exact':True,'atlas_definition_hash':canonical_hash(atlas['atlas']),
            'atlas_incidence_file_sha256':text_hash(Path('ford_vertex_equivariance.json').read_bytes()),
            'partition_file_sha256':text_hash(Path('track_b_global_partition_result.json').read_bytes()),
            'partition_definition_hash':json.loads(Path('track_b_global_partition_result.json').read_text())['partition_definition_hash'],
            'source_code_sha256':extension['code_sha256'],
            'candidate_producer_sha256':text_hash(Path(__file__).read_bytes()),
            'primary_result_definition_hash':primary['result_definition_hash'],
            'extension_result_definition_hash':extension['result_definition_hash'],
            'global_trial_frozen':False,'verified_interval_solve':False,'global_field_Dom_Delta_certified':False,
            'global_atlas_complete':False,'projected_mass_certified':False,'global_residual_certified':False,
            'hejhal_existence_certified':False,
            'blockers':['mixed-chart continuum overlap cover','same-trial value/gradient continuum residual',
                        'global gluing and Dom(Delta) compatibility','same-trial full projected mass and tails']}
    result['candidate_id']=canonical_hash(result)
    return result


def cusp_witness(candidate,segments=128):
    ctx.prec=192
    r=arb(candidate['r_decimal']); a,b=Fraction(13,10),Fraction(14,10)
    width=(b-a)/segments; integral=arb(0); leaves=[]
    for k in range(segments):
        lo=a+k*width; hi=lo+width
        yl=arb(lo.numerator)/lo.denominator; yh=arb(hi.numerator)/hi.denominator
        ball=yl.union(yh)
        value=acb(2*arb.pi()*ball).bessel_k(acb(0,r))
        if not value.is_finite() or not value.imag.contains(0) or not value.real.lower()>0:
            raise ArithmeticError('cusp witness failed to prove strictly positive K_ir')
        lower=value.real.lower()
        term=2*(arb(width.numerator)/width.denominator)*lower*lower/yh
        integral+=term
        leaves.append({'y_lower':str(lo),'y_upper':str(hi),'K_real_interval':str(value.real),
                       'mass_squared_lower_contribution':str(term.lower())})
    mass_squared=integral.lower(); mass=mass_squared.sqrt().lower()
    out={'schema':'ford-cusp-mass-witness/v1','tier':'CERTIFIED conditional one-cusp witness; not a global projected-mass certificate',
         'candidate_id':candidate['candidate_id'],'coefficient_vector_hash':candidate['coefficient_vector_hash'],
         'r_decimal':candidate['r_decimal'],'cutoff_M':8,'arb_bits':192,'segments':segments,
         'y_interval':['13/10','7/5'],'mass_squared_lower':str(mass_squared),'mass_lower':str(mass),
         'measure':'dx1 dx2 dy/y^3; infinity cross-section Z[i] modulo z->-z has half-torus measure',
         'parseval_lower_formula':'2*integral_(1.3)^(1.4) K_ir(2*pi*y)^2 dy/y',
         'hypothesis':'a global odd extension agrees exactly with the candidate infinity Fourier field on 1.3<=y<=1.4',
         'embedding_proof':'nonparabolic H element has |c|^2>=5, hence yprime<=1/(5y)<1.3; cusp slice embeds',
         'four_modes_coefficient_moduli_exactly_one':True,'direct_Arb_evaluations':segments,'fallback_count':0,
         'cusp_witness_integral_certified':bool(mass>0),'global_field_constructed':False,
         'global_projected_mass_certified':False,'global_residual_certified':False,'hejhal_existence_certified':False,
         'leaves':leaves}
    out['witness_definition_hash']=canonical_hash(out)
    return out


if __name__=='__main__':
    primary=json.loads(Path('ford_spectral_stability.json').read_text())
    extension=json.loads(Path('ford_spectral_extension.json').read_text())
    assessment=assess(primary,extension)
    Path('ford_stability_assessment.json').write_text(json.dumps(assessment,indent=2)+'\n')
    if not assessment['systems']['odd']['numerical_gate_passed']:
        raise SystemExit('no candidate passes the declared bounded numerical gate')
    candidate=freeze_candidate(primary,extension)
    Path('ford_experimental_candidate.json').write_text(json.dumps(candidate,indent=2)+'\n')
    witness=cusp_witness(candidate)
    Path('ford_cusp_mass_witness.json').write_text(json.dumps(witness,indent=2)+'\n')
    print(json.dumps({'candidate_id':candidate['candidate_id'],'mass_lower':witness['mass_lower'],
                      'global_trial_frozen':False,'global_projected_mass_certified':False},indent=2))
