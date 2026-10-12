#!/usr/bin/env python3
"""Independent physical-row replay of numerical trials and strict manifest checks.

Uses the original Track-B Arb assembler and independently reconstructed Arb
vertex jets, rather than the diagnostic's vectorized row/jet assembly.
Validation of midpoint diagnostics is not a continuum certificate.
"""
import argparse
import copy
import json
import math
from pathlib import Path

from flint import arb,ctx
import numpy as np

from certificate_hashes import text_hash
from track_b_two_cusp_data import canonical_hash
from track_b_two_cusp_hejhal import (ValidatedWhittaker,assemble_physical_rows_from_points,
                                   midpoint_matrix)
from ford_constraint_rank_arb import reconstruct

ROOT=Path(__file__).resolve().parent


def identity_checks(data):
    payload=copy.deepcopy(data); digest=payload.pop('result_definition_hash',None)
    checks={'result_hash':digest==canonical_hash(payload),
            'train_points_hash':data.get('train_points_hash')==canonical_hash(data.get('train_points')),
            'heldout_points_hash':data.get('heldout_points_hash')==canonical_hash(data.get('heldout_points')),
            'existence_flags_fail_closed':all(data.get(k) is False for k in
                 ('global_hejhal_defect_certified','rung4_certified','hejhal_existence_certified'))}
    for p,digest in data.get('code_sha256',{}).items():
        file=ROOT/p
        checks['code:'+p]=file.is_file() and text_hash(file.read_bytes())==digest
    for rec in data.get('fixed_parameter_comparison',[]):
        for name,trial in rec['systems'].items():
            if 'coefficient_vector' in trial:
                checks[f"coeff:{rec['M']}:{rec['r']}:{name}"]=canonical_hash(trial['coefficient_vector'])==trial.get('coefficient_vector_hash')
    for rec in data.get('spectral_scans',[]):
        if 'coefficient_vector' in rec:
            checks[f"scan-coeff:{rec['M']}:{rec['region']}:{rec['system']}"]=canonical_hash(rec['coefficient_vector'])==rec.get('coefficient_vector_hash')
    return checks


def replay(data):
    ctx.prec=128
    checks=identity_checks(data); reports=[]
    for rec in data['fixed_parameter_comparison']:
        M,r=rec['M'],rec['r']
        train,_,mi,mz=assemble_physical_rows_from_points(M,data['train_points'],ValidatedWhittaker(arb(r)),None)
        hold,_,_,_=assemble_physical_rows_from_points(M,data['heldout_points'],ValidatedWhittaker(arb(r)),None)
        a,h=midpoint_matrix(train),midpoint_matrix(hold)
        ctx.prec=192
        constraints,_=reconstruct(M,r)
        cmat=midpoint_matrix(constraints)
        ctx.prec=128
        # reconstruct() returns the odd-basis matrix; coefficient restriction
        # is recovered from alternating orbit columns, not a numerical inverse.
        import global_quarter_turn as gq
        columns=gq.sector_basis(M,-1)
        norm_idx=next(k for k,m in enumerate(mi) if m[:2]==(1,0))
        for name,trial in rec['systems'].items():
            if 'coefficient_vector' not in trial:
                continue
            coeff=np.array([complex(*v) for v in trial['coefficient_vector']])
            train_rms=float(np.linalg.norm(a@coeff)/math.sqrt(a.shape[0]))
            hold_rms=float(np.linalg.norm(h@coeff)/math.sqrt(h.shape[0]))
            key=f"physical:{M}:{r}:{name}"
            checks[key]=bool(abs(coeff[norm_idx]-1)<1e-10 and
                math.isclose(train_rms,trial['train_physical_residual_rms'],rel_tol=1e-7,abs_tol=1e-14) and
                math.isclose(hold_rms,trial['heldout_physical_residual_rms'],rel_tol=1e-7,abs_tol=1e-14))
            record={'M':M,'r':r,'system':name,'train_rms':train_rms,'heldout_rms':hold_rms,
                    'physical_row_assembler':'original Track B Arb',
                    'recorded_midpoint_residuals_agree':checks[key]}
            if name!='baseline':
                reduced=np.array([sum(coeff[k]*v for k,v in col.items())/4 for col in columns])
                # B has four signed unit entries per disjoint orbit.
                embedded=np.zeros(len(coeff),dtype=complex)
                for v,col in zip(reduced,columns):
                    for k,sgn in col.items(): embedded[k]=sgn*v
                checks[key+':odd']=bool(np.linalg.norm(embedded-coeff)<1e-10*max(np.linalg.norm(coeff),1))
                if name=='odd_vertex_jets':
                    err=float(np.linalg.norm(cmat@reduced))
                    record['independent_vertex_constraint_defect']=err
                    checks[key+':vertex']=err<1e-12*max(np.linalg.norm(coeff),1)
            reports.append(record)
    # Selected scan vectors receive a fresh held-out physical reconstruction.
    for rec in data['spectral_scans']:
        if rec['M']!=max(data['cutoffs']) or 'coefficient_vector' not in rec:
            continue
        coeff=np.array([complex(*v) for v in rec['coefficient_vector']])
        rows,_,_,_=assemble_physical_rows_from_points(rec['M'],data['heldout_points'],
                                                     ValidatedWhittaker(arb(rec['r_M'])),None)
        rms=float(np.linalg.norm(midpoint_matrix(rows)@coeff)/math.sqrt(len(rows)))
        key=f"scan-heldout:{rec['region']}:{rec['system']}"
        checks[key]=math.isclose(rms,rec['heldout_physical_residual_rms'],rel_tol=1e-7,abs_tol=1e-14)
        reports.append({'M':rec['M'],'r':rec['r_M'],'system':rec['system'],'region':rec['region'],
                        'heldout_rms':rms,'recorded_midpoint_residuals_agree':checks[key]})
    mutations={}
    for name in ('parameter','coefficient','geometry-source','existence-flag'):
        bad=copy.deepcopy(data)
        if name=='parameter': bad['fixed_parameter_comparison'][0]['r']='7'
        elif name=='coefficient': bad['fixed_parameter_comparison'][0]['systems']['odd']['coefficient_vector'][0][0]+=1
        elif name=='geometry-source': bad['code_sha256']['ford_vertex_equivariance.py']='0'*64
        else: bad['hejhal_existence_certified']=True
        mutations[name]=not all(identity_checks(bad).values())
    checks={k:bool(v) for k,v in checks.items()}
    return {'schema':'ford-stability-independent-verification/v1','verified':all(checks.values()) and all(mutations.values()),
            'checks':checks,'reconstructed_trials':reports,'deliberate_corruptions_rejected':mutations,
            'numerical_diagnostics_only':True,'continuum_residual_certified':False,'hejhal_existence_certified':False}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'ford_spectral_stability.json')
    parser.add_argument('--output',type=Path,default=ROOT/'ford_stability_independent_verification.json')
    ns=parser.parse_args()
    out=replay(json.loads(ns.input.read_text()))
    ns.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'verified':out['verified'],'checks':len(out['checks']),
                      'failed':[k for k,v in out['checks'].items() if not v]},indent=2))
    raise SystemExit(0 if out['verified'] else 1)
