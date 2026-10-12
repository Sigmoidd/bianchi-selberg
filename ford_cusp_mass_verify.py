#!/usr/bin/env python3
"""Independent, finer Arb replay of the conditional infinity-cusp witness."""
from fractions import Fraction
import json
from pathlib import Path
from flint import acb,arb,ctx
from track_b_two_cusp_data import canonical_hash
from certificate_hashes import text_hash


def identity(candidate,witness):
    c=dict(candidate); cid=c.pop('candidate_id',None)
    w=dict(witness); digest=w.pop('witness_definition_hash',None)
    normal=next(i for i,key in enumerate(candidate['mode_order']) if key[:3]==['infinity',1,0])
    checks={'candidate_hash':canonical_hash(c)==cid,
            'witness_hash':canonical_hash(w)==digest,
            'same_candidate':witness['candidate_id']==cid,
            'coefficient_hash':canonical_hash(candidate['coefficients_exact_rational_real_imag'])==candidate['coefficient_vector_hash']==witness['coefficient_vector_hash'],
            'same_parameter':candidate['r_decimal']==witness['r_decimal'],
            'rational_decimal_parameter':Fraction(candidate['r'])==Fraction(candidate['r_decimal']),
            'exact_lambda_formula':Fraction(candidate['lambda_exact'])==1+Fraction(candidate['r'])**2,
            'witness_height_interval':witness.get('y_interval')==['13/10','7/5'],
            'schemas':candidate.get('schema')=='ford-experimental-candidate/v1'
                and witness.get('schema')=='ford-cusp-mass-witness/v1',
            'same_cutoff':candidate['cutoff_M']==witness['cutoff_M']==8,
            'exact_normalization':candidate['coefficients_exact_rational_real_imag'][normal]==['1','0'],
            'atlas_hash':text_hash(Path('ford_vertex_equivariance.json').read_bytes())==candidate['atlas_incidence_file_sha256'],
            'partition_hash':text_hash(Path('track_b_global_partition_result.json').read_bytes())==candidate['partition_file_sha256'],
            'global_gates_false':all(witness.get(k) is False for k in ('global_field_constructed','global_projected_mass_certified','global_residual_certified','hejhal_existence_certified'))}
    # Reconstruct the crucial four mode normalization without trusting a flag.
    index={tuple(k[:3]):i for i,k in enumerate(candidate['mode_order'])}
    for k,sign in [((1,0),1),((0,1),-1),((-1,0),1),((0,-1),-1)]:
        coeff=candidate['coefficients_exact_rational_real_imag'][index[('infinity',*k)]]
        checks['unit_orbit:'+str(k)]=coeff==[str(sign),'0']
    return checks


def verify(candidate,witness):
    checks=identity(candidate,witness)
    ctx.prec=256
    r=arb(candidate['r_decimal']); total=arb(0); segments=512
    for j in range(segments):
        # Independent rational grid, four times finer than the producer.
        lo=Fraction(13,10)+Fraction(j,10*segments); hi=lo+Fraction(1,10*segments)
        yy=(arb(lo.numerator)/lo.denominator).union(arb(hi.numerator)/hi.denominator)
        kk=acb(2*arb.pi()*yy).bessel_k(acb(0,r))
        if not kk.is_finite() or not kk.imag.contains(0) or not kk.real.lower()>0:
            checks['finer_grid_positive']=False
            break
        total+=2*kk.real.lower()**2/(10*segments*(arb(hi.numerator)/hi.denominator))
    else:
        checks['finer_grid_positive']=True
    lower=total.lower().sqrt().lower()
    checks['claimed_mass_lower_independently_valid']=bool(lower>=arb(witness['mass_lower']).upper())
    return {'schema':'ford-cusp-mass-independent-verification/v1','verified':all(checks.values()),
            'checks':checks,'arb_bits':256,'segments':segments,'reconstructed_mass_lower':str(lower),
            'conditional_on_same_infinity_cusp_field':True,'global_projected_mass_certified':False,
            'hejhal_existence_certified':False}


if __name__=='__main__':
    c=json.loads(Path('ford_experimental_candidate.json').read_text())
    w=json.loads(Path('ford_cusp_mass_witness.json').read_text())
    out=verify(c,w)
    Path('ford_cusp_mass_independent_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
    raise SystemExit(0 if out['verified'] else 1)
