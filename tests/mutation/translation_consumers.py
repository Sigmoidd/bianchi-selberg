"""Trace primitive translations and mutate executing consumers in memory."""
from collections import Counter
from dataclasses import replace
import difflib
import importlib
import json
from pathlib import Path
import sys
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from flint import arb
from groups import get_group
from groups.data import EllipticClass
from core.assemble import evaluate
from core.certificate import certificate_payload
import groups.registry as registry
import groups.systoles as systoles


def ball(x):
    return {'value':x.str(40),'lower':x.lower().str(40),'upper':x.upper().str(40),'radius':x.rad().str(40)}


def profile_run(fn):
    counts=Counter()
    def profile(frame,event,arg):
        if event!='call':return
        path=Path(frame.f_code.co_filename)
        try: name=str(path.relative_to(ROOT))
        except ValueError:return
        if name.startswith(('groups/','core/','examples/')):
            counts[name+':'+frame.f_code.co_qualname]+=1
    sys.setprofile(profile)
    try: result=fn()
    finally:sys.setprofile(None)
    return result,dict(sorted(counts.items()))


def lists(G,support,truncated):
    domain=range(-1,2) if truncated else range(-6,7)
    all_rows=[];inside=[];ambiguous=[]
    for a in domain:
        for b in domain:
            pair=systoles.trace_pair(G,a,b)
            if pair[0]>9 or (b==0 and abs(a)<=2):continue
            length=systoles.cosh_length(pair).acosh()
            row={'trace':[a,b],'length':ball(length),'pair':[[v.numerator,v.denominator] for v in pair]}
            all_rows.append(row)
            if length.upper()<=support.lower():inside.append(row)
            elif not length.lower()>support.upper():ambiguous.append(row)
    return {'all':all_rows,'up_to_support':inside,'boundary_ambiguous':ambiguous}


def run(outdir, groups=(2, 7)):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    report={'profiles':{},'guard_mutations':[],'coefficient_mutations':[],'length_searches':[]}
    original_coefficient=EllipticClass.coefficient
    original_guard=registry._verify_element_witnesses
    from tests.mutation.run import CASES
    for d in groups:
        name, frac, R = CASES[d]
        def pipeline():
            group=get_group(d)
            E=evaluate(group,frac=frac,R=R,verbose=False)
            return E,certificate_payload(E)
        (base,payload),profile=profile_run(pipeline)
        frozen=json.loads((ROOT/'certificates'/name).read_text())
        if any(payload[k]!=frozen[k] for k in ('bound_ball','lower_endpoint_ball','upper_endpoint_ball','parameters')):
            raise RuntimeError('baseline mismatch')
        inventory=importlib.import_module(f'groups.d{d}_inventory')
        _,arithmetic_profile=profile_run(inventory.verify_inventory)
        report['profiles'][str(d)]={'certificate':profile,'standalone_arithmetic':arithmetic_profile}
        G=base.group
        for i in list(range(len(G.elliptic_classes)))+[None]:
            label='remove all translations' if i is None else 'remove '+G.elliptic_classes[i].label
            targets={C.primitive_translation for j,C in enumerate(G.elliptic_classes) if i is None or i==j}
            # Predict: guard consumer rejects missing witness, before B.
            # Numerical consumer: suppress matching translated class coefficient
            # at its executing consumer, B decreases by its positive NCE share.
            guard_counts={'wrapper_entries':0,'mutated_records_delivered':0,'original_guard_calls':0,'downstream_coefficient_calls':0}
            def guard(group):
                guard_counts['wrapper_entries']+=1
                modified=replace(group,elliptic_classes=tuple(replace(C,primitive_translation=None)
                     if C.primitive_translation in targets else C for C in group.elliptic_classes))
                guard_counts['mutated_records_delivered']+=sum(C.primitive_translation is None for C in modified.elliptic_classes)
                guard_counts['original_guard_calls']+=1
                return original_guard(modified)
            def count_coefficient(C):
                guard_counts['downstream_coefficient_calls']+=1
                return original_coefficient(C)
            E=None;error=None
            with patch.object(registry,'_verify_element_witnesses',guard),patch.object(EllipticClass,'coefficient',count_coefficient):
                try:E=evaluate(G,frac=frac,R=R,verbose=False);certificate_payload(E)
                except (ValueError,ArithmeticError) as exc:error=f'{type(exc).__name__}: {exc}'
            report['guard_mutations'].append({'d':d,'mutation':label,'consumer':'groups.registry._verify_element_witnesses',
                 'baseline_B':ball(base.bound),'mutated_B':None if E is None else ball(E.bound),'delta':None,
                 'check':'arithmetic replay' if error else 'none','error':error,'counts':guard_counts})
            guard_error=error
            numerical_counts={'coefficient_calls':0,'suppressed_translation_contributions':0}
            def suppress(C):
                numerical_counts['coefficient_calls']+=1
                if C.primitive_translation in targets:
                    numerical_counts['suppressed_translation_contributions']+=1
                    return arb(0)
                return original_coefficient(C)
            E=None;error=None;accepted=False
            with patch.object(EllipticClass,'coefficient',suppress):
                try:
                    E=evaluate(G,frac=frac,R=R,verbose=False);certificate_payload(E);accepted=True
                except (ValueError,ArithmeticError) as exc:error=f'{type(exc).__name__}: {exc}'
            row={'d':d,'mutation':label,'consumer':'groups.data.EllipticClass.coefficient',
                 'baseline_B':ball(base.bound),'mutated_B':None if E is None else ball(E.bound),
                 'delta':None if E is None else ball(E.bound-base.bound),'check':'none' if not error else 'arithmetic replay',
                 'error':error,'accepted':accepted,'counts':numerical_counts,
                 'zero_movement':None if E is None else bool(E.bound.mid()==base.bound.mid() and E.bound.rad()==base.bound.rad()),
                 'semantics':'Deletion of the numerical contribution associated with the selected primitive translation(s), not missing-witness removal; original arithmetic records retained.'}
            report['coefficient_mutations'].append(row)
            print('mutation',d,label,'guard',guard_counts,guard_error,'numerical',numerical_counts,'B',None if E is None else str(E.bound),flush=True)
        support=2*base.k*arb(base.delta)
        complete=lists(G,support,False);truncated=lists(G,support,True)
        a=out/f'd{d}-complete-up-to-support.json';b=out/f'd{d}-truncated-up-to-support.json'
        a.write_text(json.dumps(complete['up_to_support'],indent=2)+'\n')
        b.write_text(json.dumps(truncated['up_to_support'],indent=2)+'\n')
        diff=''.join(difflib.unified_diff(a.read_text().splitlines(True),b.read_text().splitlines(True),fromfile=a.name,tofile=b.name))
        (out/f'd{d}-up-to-support.diff').write_text(diff)
        if complete['boundary_ambiguous'] or truncated['boundary_ambiguous']:raise RuntimeError('ambiguous support classification')
        report['length_searches'].append({'d':d,'support_radius':ball(support),'complete':complete,'truncated':truncated,'diff':diff,
                  'meaning':'Trace representatives with multiplicities retained; not a primitive conjugacy-class geodesic inventory.'})
    (out/'translation-consumer-results.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    import argparse
    from tests.mutation.run import CASES
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--groups', type=int, nargs='+', choices=sorted(CASES), default=[2, 7])
    args = parser.parse_args()
    run(args.output, args.groups)
