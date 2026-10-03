"""Reachability probes and independent volume series; no production writes."""
import builtins
from contextlib import ExitStack
from dataclasses import replace
import json
from pathlib import Path
import sys
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from flint import arb,ctx
from core.assemble import evaluate
from core.certificate import certificate_payload
from core.bspline import g0,gbspline
from groups import get_group
from groups.data import GroupData,EllipticClass
from groups.registry import InventoryBackend
import groups.systoles as systoles


def ball(x):
    return {'value':x.str(40),'radius':x.rad().str(40),'lower':x.lower().str(40),'upper':x.upper().str(40)}


def independent_volume(G):
    # Independent character tables, not production character_data/zetaK2.
    q=abs(G.field.D)
    residues={8:{1:1,3:1,5:-1,7:-1},7:{1:1,2:1,3:-1,4:1,5:-1,6:-1},
              11:{a:(1 if a in (1,3,4,5,9) else -1) for a in range(1,11)}}[q]
    N=q*(100000//q)
    prefixes=[sum(residues.get(n%q,0) for n in range(1,j+1)) for j in range(q+1)]
    if prefixes[-1] != 0: raise RuntimeError('nonzero period mean')
    finite=sum((arb(residues.get(n%q,0))/n**2 for n in range(1,N+1)),arb(0))
    # Abel summation: prefix_min/(N+1)^2 <= tail <= prefix_max/(N+1)^2.
    tail=(arb(min(prefixes))/(N+1)**2).union(arb(max(prefixes))/(N+1)**2)
    L=finite+tail
    # zeta(2)=pi^2/6, independently from production Hurwitz zeta calls.
    zetaK=arb.pi()**2/6*L
    V=arb(q)*arb(q).sqrt()*zetaK/(4*arb.pi()**2)
    production=G.analytic_data(require_inventory=False)['vol']
    return {'D':G.field.D,'N':N,'character_residues':residues,'prefix_bounds':[min(prefixes),max(prefixes)],
            'L2':ball(L),'zetaK2':ball(zetaK),'independent_volume':ball(V),'production_volume':ball(production),
            'difference':ball(production-V),'overlap':bool(production.overlaps(V)),
            'production_contained_in_independent':bool(V.contains(production)),
            'volume_plus1_percent_overlaps':bool((production*arb('1.01')).overlaps(V)),
            'volume_minus1_percent_overlaps':bool((production*arb('0.99')).overlaps(V))}


def run(output, groups=(2, 7)):
    report={'zero_mutations':[],'inside_support':[],'volume':[]}
    from tests.mutation.run import CASES
    for d in groups:
        _, frac, R = CASES[d]
        G=get_group(d);base=evaluate(G,frac=frac,R=R,verbose=False)
        support=2*base.k*arb(base.delta)
        coefficient=EllipticClass.coefficient
        verify=InventoryBackend.verify
        # Translation mutations already tested: target witness-to-NCE dependency.
        for i in list(range(len(G.elliptic_classes)))+[None]:
            group=replace(G,elliptic_classes=tuple(replace(C,primitive_translation=None)
                       if i is None or i==j else C for j,C in enumerate(G.elliptic_classes)))
            counts={'native_inventory_calls':0,'native_coefficient_calls':0,
                    'diagnostic_coefficient_calls':0,'diagnostic_missing_witness_coefficient_calls':0}
            phase='native'
            def counted_verify(self, group):
                counts['native_inventory_calls']+=1
                return verify(self,group)
            def counted_coefficient(C):
                counts[phase+'_coefficient_calls']+=1
                if phase=='diagnostic' and C.primitive_translation is None:
                    counts['diagnostic_missing_witness_coefficient_calls']+=1
                return coefficient(C)
            error=None
            with patch.object(InventoryBackend,'verify',counted_verify),patch.object(EllipticClass,'coefficient',counted_coefficient):
                try:
                    E=evaluate(group,frac=frac,R=R,verbose=False);certificate_payload(E)
                except (ArithmeticError,ValueError) as exc:
                    error=f'{type(exc).__name__}: {exc}'
                phase='diagnostic'
                with patch.object(GroupData,'require_inventory',lambda self,registry=None: None):
                    E=evaluate(group,frac=frac,R=R,verbose=False)
            target=base.terms['NCE'] if i is None else g0(base.k,base.delta)*coefficient(G.elliptic_classes[i])
            report['zero_mutations'].append({'d':d,'mutation':'translation_all' if i is None else 'translation:'+str(i)+':'+G.elliptic_classes[i].label,
                  'target_term':'NCE' if i is None else 'NCE class contribution', 'baseline_term':ball(target),
                  'classification':3,'counts':counts,'native_error':error,'diagnostic_B':ball(E.bound),
                  'same_B':bool(E.bound.mid()==base.bound.mid() and E.bound.rad()==base.bound.rad()),
                  'reason':'Arithmetic gate reached and rejects; numerical coefficient consumes stored norm, not the missing translation witness. No witness-to-norm recomputation path exists in coefficient().'})
        # Enumerate removed traces with production length formula, independently
        # of monkeypatched loop, to identify the shortest omitted target.
        traces=[];retained_pairs=set()
        for a in range(-6,7):
            for b in range(-6,7):
                pair=systoles.trace_pair(G,a,b)
                if pair[0]>9 or (b==0 and abs(a)<=2): continue
                length=systoles.cosh_length(pair).acosh()
                if -1<=a<=1 and -1<=b<=1: retained_pairs.add(pair)
                else: traces.append(((a,b),pair,length))
        traces.sort(key=lambda row:float(row[2].mid()))
        shortest=traces[0]
        target=next(row for row in traces if row[1] not in retained_pairs)
        counts={'range_replacements':0,'systole_verifications':0}
        original_verify=systoles.verify_systole
        def shrunk_range(*args):
            if args==(-6,7):
                counts['range_replacements']+=1
                return builtins.range(-1,2)
            return builtins.range(*args)
        def counted_systole(group):
            counts['systole_verifications']+=1
            return original_verify(group)
        # LevelOneBackend imports verify_systole directly; patch both bindings.
        with patch.object(systoles,'range',shrunk_range,create=True),patch('core.backends.level_one.verify_systole',counted_systole):
            E=evaluate(G,frac=frac,R=R,verbose=False);certificate_payload(E)
        report['zero_mutations'].append({'d':d,'mutation':'trace_box','target_term':'omitted loxodromic/geodesic sum (not an Evaluation.terms entry)',
             'baseline_term':ball(arb(0)),'baseline_term_status':'zero by proved support, not a stored computed sum',
             'classification':1,'support_radius':ball(support),'systole':ball(G.systole()),
             'shortest_removed_trace':shortest[0],'shortest_removed_length':ball(shortest[2]),
             'removed_length_provably_outside':bool(shortest[2].lower()>support.upper()),
             'inside_probe_target_trace':target[0],'inside_probe_original_length':ball(target[2]),
             'counts':counts,'same_B':bool(E.bound.mid()==base.bound.mid() and E.bound.rad()==base.bound.rad())})
        # Before executing: a target shortened to support/2 should make the
        # complete verifier reject. A truncated verifier cannot visit it.
        inside_length=support/2
        original_cosh=systoles.cosh_length
        for truncated in [False,True]:
            counts={'shortened_length_calls':0,'range_replacements':0}
            def shortened(pair):
                if pair==target[1]:
                    counts['shortened_length_calls']+=1
                    return inside_length.cosh()
                return original_cosh(pair)
            def inside_range(*args):
                if args==(-6,7):
                    counts['range_replacements']+=1
                    return builtins.range(-1,2)
                return builtins.range(*args)
            E=None;error=None;accepted=False
            with ExitStack() as stack:
                stack.enter_context(patch.object(systoles,'cosh_length',shortened))
                if truncated: stack.enter_context(patch.object(systoles,'range',inside_range,create=True))
                try:
                    E=evaluate(G,frac=frac,R=R,verbose=False);certificate_payload(E);accepted=True
                except (ValueError,ArithmeticError) as exc: error=f'{type(exc).__name__}: {exc}'
            report['inside_support'].append({'d':d,'truncated':truncated,'prediction':'reject shorter trace' if not truncated else 'not visited; numerical sum unavailable',
                  'target_trace':target[0],'synthetic_length':ball(inside_length),'support_radius':ball(support),
                  'kernel_at_inside_length':ball(gbspline(inside_length,base.k,base.delta)),
                  'counts':counts,'error':error,'native_export_accepted':accepted,'B':None if E is None else ball(E.bound),
                  'same_B':None if E is None else bool(E.bound.mid()==base.bound.mid() and E.bound.rad()==base.bound.rad()),
                  'caveat':'Synthetic inconsistent trace-length fixture; no actual group geodesic or orbital weight is claimed.'})
        report['volume'].append({'d':d,**independent_volume(G)})
    Path(output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    import argparse
    from tests.mutation.run import CASES
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--groups', type=int, nargs='+', choices=sorted(CASES), default=[2, 7])
    args = parser.parse_args()
    run(args.output, args.groups)
