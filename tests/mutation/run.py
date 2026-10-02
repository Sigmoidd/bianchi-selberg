"""In-memory sensitivity audit. No production or frozen-artifact writes."""
from contextlib import ExitStack
from dataclasses import replace
import builtins
import json
from pathlib import Path
import sys
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from flint import arb
from core.assemble import evaluate
from core.certificate import certificate_payload
from core.backends.level_one import LevelOneBackend
from groups import get_group
from groups.data import GroupData, EllipticClass
import groups.systoles as systoles


def interval(x):
    return dict(ball=x.str(40), lower=x.lower().str(40), upper=x.upper().str(40))


CASES = {2: ('d2-k2-v2.json', .999, 40), 7: ('d7-k2.json', 1., 256),
         11: ('d11-k2.json', 1., 256), 19: ('d19-k2.json', 1., 256)}


def run(output, groups=(2, 7)):
    rows = []
    for d in groups:
        filename, frac, R = CASES[d]
        G = get_group(d)
        baseline = evaluate(G, frac=frac, R=R, verbose=False)
        frozen = json.loads((ROOT/'certificates'/filename).read_text())
        payload = certificate_payload(baseline)
        for key in ('parameters','bound_ball','lower_endpoint_ball','upper_endpoint_ball'):
            if payload[key] != frozen[key]:
                raise RuntimeError(f'BASELINE REPLAY FAILED: d={d} {key}: {payload[key]!r} != {frozen[key]!r}')
        cases = [('drop', i, -1) for i in range(len(G.elliptic_classes))]
        cases += [('centralizer',2,-1),('centralizer',.5,1),
                  ('centralizer_record',2,-1),('centralizer_record',.5,1)]
        cases += [('translation',i,0) for i in range(len(G.elliptic_classes))]
        cases += [('translation_all',None,0),('volume',1.01,1),('volume',.99,-1),
                  ('trace_box',None,0),('cusp',0,1),('cusp',-1,1)]
        for kind, value, prediction in cases:
            label = f'{kind}:{value}'
            if kind in ('drop','translation'):
                label += ':'+G.elliptic_classes[value].label
            group = G
            with ExitStack() as stack:
                if kind == 'drop':
                    group = replace(G, elliptic_classes=tuple(C for i,C in enumerate(G.elliptic_classes) if i != value))
                elif kind in ('translation','translation_all'):
                    group = replace(G, elliptic_classes=tuple(replace(C, primitive_translation=None)
                         if kind == 'translation_all' or i == value else C for i,C in enumerate(G.elliptic_classes)))
                elif kind == 'centralizer_record':
                    group = replace(G, elliptic_classes=tuple(replace(C,
                         finite_centralizer_order=C.finite_centralizer_order*value)
                         for C in G.elliptic_classes))
                elif kind == 'centralizer':
                    original = EllipticClass.coefficient
                    stack.enter_context(patch.object(EllipticClass, 'coefficient', lambda C: original(C)/arb(value)))
                elif kind == 'volume':
                    original = LevelOneBackend.geometry
                    def geometry(self, group):
                        result = original(self, group)
                        return replace(result, volume=result.volume*arb(str(value)))
                    stack.enter_context(patch.object(LevelOneBackend, 'geometry', geometry))
                elif kind == 'trace_box':
                    def shrunk_range(*args):
                        return builtins.range(-1,2) if args == (-6,7) else builtins.range(*args)
                    stack.enter_context(patch.object(systoles, 'range', shrunk_range, create=True))
                elif kind == 'cusp':
                    original = LevelOneBackend.terms
                    def terms(self, *args):
                        result = original(self, *args)
                        for name in ('CE','Ch0','PARg0','PSI','PHIINT'):
                            result[name] *= value
                        return result
                    stack.enter_context(patch.object(LevelOneBackend, 'terms', terms))
                E = None
                error = None
                caught = 'none'
                mode = 'native evaluation'
                export_ok = False
                try:
                    E = evaluate(group, frac=frac, R=R, verbose=False)
                    certificate_payload(E)
                    export_ok = True
                except (ValueError, ArithmeticError) as exc:
                    error = f'{type(exc).__name__}: {exc}'
                    if E is None:
                        caught = 'arithmetic replay' if kind in ('drop','translation','translation_all','centralizer','centralizer_record') else 'completeness check' if kind == 'trace_box' else 'positivity check'
                    else:
                        caught = 'positivity check' if 'B < 1' in str(exc) else 'arithmetic replay'
                if E is None:
                    mode = 'diagnostic evaluation: inventory gate bypassed; no certificate'
                    with patch.object(GroupData, 'require_inventory', lambda self, registry=None: None):
                        E = evaluate(group, frac=frac, R=R, verbose=False)
                delta = E.bound-baseline.bound
                # Marginal interval subtraction is conservative; for identical
                # balls exact equality establishes zero deterministic movement.
                same_ball = (E.bound.mid() == baseline.bound.mid()
                             and E.bound.rad() == baseline.bound.rad())
                midpoint_change = E.bound.mid()-baseline.bound.mid()
                tiny = same_ball or bool(abs(midpoint_change).upper() < baseline.bound.rad())
                opposite = bool(midpoint_change < 0) if prediction > 0 else bool(midpoint_change > 0) if prediction < 0 else not same_ball
                flags = []
                if tiny:
                    flags.append('movement below baseline Arb radius (zero)' if same_ball else 'midpoint movement below baseline Arb radius')
                if opposite:
                    flags.append('opposite prediction' if prediction else 'unexpected nonzero movement')
                if caught == 'none' and E.bound.upper() < 1:
                    flags.append('uncaught and B < 1')
                row = dict(d=d, mutation=label, prediction=prediction, baseline=interval(baseline.bound),
                           baseline_radius=baseline.bound.rad().str(40), mutated=interval(E.bound),
                           delta=interval(arb(0) if same_ball else delta), midpoint_delta=midpoint_change.str(40),
                           check=caught, error=error, interval_source=mode,
                           native_export_accepted=export_ok, flags=flags)
                if kind == 'trace_box':
                    row['truncated_trace_check'] = {k:v for k,v in systoles.verify_systole(group).items() if k not in ('systole',)}
                rows.append(row)
                print(f'd={d} {label}: {E.bound}; check={caught}; flags={flags}', flush=True)
    Path(output).write_text(json.dumps(rows,indent=2)+'\n')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--groups', type=int, nargs='+', choices=sorted(CASES), default=[2, 7])
    args = parser.parse_args()
    run(args.output, args.groups)
