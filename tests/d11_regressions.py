"""Read-only d=11 diagnostics; exact frozen replay and independent checks."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from flint import arb, ctx
from core.assemble import evaluate
from core.certificate import certificate_payload
from groups import get_group
from tests.test_d11_inventory import trace_coefficient

BASE = 'bba4a2cd07600c0cb460d9df02ad076ccc5fc5bb'


def replay():
    frozen = json.loads((ROOT/'certificates/d11-k2.json').read_text())
    actual = certificate_payload(evaluate(11, frac=1., R=256, verbose=False))
    ok = True
    for key in ('parameters', 'bound_ball', 'lower_endpoint_ball', 'upper_endpoint_ball'):
        equal = actual[key] == frozen[key]
        print(key, 'PASS' if equal else 'FAIL', actual[key])
        ok &= equal
    return ok


def coefficients():
    ctx.prec = 100
    G = get_group(11)
    production = [C.coefficient() for C in G.elliptic_classes]
    ctx.prec = 200
    independent = [trace_coefficient(C) for C in G.elliptic_classes]
    ok = True
    for C, p, v in zip(G.elliptic_classes, production, independent):
        contained = bool(p.contains(v))
        print(json.dumps(dict(label=C.label, production=p.str(40), radius=p.rad().str(40),
                              matrix_trace_value=v.str(60), delta=(v-p).str(40),
                              contained_in_production=contained)))
        ok &= contained
    return ok


def volume():
    ctx.prec = 100
    # Independent table for chi_{-11}, rather than production character_data.
    residues = {a: (1 if a in (1,3,4,5,9) else -1) for a in range(1,11)}
    N = 11*(100000//11)
    prefixes = [sum(residues.get(n%11,0) for n in range(1,j+1)) for j in range(12)]
    if prefixes[-1] != 0:
        raise ArithmeticError('nonzero period mean')
    finite = sum((arb(residues.get(n%11,0))/n**2 for n in range(1,N+1)), arb(0))
    # Abel summation bounds the tail at a complete period by prefix extrema.
    tail = (arb(min(prefixes))/(N+1)**2).union(arb(max(prefixes))/(N+1)**2)
    L2 = finite + tail
    zetaK2 = arb.pi()**2/6*L2
    independent = arb(11)*arb(11).sqrt()*zetaK2/(4*arb.pi()**2)
    production = get_group(11).analytic_data()['vol']
    contained = bool(independent.contains(production))
    print(json.dumps(dict(D=-11,N=N,character_residues=residues,prefix_bounds=[min(prefixes),max(prefixes)],
                          L2=L2.str(40),zetaK2=zetaK2.str(40),independent_volume=independent.str(40),
                          independent_radius=independent.rad().str(40),production_volume=production.str(40),
                          difference=(production-independent).str(40),production_contained=contained,
                          plus_one_percent_overlaps=bool((production*arb('1.01')).overlaps(independent)),
                          minus_one_percent_overlaps=bool((production*arb('0.99')).overlaps(independent))),indent=2))
    return contained


def protected():
    files = subprocess.check_output(['git','ls-tree','-r','--name-only',BASE],cwd=ROOT,text=True).splitlines()
    names = [n for n in files if n.endswith(('.json','.jsonl')) or n == 'old_RIGOR_GAPS.md']
    changed = [n for n in names if subprocess.check_output(['git','show',BASE+':'+n],cwd=ROOT) != (ROOT/n).read_bytes()]
    print('base',BASE,'protected files',len(names),'changed',changed)
    for n in ('certificates/d2-k2.json','certificates/d2-k2-v2.json','certificates/d7-k2.json'):
        print(n,hashlib.sha256((ROOT/n).read_bytes()).hexdigest())
    return not changed


def cli(directory):
    frozen = json.loads((ROOT/'certificates/d11-k2.json').read_text())
    ok = True
    for form in ('script','module'):
        equal = json.loads((directory/f'd11-{form}.json').read_text()) == frozen
        print('d11',form,'full frozen JSON equality:',equal)
        ok &= equal
    return ok


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('check',choices=['replay','coefficients','volume','protected','cli'])
    parser.add_argument('--output-dir',type=Path)
    args = parser.parse_args()
    if args.check == 'cli' and args.output_dir is None:
        parser.error('cli requires --output-dir')
    ok = cli(args.output_dir) if args.check == 'cli' else globals()[args.check]()
    raise SystemExit(0 if ok else 1)
