"""Read-only regression diagnostics; run from the repository root."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.request import Request, urlopen
from concurrent.futures import ThreadPoolExecutor

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def replay():
    from core.assemble import evaluate
    from core.certificate import certificate_payload
    ok = True
    for d, filename, frac, R in [(2, 'd2-k2-v2.json', .999, 40), (7, 'd7-k2.json', 1., 256)]:
        frozen = json.loads((ROOT/'certificates'/filename).read_text())
        payload = certificate_payload(evaluate(d, frac=frac, R=R, verbose=False))
        for key in ('parameters', 'bound_ball', 'lower_endpoint_ball', 'upper_endpoint_ball'):
            equal = payload[key] == frozen[key]
            ok &= equal
            print(d, key, 'PASS' if equal else 'FAIL', 'actual:', payload[key], 'frozen:', frozen[key])
    return ok


def cli_outputs(directory):
    ok = True
    for d, filename in [(2, 'd2-k2-v2.json'), (7, 'd7-k2.json')]:
        frozen = json.loads((ROOT/'certificates'/filename).read_text())
        for form in ('script', 'module'):
            result = json.loads((Path(directory)/f'd{d}-{form}.json').read_text())
            equal = result == frozen
            ok &= equal
            print(d, form, 'full frozen JSON equality:', equal)
    return ok


def triage():
    import flint
    from track_b_two_cusp_verify import verify_from_paths
    from track_b_two_cusp_data import canonical_hash, cusp_data
    print('Python:', sys.version, 'python-flint:', flint.__version__, 'FLINT (includes Arb):', flint.__FLINT_VERSION__)
    result = ROOT/'track_b_two_cusp_result.json'
    ledger = ROOT/'track_b_hejhal_rows.jsonl'
    for p in (result, ledger):
        tracked = subprocess.run(['git','ls-files','--error-unmatch',str(p.relative_to(ROOT))], cwd=ROOT, capture_output=True).returncode == 0
        print('file:', p, 'exists:', p.is_file(), 'git tracked:', tracked)
    print(json.dumps(verify_from_paths(result, ledger), indent=2))
    r = json.loads(result.read_text())
    hashes = dict(assembly_hash=canonical_hash(r['assembly_definition']),
                  cusp_data_hash=cusp_data()['cusp_data_hash'],
                  coefficient_vector_hash=canonical_hash(r['physical_coefficient_interval_vector']),
                  collocation_ledger_hash=hashlib.sha256(ledger.read_bytes()).hexdigest())
    for k, v in hashes.items():
        print(k, 'PASS' if v == r[k] else 'FAIL', 'actual:', v, 'stored:', r[k])
    print('artifact cusp data:', json.dumps(r['cusp_data'], indent=2))
    import unittest
    for start in ('.', 'tests'):
        def flatten(suite):
            for test in suite:
                if isinstance(test, unittest.TestSuite):
                    yield from flatten(test)
                else:
                    yield test.id()
        ids = list(flatten(unittest.defaultTestLoader.discover(start)))
        print('discovery:', start, 'count:', len(ids), 'two_cusp:', [i for i in ids if 'two_cusp' in i])
    return True


def protected():
    base = '70c01bb5f0542141fe0fa61451921273be4fbbe5'
    paths = ['certificates/d2-k2.json', 'certificates/d2-k2-v2.json',
             'certificates/legacy-pre-m0.json', 'certificates/m0-regression.json',
             'certificates/mechanical-k2.json', 'old_RIGOR_GAPS.md']
    ok = True
    for name in paths:
        old = subprocess.check_output(['git', 'show', base+':'+name], cwd=ROOT)
        current = (ROOT/name).read_bytes()
        ok &= old == current
        print('PASS' if old == current else 'FAIL', name, hashlib.sha256(current).hexdigest())
    # All pre-existing tracked JSON/JSONL artifacts unchanged by the d=7 commit.
    files = subprocess.check_output(['git','ls-tree','-r','--name-only',base],cwd=ROOT,text=True).splitlines()
    historical = [n for n in files if n.endswith(('.json','.jsonl'))]
    bad = [n for n in historical if subprocess.check_output(['git','show',base+':'+n],cwd=ROOT) != (ROOT/n).read_bytes()]
    print('all pre-existing JSON/JSONL artifacts:', len(historical), 'changed:', bad)
    return ok and not bad


def links(external):
    files = sorted((ROOT/'docs').glob('*.md')) + [ROOT/n for n in ('README.md','README_STRUCTURE.md','RIGOR_GAPS.md','certificates/README.md')]
    bad, urls = [], set()
    count = 0
    for p in files:
        content = re.sub(r'```.*?```', '', p.read_text(), flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)', content):
            target = target.split()[0].strip('<>')
            if re.match(r'\w+://', target):
                urls.add(target)
                continue
            local = target.split('#')[0]
            if not local:
                continue
            count += 1
            if not (p.parent/local).exists():
                bad.append((str(p.relative_to(ROOT)), target))
    print('local path targets:', count, 'broken:', bad)
    print('scope: Markdown inline-link file targets; heading fragments not validated')
    if external:
        def check(url):
            try:
                with urlopen(Request(url, headers={'User-Agent':'Mozilla/5.0'}), timeout=15) as response:
                    return url, response.status, response.url
            except Exception as exc:
                return url, 'UNVERIFIED', str(exc)
        for row in ThreadPoolExecutor(max_workers=8).map(check, sorted(urls)):
            print('external:', *row)
            if row[1] == 'UNVERIFIED':
                bad.append(row)
    else:
        print('external URLs not checked:', len(urls))
    return not bad


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('check', choices=['triage','replay','protected','links','cli'])
    parser.add_argument('--external', action='store_true')
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    if args.check == 'cli' and args.output_dir is None:
        parser.error('cli requires --output-dir')
    success = (links(args.external) if args.check == 'links' else
               cli_outputs(args.output_dir) if args.check == 'cli' else globals()[args.check]())
    raise SystemExit(0 if success else 1)
