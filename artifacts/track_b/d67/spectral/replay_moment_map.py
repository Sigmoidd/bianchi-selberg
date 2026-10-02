"""Bind exact geometry/plan to a fresh independent moment-map replay.

Generated binaries are reproducible intermediates. This verifies the global
face-moment prolongation only; it does not certify threshold-matrix positivity.
"""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from prepare_moment_input import prepare
from reference_mesh import ROOT
from threshold_probe import planar
from verify_adaptive_plan import require, verify_plan
from fractions import Fraction as F


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def replay(directory, rebuild=False):
    require(sys.byteorder == 'little', 'binary replay requires little-endian host')
    coverage = verify_plan()  # includes exact geometry and reference replay
    _, produced = planar(0)
    mesh = json.loads((ROOT / 'reference_mesh.json').read_text())
    stored = [(tuple(map(tuple, (r['c'], r['l']))),
               tuple(tuple(map(F, v)) for v in r['vertices']))
              for r in mesh['triangles']]
    require(produced == stored, 'adaptive initial-root ordering mismatch')
    input_meta = prepare(directory)  # rebuild pairings from exact Ford maps
    prefix = directory / 'full'
    for name in ('moment_map', 'verify_moment_map'):
        subprocess.run(['g++', '-std=c++17', '-O2', '-Wall', '-Wextra',
                        str(ROOT / (name + '.cpp')), '-o', str(directory / name)],
                       check=True)
    if rebuild or not prefix.with_suffix('.topology.bin').exists():
        subprocess.run([str(directory / 'moment_map'),
                        str(directory / 'moment_input.txt'), str(prefix)], check=True)
    subprocess.run([str(directory / 'verify_moment_map'),
                    str(directory / 'moment_input.txt'), str(prefix)], check=True)
    report = json.loads(prefix.with_suffix('.verification.json').read_text())
    require(report['spectral_exclusion_certified'] is False, 'unearned spectral claim')
    require(report['leaves'] == coverage['leaf_tetrahedra'], 'replay leaf count')
    report['source_geometry_binding_requires_wrapper'] = False
    report['geometry_and_pairing_binding_verified'] = True
    report['input_sha256'] = input_meta['input_sha256']
    report['plan_sha256'] = input_meta['plan_sha256']
    report['topology_sha256'] = digest(prefix.with_suffix('.topology.bin'))
    report['rows_sha256'] = digest(prefix.with_suffix('.rows.bin'))
    report['producer_source_sha256'] = digest(ROOT / 'moment_map.cpp')
    report['verifier_source_sha256'] = digest(ROOT / 'verify_moment_map.cpp')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('--rebuild', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = replay(args.directory, args.rebuild)
    data = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(data)
    print(data, end='')
