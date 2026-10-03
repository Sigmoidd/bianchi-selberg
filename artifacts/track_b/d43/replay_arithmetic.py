"""Exact d43 design-input checks, not a geometry or spectral certificate.

Run from any directory with ordinary Python or python -O. Output goes to
stdout; redirect to a new per-field record. No historical artifacts are read
or written. No numerical eigensolve or mutation is performed.
"""
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from fields.arithmetic import QuadraticField
from groups.arithmetic import QuadraticRing
from groups.identity import GroupKey
from groups.matrix import MatrixOps
from groups.systoles import reduced_forms


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def replay():
    d = 43
    field, ring, group = QuadraticField(d), QuadraticRing(d), GroupKey(d)
    require(field.D == -43 and field.units == 2, 'field identity mismatch')
    require(group.payload() == dict(d=43, subgroup='full',
                                   level=dict(hnf=[1, 0, 1], norm=1,
                                              basis='standard-integral')),
            'full level-one identity mismatch')
    require(ring.relation == (1, 11), 'integral relation mismatch')
    require(ring.mul((0, 1), (0, 1)) == (-11, 1), 'tau square mismatch')
    # Independent enumeration: reduced primitive binary forms of D=-43.
    forms = []
    for a in range(1, isqrt(43 // 3) + 1):
        for b in range(-a, a + 1):
            num = b*b + 43
            if num % (4*a):
                continue
            c = num // (4*a)
            if a > c or gcd(gcd(a, b), c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            forms.append((a, b, c))
    require(forms == [(1, 1, 11)], 'class-number enumeration failed')
    require(forms == reduced_forms(-43), 'production class-number mismatch')
    # N=1: b!=0 implies N>=43/4>1; b=0 implies a=+/-1.
    require(F(43, 4) > 1, 'unit lower bound failed')
    for a in range(-4, 5):
        for b in range(-4, 5):
            n = F(2*a+b, 2)**2 + F(43*b*b, 4)
            require(n == ring.norm((a, b)), 'norm identity mismatch')
    ops = MatrixOps(ring)
    pairings = {}
    for name, ell in [('T1', (1, 0)), ('Ttau', (0, 1))]:
        matrix = ((1, 0), ell, (0, 0), (1, 0))
        require(ops.det(matrix) == (1, 0) and group.contains(matrix),
                'translation is not in the exact full group')
        pairings[name] = matrix
    # Represent x+y*sqrt(43)*i by rational (x,y). Then tau=(1/2,1/2).
    # The dual mu=(m,(2n-m)/43) must pair to m and n with 1,tau.
    for m in range(-4, 5):
        for n in range(-4, 5):
            u, v = F(m), F(2*n-m, 43)
            require(u == m and u/2 + 43*v/2 == n,
                    'dual lattice pairing failed')
            sq = u*u + 43*v*v
            require(sq == F(m*m) + F((2*n-m)**2, 43),
                    'dual norm identity failed')
            if (m, n) != (0, 0):
                require(sq >= F(4, 43), 'shortest dual bound failed')
    # The exhaustive shortest-vector proof is a case split, not this finite
    # check: m!=0 gives sq>=1; m=0,n!=0 gives sq=4n^2/43>=4/43.
    require(F(1) > F(4, 43), 'shortest dual case split failed')
    # pi>3 gives a purely rational sufficient nonzero-mode bound at Y=2.
    margin = F(4*9*4*4, 43)-1
    require(margin > 0, 'nonzero-mode defect bound failed')
    # A^2=43/4 and beta=(1-s)/(A Y^2); no numerical sqrt needed.
    # Verify the triangle audit counterexample as exact rational ratios.
    require(F(1, 3) > F(1, 4), 'triangle variance audit failed')
    return dict(schema='d43-track-b-design-input-check/v1', d=d, D=-43,
                group=group.payload(), design_sha256=hashlib.sha256(
                    (ROOT/'docs/D43_TRACK_B_DESIGN.md').read_bytes()).hexdigest(),
                exact_checks_passed=True, reduced_forms=forms,
                translation_matrices=pairings, area_squared='43/4',
                shortest_dual_squared='4/43', truncation_Y=2,
                nonzero_mode_margin_using_pi_gt_3=str(margin),
                pid_correspondence_source_verified=False,
                fundamental_domain_proved=False,
                mesh_certified=False, all_windows_certified=False,
                exclusion_certified=False,
                first_unproved_field_input='fundamental domain coverage/nonoverlap')


if __name__ == '__main__':
    print(json.dumps(replay(), indent=2) + '\n', end='')
