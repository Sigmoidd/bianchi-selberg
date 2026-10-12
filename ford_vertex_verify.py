#!/usr/bin/env python3
"""Independent SymPy reconstruction of the serialized exact Ford jet ledger.

Does not import the producer, its Gaussian arithmetic, its coset action, or
its linear algebra. Enumerates stabilizers afresh from the proved entry bound.
Only exact finite geometry is verified; no analytic certificate is emitted.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import itertools
import json
from pathlib import Path

import sympy as s

POINTS = {'v_00':(0,0,1), 'v_0h':(0,s.Rational(1,2),s.Rational(3,4)),
          'v_mh0':(-s.Rational(1,2),0,s.Rational(3,4)),
          'v_ph0':(s.Rational(1,2),0,s.Rational(3,4)),
          'v_mhh':(-s.Rational(1,2),s.Rational(1,2),s.Rational(1,2)),
          'v_phh':(s.Rational(1,2),s.Rational(1,2),s.Rational(1,2))}
QMAT = s.diag(s.I,1)
T1 = s.Matrix([[1,1],[0,1]])
TI = s.Matrix([[1,s.I],[0,1]])
R = s.diag(s.I,-s.I)
REDUCTIONS = {'v_00':('v_00',s.eye(2)), 'v_0h':('v_mh0',s.eye(2)),
              'v_mh0':('v_0h',R), 'v_ph0':('v_0h',s.eye(2)),
              'v_mhh':('v_mhh',TI), 'v_phh':('v_mhh',s.eye(2))}


def decode(raw):
    return s.Matrix([[s.Integer(v[0])+s.I*s.Integer(v[1]) for v in row] for row in raw])


def projective_key(m):
    def flatten(mm):
        return tuple((int(s.re(v)),int(s.im(v))) for v in mm)
    return min(flatten(m),flatten(-m))


def independent_action_jacobian(m,p):
    x,z,t=s.symbols('x z t', real=True)
    a,b,c,d=list(m)
    zz=x+s.I*z
    den=s.expand((c*zz+d)*s.conjugate(c*zz+d)+c*s.conjugate(c)*t)
    num=s.expand((a*zz+b)*s.conjugate(c*zz+d)+a*s.conjugate(c)*t)
    out=s.Matrix([s.re(num)/den,s.im(num)/den,t/den**2])
    substitutions=dict(zip((x,z,t),p))
    return tuple(s.cancel(v.subs(substitutions)) for v in out), out.jacobian((x,z,t)).subs(substitutions).applyfunc(s.cancel)


def mod_pi(v):
    return (int(s.re(v))+3*int(s.im(v)))%5


def perm(m):
    pairs=[(0,1)]+[(1,k) for k in range(5)]
    index={p:i for i,p in enumerate(pairs)}
    def normalize(a,b):
        return (0,1) if a==0 else (1,(b*pow(a,-1,5))%5)
    out=[]
    for a,b in pairs:
        out.append(index[normalize((a*mod_pi(m[0,0])+b*mod_pi(m[1,0]))%5,
                                   (a*mod_pi(m[0,1])+b*mod_pi(m[1,1]))%5)])
    return out


def qperm():
    # Q conjugation acts on the bottom row by (c,d)->(-ic,d).
    return [0]+[1+(3*k)%5 for k in range(5)]


def inverse(p):
    return [p.index(k) for k in range(6)]


def transport(p,j):
    out=s.zeros(24)
    for c in range(6):
        out[c,p[c]]=1
        out[6+3*c:9+3*c,6+3*p[c]:9+3*p[c]]=j.T
    return out


def independent_stabilizer(p):
    """Independent complex-rational entry enumeration, not a word search."""
    x,y,t=p; z=x+s.I*y
    gis=[s.Integer(a)+s.I*b for a,b in itertools.product(range(-2,3),repeat=2)]
    elements={}
    for c in gis:
        normc=s.expand(c*s.conjugate(c))
        if normc*t>1:
            continue
        rem=1-normc*t
        aa=[a for a in gis if s.expand((a-c*z)*s.conjugate(a-c*z))==rem]
        dd=[d for d in gis if s.expand((c*z+d)*s.conjugate(c*z+d))==rem]
        for a,d in itertools.product(aa,dd):
            b=s.expand_complex(s.cancel((a*d-1)/c)) if c else s.expand(d*z-a*z)
            if not (s.re(b).is_Integer and s.im(b).is_Integer) or s.expand(a*d-b*c-1)!=0:
                continue
            num=s.expand((a*z+b)*s.conjugate(c*z+d)+a*s.conjugate(c)*t)
            if s.expand(num-z)!=0:
                continue
            m=s.Matrix([[a,b],[c,d]])
            elements[projective_key(m)]=m
    return elements


def verify(data):
    checks={}; groups={}; projectors={}
    constraints=[]; structural=[]; names=list(POINTS); offsets={n:24*k for k,n in enumerate(names)}
    def relation(src,tgt,a,sign,struct=False):
        row=s.zeros(24,144)
        row[:,offsets[src]:offsets[src]+24]=s.eye(24)
        row[:,offsets[tgt]:offsets[tgt]+24]+=sign*a
        constraints.append(row)
        if struct:
            structural.append(row)
    for name,p in POINTS.items():
        rec=data['vertices'][name]
        group=independent_stabilizer(p)
        groups[name]=group
        recorded={projective_key(decode(r['matrix'])) for r in rec['stabilizer_elements']}
        checks[name+':complete_stabilizer']=set(group)==recorded and len(group)==rec['stabilizer_order']
        checks[name+':closure']=all(projective_key(a*b) in group for a in group.values() for b in group.values())
        rp=s.zeros(24)
        for row in rec['stabilizer_elements']:
            m=decode(row['matrix']); mapped,j=independent_action_jacobian(m,p)
            checks[name+':element:'+str(projective_key(m))]=(mapped==p and perm(m)==row['forward_perm']
                 and j==s.Matrix(row['jacobian_t']))
            rp+=transport(inverse(perm(m)),j)/len(group)
        projectors[name]=rp
        checks[name+':projector']=rp==s.Matrix(rec['stabilizer_jet_projector']) and rp*rp==rp
        relation(name,name,rp,-1)
        pg=[m for m in group.values() if m[1,0]==0]
        pp=s.zeros(24)
        for m in pg:
            _,j=independent_action_jacobian(m,p)
            pp+=transport(inverse(perm(m)),j)/len(pg)
        # Structural rows are not additional full-system equations.
        row=s.zeros(24,144)
        row[:,offsets[name]:offsets[name]+24]=s.eye(24)-pp
        structural.append(row)
        target,d=REDUCTIONS[name]
        image,j=independent_action_jacobian(d*QMAT,p)
        per=[inverse(perm(d))[k] for k in qperm()]
        qt=transport(per,j)
        checks[name+':Q_reduction']=(image==POINTS[target] and decode(rec['reducing_PSL_matrix'])==d
                 and qt==s.Matrix(rec['Q_odd_jet_transport'])
                 and tuple(map(s.Rational,rec['Q_image_x1_x2_y_squared']))==(-p[1],p[0],p[2]))
        relation(name,target,qt,1,True)
    for source,target in [('v_mh0','v_ph0'),('v_mhh','v_phh')]:
        image,j=independent_action_jacobian(T1,POINTS[source])
        checks[source+':translation_alias']=image==POINTS[target]
        relation(source,target,transport(inverse(perm(T1)),j),-1,True)
    c=s.Matrix.vstack(*constraints)
    rank=c.rank()
    checks['simultaneous_exact_dimension']=144-rank==6
    par=s.Matrix.vstack(*structural)
    par_dimension=144-par.rank()
    checks['structural_jet_dimension20']=par_dimension==20
    checks['genuinely_added_constraint_rank14']=par_dimension-(144-rank)==14
    reduced=s.zeros(14,144)
    for i,row in enumerate(data['reduced_added_constraint_rows_14_by_144']):
        for k,v in row.items():
            reduced[i,int(k)]=s.Rational(v)
    combined=s.Matrix.vstack(par,reduced)
    checks['reduced_constraints_equivalent']=combined.rank()==rank and s.Matrix.vstack(combined,c).rank()==rank
    for component in data['simultaneous_components']:
        cn=component['vertices']; p=s.Matrix(component['projector']); idx=[]
        for name in cn:
            idx.extend(range(offsets[name],offsets[name]+24))
        subset=c[:,idx]
        # Rows coupling outside this component have zero coefficients in it.
        checks['component:'+','.join(cn)]=(subset*p==s.zeros(subset.rows,p.cols) and p*p==p
               and p.rank()==int(component['simultaneous_odd_jet_dimension'])
               and p.rank()==len(idx)-subset.rank())
    for row in data['atlas']['transition_words']:
        m=decode(row['matrix']); qc=QMAT*m*QMAT.inv()
        checks['atlas:'+row['id']]=(m.det()==1 and perm(m)==row['forward_perm']
           and projective_key(qc)==projective_key(decode(row['Q_conjugate_matrix'])) and perm(qc)==row['Q_conjugate_forward_perm'])
    checks['spectral_gates_preserved']=all(data.get(k) is False for k in
        ('spectral_certificate','global_hejhal_defect_certified','rung4_certified','hejhal_existence_certified'))
    return {'schema':'ford-vertex-independent-verification/v1','verified':all(checks.values()),
            'method':'independent SymPy symbolic derivatives, entry enumeration, P1(F5) action and rational ranks',
            'producer_imported':False,'checks':checks,'simultaneous_allowed_jet_dimension':144-rank,
            'structural_jet_dimension':par_dimension,'exact_added_constraint_rank':14,
            'analytic_bounds_certified':False,'hejhal_existence_certified':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path('ford_vertex_equivariance.json'))
    parser.add_argument('--output',type=Path,default=Path('ford_vertex_independent_verification.json'))
    ns=parser.parse_args()
    result=verify(json.loads(ns.input.read_text()))
    ns.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'verified':result['verified'],'checks':len(result['checks']),
                      'failed':[k for k,v in result['checks'].items() if not v]},indent=2))
    return 0 if result['verified'] else 1


if __name__=='__main__': raise SystemExit(main())
