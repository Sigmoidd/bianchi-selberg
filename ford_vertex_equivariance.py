#!/usr/bin/env python3
"""Exact Ford vertex jets and Q transport; no spectral certification.

Coordinates are (x1,x2,t=y^2). A jet is [F_0,...,F_5,
 d_x1 F_0,d_x2 F_0,d_t F_0,...,d_x1 F_5,d_x2 F_5,d_t F_5].
All matrices act on column jets. Fractions throughout; no point sampling.
"""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache
import json
import hashlib
from pathlib import Path
from global_cusp_gaussian import (ONE, ZERO, I, QMAT, TI, T1, R, g, gadd,
 gmul, gconj, gnorm, minv, mmul, mdet)
from global_cusp_character_algebra import forward_perm_from_matrix, quarter_turn_perm, inverse
import track_b_partition_geometry as geo

VERTICES = dict(geo.SPECIAL)
NAMES = tuple(VERTICES)
IDENTITY = (ONE, ZERO, ZERO, ONE)
Q_REDUCTIONS = {
 'v_00': ('v_00', IDENTITY), 'v_0h': ('v_mh0', IDENTITY),
 'v_mh0': ('v_0h', R), 'v_ph0': ('v_0h', IDENTITY),
 'v_mhh': ('v_mhh', TI), 'v_phh': ('v_mhh', IDENTITY),
}
ALIASES = (('v_mh0','v_ph0',T1), ('v_mhh','v_phh',T1))
SOURCE_COMMIT = 'b5f29bd1b314e4d385e5af01805729ba27e3f3e2'
REDUCING_WORDS = {n:('R' if n=='v_mh0' else 'Ti' if n=='v_mhh' else 'I') for n in NAMES}

def evaluate_group_word(word):
 out=geo.ID
 if word=='I': return out
 for letter in word.split(): out=geo.mmul(out,geo.GEN[letter])
 return out

@lru_cache(None)
def stabilizer_words():
 """Recover and independently evaluate words from the committed BFS ledger."""
 source=json.loads(Path(__file__).with_name('track_b_partition_result.json').read_text())
 out={}
 for name,row in source['stabilizers'].items():
  out[name]={}
  for element in row['elements']:
   matrix=geo.canon(tuple(tuple(e) for r in element['matrix'] for e in r))
   word=element['word']
   if evaluate_group_word(word)!=matrix: raise ValueError('stabilizer word/matrix mismatch')
   out[name][matrix]=word
 return out

def source_provenance():
 """Both raw evidence and explicit canonical-LF hashes; no old digest changes."""
 from certificate_hashes import text_hash, TEXT_HASH_POLICY
 files=('track_b_partition_geometry.py','track_b_partition_result.json',
        'track_b_global_partition_arb.py','track_b_global_partition_result.json',
        'track_b_partition_cells.jsonl','track_b_partition_cells_refined.jsonl')
 hashes={}
 for name in files:
  raw=Path(__file__).with_name(name).read_bytes()
  hashes[name]={'raw_sha256':hashlib.sha256(raw).hexdigest(),
                'raw_hash_policy':'sha256-raw-bytes/v1',
                'canonical_sha256':text_hash(raw),'canonical_hash_policy':TEXT_HASH_POLICY}
 return {'source_commit':SOURCE_COMMIT,'artifacts':hashes,'historical_digests_rewritten':False}

def eye(n):
 return [[F(i==j) for j in range(n)] for i in range(n)]

def matmul(a,b):
 bt=list(zip(*b)); return [[sum((x*y for x,y in zip(r,c) if x and y),F(0)) for c in bt] for r in a]

def transpose(a): return [list(r) for r in zip(*a)]
def madd(a,b): return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def scale(a,s): return [[s*x for x in r] for r in a]
def trace(a): return sum(a[k][k] for k in range(len(a)))
def exact_rank(a): return len(rref_sparse([{j:x for j,x in enumerate(r) if x} for r in a]))

def action_t_jacobian(m,p):
 """Poincare action and exact derivative in rational squared-height coordinates."""
 a,b,c,d=m; x,v,t=map(F,p); z=(x,v)
 q=gadd(gmul(c,z),d); azb=gadd(gmul(a,z),b)
 ac=gmul(a,gconj(c)); n=gadd(gmul(azb,gconj(q)),gmul(ac,g(t)))
 den=gnorm(q)+gnorm(c)*t
 if gnorm(mdet(m))!=1: raise ValueError('requires determinant of norm one')
 dc=(c,gmul(I,c)); da=(a,gmul(I,a))
 dd=[2*gmul(e,gconj(q))[0] for e in dc]+[gnorm(c)]
 dn=[gadd(gmul(da[k],gconj(q)),gmul(azb,gconj(dc[k]))) for k in range(2)]+[ac]
 jac=[[(dn[k][j]*den-n[j]*dd[k])/den**2 for k in range(3)] for j in range(2)]
 jac.append([F(k==2)/den**2-2*t*dd[k]/den**3 for k in range(3)])
 return (n[0]/den,n[1]/den,t/den**2),jac

def jet_transport(perm,jac):
 """Pull back a target jet: component c <- perm[c], covector <- jac^T."""
 out=[[F(0)]*24 for _ in range(24)]
 for c in range(6):
  out[c][perm[c]]=F(1)
  for k in range(3):
   for j in range(3): out[6+3*c+k][6+3*perm[c]+j]=jac[j][k]
 return out

@lru_cache(None)
def stabilizer_data(name):
 p=VERTICES[name]; group=sorted(geo.exhaustive_stabilizer(p[:2],p[2])); rows=[]
 for h in group:
  hp,jac=action_t_jacobian(h,p)
  assert hp==p
  perm=forward_perm_from_matrix(h)
  rows.append(dict(matrix=h,forward_perm=perm,inverse_perm=inverse(perm),jacobian=jac,
                   jet=jet_transport(inverse(perm),jac)))
 return rows

@lru_cache(None)
def vertex_projector(name):
 """24x24 exact stabilizer Reynolds projector (before cross-vertex oddness)."""
 rows=stabilizer_data(name); out=[[F(0)]*24 for _ in range(24)]
 for r in rows: out=madd(out,r['jet'])
 return scale(out,F(1,len(rows)))

@lru_cache(None)
def parabolic_vertex_projector(name):
 """Average only the vertex stabilizers with lower-left entry c=0.

 These relations are already satisfied by the exact two-cusp Fourier map.
 """
 rows=[r for r in stabilizer_data(name) if r['matrix'][2]==ZERO]
 out=[[F(0)]*24 for _ in range(24)]
 for r in rows: out=madd(out,r['jet'])
 return scale(out,F(1,len(rows)))

def q_relation(name):
 """Return (target,T) for odd equation J_source + T J_target = 0."""
 target,d=Q_REDUCTIONS[name]; p=VERTICES[name]
 qp,_=action_t_jacobian(QMAT,p); rp,jac=action_t_jacobian(mmul(d,QMAT),p)
 assert rp==VERTICES[target]
 invd=inverse(forward_perm_from_matrix(d)); q=quarter_turn_perm()
 perm=tuple(invd[q[c]] for c in range(6))
 return target,jet_transport(perm,jac)

def alias_relation(source,target,m):
 """Return T for Gamma-equation J_source - T J_target = 0."""
 mp,jac=action_t_jacobian(m,VERTICES[source]); assert mp==VERTICES[target]
 return jet_transport(inverse(forward_perm_from_matrix(m)),jac)

def constraint_rows(names=NAMES,include_odd=True,include_aliases=True,parabolic_only=False):
 """Sparse exact rows C, where C*[J_v for v in names]=0.

 Includes complete stabilizer equations and each Q edge whose endpoints
 occur in names. Output list of sparse dicts (column -> Fraction).
 """
 offsets={n:24*k for k,n in enumerate(names)}; out=[]
 def relation(src,tgt,a,sgn):
  for i in range(24):
   row={offsets[src]+i:F(1)}
   for j,x in enumerate(a[i]):
    if x:
     k=offsets[tgt]+j; row[k]=row.get(k,F(0))+sgn*x
   row={j:x for j,x in row.items() if x}
   if row: out.append(row)
 projector=parabolic_vertex_projector if parabolic_only else vertex_projector
 for name in names: relation(name,name,projector(name),-1)
 if include_odd:
  for name in names:
   target,a=q_relation(name)
   if target in offsets: relation(name,target,a,1)
 if include_aliases:
  for src,tgt,m in ALIASES:
   if src in offsets and tgt in offsets: relation(src,tgt,alias_relation(src,tgt,m),-1)
 return out

def rref_sparse(rows):
 piv={}
 for original in rows:
  row={k:F(v) for k,v in original.items() if v}
  while row:
   j=min(row)
   if j not in piv:
    d=row[j]; piv[j]={k:v/d for k,v in row.items()}; break
   d=row[j]
   for k,v in piv[j].items():
    row[k]=row.get(k,F(0))-d*v
    if not row[k]: del row[k]
 for j in sorted(piv,reverse=True):
  for k in piv:
   if k>=j or j not in piv[k]: continue
   d=piv[k][j]
   for c,v in piv[j].items():
    piv[k][c]=piv[k].get(c,F(0))-d*v
    if not piv[k][c]: del piv[k][c]
 return piv

def nullspace(rows,n):
 piv=rref_sparse(rows); basis=[]
 for j in range(n):
  if j in piv: continue
  v=[F(0)]*n; v[j]=F(1)
  for k,row in piv.items(): v[k]=-row.get(j,F(0))
  basis.append(v)
 return basis

def invert(a):
 n=len(a); rows=[{j:F(v) for j,v in enumerate(r+e) if v} for r,e in zip(a,eye(n))]
 piv=rref_sparse(rows)
 if any(i not in piv for i in range(n)): raise ValueError('singular')
 return [[piv[i].get(n+j,F(0)) for j in range(n)] for i in range(n)]

def jet_metric(name):
 t=VERTICES[name][2]; return [F(1)]*6+[t,t,4*t*t]*6

@lru_cache(None)
def simultaneous_projector(names=NAMES):
 """Metric-orthogonal exact projector onto all Gamma/Q-odd first-jet constraints."""
 ns=nullspace(constraint_rows(names),24*len(names)); n=24*len(names)
 if not ns: return [[F(0)]*n for _ in range(n)]
 b=transpose(ns); weights=sum((jet_metric(name) for name in names),[])
 btw=[[v*w for v,w in zip(row,weights)] for row in ns]
 return matmul(matmul(b,invert(matmul(btw,b))),btw)

@lru_cache(None)
def allowed_odd_vertex_projector(name):
 """24x24 metric-orthogonal projector onto the possible jet at one vertex.

 Cross-vertex compatibility still requires q_relation/alias_relation; the
 local images alone do not encode those relationships.
 """
 names=next(c for c in (('v_00',),('v_0h','v_mh0','v_ph0'),('v_mhh','v_phh')) if name in c)
 k=names.index(name)
 basis=nullspace(constraint_rows(names),24*len(names))
 restricted=[v[24*k:24*(k+1)] for v in basis]
 # RREF selects an independent spanning row basis of the restricted image.
 reduced=rref_sparse([{i:x for i,x in enumerate(v) if x} for v in restricted])
 ns=[[row.get(i,F(0)) for i in range(24)] for row in reduced.values()]
 if not ns: return [[F(0)]*24 for _ in range(24)]
 b=transpose(ns); weights=jet_metric(name)
 btw=[[v*w for v,w in zip(row,weights)] for row in ns]
 return matmul(matmul(b,invert(matmul(btw,b))),btw)

def serial(a):
 if isinstance(a,F): return str(a)
 if isinstance(a,dict): return {str(k):serial(v) for k,v in a.items()}
 if isinstance(a,(tuple,list)): return [serial(x) for x in a]
 return a

def matrix_json(m): return [[[int(v) for v in m[k]] for k in (0,1)],[[int(v) for v in m[k]] for k in (2,3)]]

def fiber_orbits(rows):
 todo=set(range(6)); out=[]
 while todo:
  orbit={rows[k]['forward_perm'][min(todo)] for k in range(len(rows))}
  out.append(sorted(orbit)); todo-=orbit
 return out

def checks_projector(p,weights):
 n=len(p)
 return {'idempotent':matmul(p,p)==p,
         'metric_self_adjoint':all(weights[i]*p[i][j]==weights[j]*p[j][i] for i in range(n) for j in range(n))}

def redundant_fourier_jet_certificate():
 """Universal rank bound for genuinely added vertex rows on an odd trial.

 A Fourier jet map lands in the parabolic/odd/alias nullspace (dimension20).
 Full isotropy cuts this to dimension6. Consequently the additional vertex
 operator has rank at most14 for every cutoff and every spectral parameter.
 """
 rows=constraint_rows(); ns=nullspace(constraint_rows(parabolic_only=True),144)
 restriction=[[sum((x*b[k] for k,x in row.items()),F(0)) for b in ns] for row in rows]
 r=exact_rank(restriction)
 components=[]
 for names in (('v_00',),('v_0h','v_mh0','v_ph0'),('v_mhh','v_phh')):
  par=len(nullspace(constraint_rows(names,parabolic_only=True),24*len(names)))
  full=len(nullspace(constraint_rows(names),24*len(names)))
  components.append({'vertices':names,'parabolic_odd_jet_dimension':par,
                     'full_odd_jet_dimension':full,'additional_rank':par-full})
 return {'parabolic_odd_ambient_jet_dimension':len(ns),'full_odd_jet_dimension':6,
         'exact_additional_rank_on_ambient_jets':r,'universal_added_coefficient_row_rank_upper':r,
         'components':components,'exact_checks_pass':len(ns)==20 and r==14,
         'scope':'upper bound only for the existing exact parabolic-equivariant odd Fourier map'}

@lru_cache(None)
def reduced_vertex_constraint_rows():
 """Exactly 14 independent added rows in the six-vertex 144-jet ordering.

 Rows are selected from full isotropy equations by exact independence AFTER
 restriction to the parabolic/Q/alias jet space. No zero/redundant numerical
 row is equilibrated. On a map satisfying those structural identities these
 rows are equivalent to every full vertex constraint.
 """
 basis=nullspace(constraint_rows(parabolic_only=True),144)
 candidates=constraint_rows(include_odd=False,include_aliases=False)
 chosen=[]; restricted=[]; rank=0
 for row in candidates:
  r={j:sum((x*b[k] for k,x in row.items()),F(0)) for j,b in enumerate(basis)}
  r={j:x for j,x in r.items() if x}
  newrank=len(rref_sparse(restricted+[r]))
  if newrank>rank:
   chosen.append(row); restricted.append(r); rank=newrank
 assert rank==14
 return chosen

def zeta10_power(exponent):
 """Exact Q(zeta_10) coefficients on 1,z,z²,z³; Phi10=z4-z3+z2-z+1."""
 terms=[F(0)]*10; terms[exponent%10]=F(1)
 for k in range(9,3,-1):
  a=terms[k]
  for shift,sign in ((1,1),(2,-1),(3,1),(4,-1)): terms[k-shift]+=sign*a
 return terms[:4]

def structurally_zero_odd_jet_columns(cutoff):
 """Odd-sector column indices whose full jets vanish at every listed vertex.

 This is an exact cyclotomic check of phase and horizontal derivative sums.
 Every orbit shares a radial factor and radial derivative, so these zeros
 hold at every spectral parameter and need no Bessel evaluations.
 """
 from global_quarter_turn import cusp_modes, sector_basis
 mi,mz=cusp_modes(cutoff); cols=sector_basis(cutoff,-1); zero=[]
 for col_index,col in enumerate(cols):
  is_inf=min(col)<len(mi); valid=True
  for x,y,t in VERTICES.values():
   for copy in ((0,) if is_inf else range(1,6)):
    sums=[[F(0)]*4 for _ in range(3)]
    for k,sign in col.items():
     a,b,_=(mi[k] if is_inf else mz[k-len(mi)])
     u,w=(F(a),F(b)) if is_inf else (F(2*a-b,5),F(a+2*b,5))
     exponent=10*(u*(x+(0 if copy==0 else copy-1))+w*y)
     assert exponent.denominator==1
     phase=zeta10_power(int(exponent))
     for j,factor in enumerate((F(1),u,w)):
      sums[j]=[old+sign*factor*z for old,z in zip(sums[j],phase)]
    if any(any(s) for s in sums): valid=False; break
   if not valid: break
  if valid: zero.append(col_index)
 return zero

def cutoff_structural_rank_bounds():
 from global_quarter_turn import sector_basis
 rows=[]
 for cutoff in range(2,7):
  n=len(sector_basis(cutoff,-1)); zeros=structurally_zero_odd_jet_columns(cutoff)
  rows.append({'cutoff':cutoff,'odd_unknowns':n,'identically_zero_jet_columns':zeros,
               'exact_vertex_constraint_rank_upper':min(14,n-len(zeros)),
               'proof':'exact Q(zeta_10) value/horizontal-gradient phase sums; common radial factors'})
 return rows

def build_certificate():
 vertices={}
 for name,p in VERTICES.items():
  rows=stabilizer_data(name); projector=vertex_projector(name)
  target,d=Q_REDUCTIONS[name]; qp,_=action_t_jacobian(QMAT,p); target2,transport=q_relation(name)
  group={geo.canon(r['matrix']) for r in rows}
  vertices[name]={'coordinates_x1_x2_y_squared':p,'Q_image_x1_x2_y_squared':qp,
   'Q_reduced_target':target,'reducing_PSL_matrix':matrix_json(d),
   'reducing_group_word':REDUCING_WORDS[name],
   'reducing_transformation_in_H_pi':d[2]==ZERO,
   'reducing_word_matrix_exact':evaluate_group_word(REDUCING_WORDS[name])==geo.canon(d),
   'reducing_forward_fiber_permutation':forward_perm_from_matrix(d),
   'Q_reduced_fiber_pullback_permutation':tuple(inverse(forward_perm_from_matrix(d))[c] for c in quarter_turn_perm()),
   'Q_odd_jet_transport':transport,'stabilizer_order':len(rows),'fiber_orbits':fiber_orbits(rows),
   'stabilizer_elements':[{'matrix':matrix_json(r['matrix']),'group_word':stabilizer_words()[name][geo.canon(r['matrix'])],
      'group_word_source':'track_b_partition_result.json, exact evaluation rechecked', 'forward_perm':r['forward_perm'],
      'inverse_perm':r['inverse_perm'],'jacobian_t':r['jacobian']} for r in rows],
   'stabilizer_group_closure':all(geo.mmul(a,b) in group for a in group for b in group),
   'stabilizer_value_dimension':trace([r[:6] for r in projector[:6]]),
   'stabilizer_gradient_dimension':trace([r[6:] for r in projector[6:]]),
   'stabilizer_jet_projector':projector,'projector_checks':checks_projector(projector,jet_metric(name))}
 components=[]
 for names in (('v_00',),('v_0h','v_mh0','v_ph0'),('v_mhh','v_phh')):
  p=simultaneous_projector(names); ns=nullspace(constraint_rows(names),24*len(names))
  value_dim=sum(p[24*k+i][24*k+i] for k in range(len(names)) for i in range(6))
  checks=checks_projector(p,sum((jet_metric(name) for name in names),[]))
  residual_zero=all(all(sum(v*p[j][k] for j,v in row.items())==0 for k in range(len(p))) for row in constraint_rows(names))
  components.append({'vertices':names,'simultaneous_odd_value_dimension':value_dim,
    'simultaneous_odd_gradient_dimension':len(ns)-value_dim,'simultaneous_odd_jet_dimension':len(ns),
    'projector':p,'basis_vectors':ns,'projector_checks':checks,'constraint_times_projector_zero':residual_zero})
  for k,name in enumerate(names):
   restricted=[[v[24*k+i] for i in range(24)] for v in ns]
   vertices[name]['allowed_odd_value_dimension']=exact_rank([r[:6] for r in restricted])
   vertices[name]['allowed_odd_gradient_dimension']=exact_rank([r[6:] for r in restricted])
   local=allowed_odd_vertex_projector(name)
   vertices[name]['allowed_odd_jet_projector']=local
   vertices[name]['allowed_odd_projector_checks']=checks_projector(local,jet_metric(name))
   vertices[name]['component_zero_value_and_full_gradient_vanish']=all(not v[24*k] and not any(v[24*k+6:24*k+9]) for v in ns)
 atlas=atlas_ledger()
 redundancy=redundant_fourier_jet_certificate()
 checks=all(v['stabilizer_group_closure'] and all(v['projector_checks'].values()) and all(v['allowed_odd_projector_checks'].values()) for v in vertices.values()) and all(all(c['projector_checks'].values()) and c['constraint_times_projector_zero'] for c in components) and redundancy['exact_checks_pass']
 return serial({'schema':'ford-vertex-equivariance/v1','tier':'PROVED: finite rational algebra; atlas projection theorem conditional on a global input section',
  'provenance':source_provenance(),
  'jet_order':'F[0:6], then (d_x1,d_x2,d_t) per fiber; t=y^2',
  'coordinate_gradient_metric':'diag(t,t,4t^2); d_t=d_y/(2y)',
  'vertices':vertices,'simultaneous_components':components,
  'redundant_fourier_jet_certificate':redundancy,
  'reduced_added_constraint_rows_14_by_144':reduced_vertex_constraint_rows(),
  'cutoff_structural_rank_bounds':cutoff_structural_rank_bounds(),
  'fourier_redundancy':{
   'Q_edge_equations':'implied by the existing odd Fourier sector',
   'T1_alias_equations':'implied by the exact two-cusp Fourier lattices and fiber translations',
   'scalar_fiber_zero_nodal_jets':'implied mode by mode at 0 and both corners',
   'full_stabilizer_constraints':'not implied; inversion-containing isotropy supplies genuine physical rows',
   'corner_counterexample':'zero-cusp beta=(1,0) odd orbit at v_mhh, fiber 1: value=sqrt(5)*radial != identically zero'},
  'translation_aliases':[{'source':s,'target':t,'group_word':'T1','matrix':matrix_json(m),'jet_transport':alias_relation(s,t,m)} for s,t,m in ALIASES],
  'atlas':atlas,'all_exact_checks_pass':checks,'spectral_certificate':False,
  'global_hejhal_defect_certified':False,'rung4_certified':False,'hejhal_existence_certified':False})

def atlas_ledger():
 """Finite incidence/word ledger, explicitly NOT a continuum-overlap cover."""
 incidence=geo.global_singular_strata_incidence()
 faces={'x1m':geo.T1,'x1p':geo.T1I,'x2m':geo.R,'x2p':geo.TIR,'floor':geo.S}
 words={geo.canon(matrix):word for word,matrix in geo.GEN.items()}
 words[geo.ID]='I'
 for entries in stabilizer_words().values():
  for matrix,word in entries.items():
   if matrix not in words or len(word)<len(words[matrix]): words[matrix]=word
 qword={'I':'I','T1':'Ti','T1^-1':'Ti^-1','Ti':'T1^-1','Ti^-1':'T1','R':'R','S':'R S'}
 def record(m):
  qconj=mmul(mmul(QMAT,m),minv(QMAT))
  word=words[geo.canon(m)]
  conjugate_word=' '.join(qword[letter] for letter in word.split())
  assert evaluate_group_word(word)==geo.canon(m)
  assert evaluate_group_word(conjugate_word)==geo.canon(qconj)
  p=forward_perm_from_matrix(m); qp=forward_perm_from_matrix(qconj); q=quarter_turn_perm()
  # rho_q rho(QmQ^-1) = rho(m) rho_q, in the pullback convention.
  intertwining=tuple(qp[q[c]] for c in range(6))==tuple(q[p[c]] for c in range(6))
  assert intertwining
  return {'matrix':matrix_json(m),'group_word':word,'forward_perm':p,
          'inverse_perm':inverse(forward_perm_from_matrix(m)),
          'Q_conjugate_matrix':matrix_json(qconj),'Q_conjugate_forward_perm':qp,
          'Q_conjugate_group_word':conjugate_word,
          'Q_conjugate_intertwining_exact':intertwining}
 rows=[{'id':'face:'+n,'kind':'face',**record(m)} for n,m in faces.items()]
 for stratum in incidence['one_dimensional_strata']:
  for k,raw in enumerate(stratum['elements']):
   m=tuple(tuple(e) for row in raw for e in row)
   rows.append({'id':stratum['id']+':'+str(k),'kind':'singular_stratum_isotropy',**record(m)})
 for name in NAMES:
  for k,r in enumerate(stabilizer_data(name)):
   rows.append({'id':name+':'+str(k),'kind':'vertex_isotropy',**record(r['matrix'])})
 return {'schema':'ford-q-atlas-incidence/v1','method':'global projection P_- W=(W-U_Q W)/2',
  'finite_incidence_complete':incidence['certified'],'incidence':incidence,'transition_words':rows,
  'global_projection_requires_pointwise_partition_equivariance':False,
  'doubled_chart_weight':'chi_j/2 on U_j and (chi_j o Q)/2 on Q^-1 U_j',
  'doubled_partition_sum':'(sum chi_j + sum chi_j o Q)/2=1',
  'original_labels':33,'doubled_labels':66,'each_doubled_weight_factor':'1/2',
  'B0_B1_transfer':'B_double(x)=(B(x)+B(Qx))/2; uniform and L2 bounds preserved',
  'Reynolds_multiplicity':'convex 1/|G| averages; no additional integration leaves or stabilizer-order multiplier',
  'isometric_gradient_laplacian_bounds_preserved':True,
  'existing_partition_derivative_certificate_modified':False,
  'complete_continuum_transition_cover':False,'all_transition_defects_covered':False,
  'odd_trial_mass_certified':False,'odd_trial_residual_certified':False,
  'unresolved':['complete overlap domains for original and Q-transformed charts',
                'mixed original/Q chart overlaps and their value/first-gradient defects',
                'same-trial projected mass and continuum residual']}

def require_complete_transition_coverage(ledger):
 """A finite incidence ledger alone must never enable a continuum consumer."""
 if ledger.get('complete_continuum_transition_cover') is not True or ledger.get('all_transition_defects_covered') is not True:
  raise ValueError('incomplete continuum transition coverage; fail closed')

if __name__=='__main__':
 result=build_certificate()
 Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'all_exact_checks_pass':result['all_exact_checks_pass'],
  'vertices':{n:{k:v[k] for k in ('stabilizer_order','allowed_odd_value_dimension','allowed_odd_gradient_dimension','component_zero_value_and_full_gradient_vanish')} for n,v in result['vertices'].items()}},indent=2))
 raise SystemExit(0 if result['all_exact_checks_pass'] else 1)
