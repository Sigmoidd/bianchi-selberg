#!/usr/bin/env python3
"""EXPERIMENTAL oversampled physical collocation; never an existence certificate.

Exact rational sample coordinates, repository mode labels and fiber maps are
used. Arb evaluates K_ir; dense constraints, SVDs and scans use midpoints.
All comparisons fix a_infinity,(1,0)=1, with no fitted phase or scale.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import random
import subprocess
import time
from functools import lru_cache

from flint import acb, arb, ctx
import numpy as np
from scipy.linalg import svd
from scipy.optimize import minimize_scalar

from certificate_hashes import text_hash
from global_cusp_symmetry_experiment import sector_matrix
from track_b_collocation_family import exact_h3_action
from track_b_two_cusp_data import GENERATORS, GLUE, canonical_hash, cusp_modes
from track_b_two_cusp_hejhal import ValidatedWhittaker, component_row, midpoint_matrix

ROOT = Path(__file__).resolve().parent
REGIONS = {"picard_region": ("6.55", "6.69"),
           "prior_region": ("6.69", "6.81")}
FIXED_R = ("6.62208", "6.7439020359331625")
BITS = 128


def sample_points(count: int, seed: int) -> list[dict]:
    """Distinct exact interior points of the *actual* curved half-cell."""
    rng = random.Random(seed)
    points, seen = [], set()
    while len(points) < count:
        x = Fraction(rng.randint(-47, 47), 100)
        z = Fraction(rng.randint(1, 47), 100)
        y = Fraction(rng.randint(74, 137), 100)
        key = (x, z, y)
        if x*x+z*z+y*y <= 1 or key in seen or x in (0, z, -z):
            continue
        seen.add(key)
        points.append({"point_id": f"s{seed}-{len(points):03}",
                       "x1": str(x), "x2": str(z), "y": str(y)})
    return points


def mode_keys(cutoff: int) -> list[tuple]:
    mi, mz = cusp_modes(cutoff)
    return [("infinity", a, b) for a, b, _ in mi] + [("zero", a, b) for a, b, _ in mz]


def restriction_indices(smaller: int, larger: int) -> list[int]:
    large = {k: i for i, k in enumerate(mode_keys(larger))}
    return [large[k] for k in mode_keys(smaller)]


def drift(small: np.ndarray, big: np.ndarray, lo: int, hi: int) -> dict:
    idx = restriction_indices(lo, hi)
    d = big[idx]-small
    ni = len(cusp_modes(lo)[0])
    return {"from_M": lo, "to_M": hi, "delta_l2": float(np.linalg.norm(d)),
            "infinity_delta_l2": float(np.linalg.norm(d[:ni])),
            "zero_delta_l2": float(np.linalg.norm(d[ni:])),
            "smaller_coefficient_l2": float(np.linalg.norm(small)),
            "relative_delta": float(np.linalg.norm(d)/max(np.linalg.norm(small), 1e-300)),
            "phase_alignment_used": False, "rescaling_used": False}


class PhysicalOperator:
    """Vectorized reconstruction of the original physical row operator."""
    def __init__(self, points: list[dict], cutoff: int):
        self.points = points
        self.cutoff = cutoff
        self.mi, self.mz = cusp_modes(cutoff)
        self.ni = len(self.mi)
        self.n = self.ni+len(self.mz)
        self.samples = []
        self.rows = []
        index = {}
        for record in points:
            p = tuple(Fraction(record[k]) for k in ("x1", "x2", "y"))
            for name, g in GENERATORS.items():
                target = exact_h3_action(g, p)
                for c in range(6):
                    t = GLUE[name][c]
                    ids = []
                    for copy, point in ((c, p), (t, target)):
                        key = (copy, *point)
                        if key not in index:
                            index[key] = len(self.samples)
                            self.samples.append((copy, point))
                        ids.append(index[key])
                    self.rows.append(tuple(ids))
        self.phases = []
        for copy, p in self.samples:
            modes = self.mi if copy == 0 else self.mz
            u, v = self.frequencies(copy, modes)
            self.phases.append(np.exp(2j*np.pi*(u*float(p[0]+(0 if copy == 0 else copy-1))+v*float(p[1]))))

    @staticmethod
    def frequencies(copy, modes):
        a, b = np.array([(m[0], m[1]) for m in modes], dtype=float).T
        return (a, b) if copy == 0 else ((2*a-b)/5, (a+2*b)/5)

    def assemble(self, r: str | float) -> tuple[np.ndarray, dict]:
        backend = ValidatedWhittaker(arb(str(r)))
        vals = np.zeros((len(self.samples), self.n), dtype=complex)
        shells = {}
        for j, (copy, p) in enumerate(self.samples):
            modes = self.mi if copy == 0 else self.mz
            off = 0 if copy == 0 else self.ni
            radial = []
            for a, b, norm in modes:
                key = (norm, copy == 0, p[2])
                if key not in shells:
                    mag = (arb(norm)/(1 if copy == 0 else 5)).sqrt()
                    y = arb(p[2].numerator)/p[2].denominator
                    val, _ = backend.radial(mag, y)
                    shells[key] = float(val.real.mid())
                radial.append(shells[key])
            vals[j, off:off+len(modes)] = np.array(radial)*self.phases[j]
        ii, jj = np.array(self.rows).T
        return vals[ii]-vals[jj], backend.summary()

    def independently_check(self, r: str, selected: tuple[int, ...] = (0, 1, 13, 27)) -> dict:
        fast, _ = self.assemble(r)
        backend = ValidatedWhittaker(arb(r))
        maxerr = 0.0
        for k in selected:
            i, j = self.rows[k]
            cc, pp = self.samples[i]
            dd, qq = self.samples[j]
            toarb = lambda p: tuple(arb(t.numerator)/t.denominator for t in p)
            left = component_row(cc, toarb(pp), self.mi, self.mz, backend)
            right = component_row(dd, toarb(qq), self.mi, self.mz, backend)
            slow = midpoint_matrix([[a-b for a,b in zip(left,right)]])[0]
            maxerr = max(maxerr, float(np.max(np.abs(slow-fast[k]))))
        return {"checked_row_indices": list(selected), "max_absolute_midpoint_difference": maxerr,
                "agrees_to_1e-16": maxerr < 1e-16, "rank_certified": False}


def column_metrics(a: np.ndarray) -> dict:
    if a.shape[1] == 0:
        return {"unknowns": 0, "numerical_rank": 0, "condition_equilibrated": None,
                "sigma_min_equilibrated": None, "rank_certified": False}
    norms = np.linalg.norm(a, axis=0)
    if np.any(norms == 0):
        return {"unknowns": a.shape[1], "numerical_rank": 0, "condition_equilibrated": None,
                "sigma_min_equilibrated": 0, "rank_certified": False}
    s = svd(a/norms, compute_uv=False, check_finite=True)
    raw = svd(a, compute_uv=False, check_finite=True)
    tol = max(a.shape)*np.finfo(float).eps*s[0]
    return {"unknowns": a.shape[1], "numerical_rank": int(np.sum(s>tol)),
            "condition_equilibrated": float(s[0]/s[-1]),
            "sigma_min_equilibrated": float(s[-1]), "sigma_max_equilibrated": float(s[0]),
            "condition_physical":float(raw[0]/raw[-1]),"sigma_min_physical":float(raw[-1]),
            "rank_certified": False}


def normalize_singular_vector(a: np.ndarray, embedding: np.ndarray, norm_idx: int):
    if a.shape[1] == 0:
        return None
    # embedding columns have one common norm (1 baseline, 2 in odd sectors).
    # Hence this SVD uses a fixed Euclidean metric on full coefficients and
    # is invariant under any orthogonal/unitary rotation of a kernel basis.
    en=np.linalg.norm(embedding,axis=0)
    if not np.allclose(en,en[0],rtol=1e-12,atol=1e-12):
        raise ValueError('coefficient embedding is not uniformly orthogonal')
    _u, _s, vh = svd(a/en[0], full_matrices=False)
    reduced = vh[-1].conj()/en[0]
    full = embedding@reduced
    if abs(full[norm_idx]) < 1e-12*np.linalg.norm(full):
        return None
    # This single prescribed normalization fixes scale and phase for every M.
    return full/full[norm_idx]


def normalized_lstsq(a: np.ndarray, embedding: np.ndarray, norm_idx: int):
    n = embedding[norm_idx]
    if a.shape[1] == 0 or np.linalg.norm(n) < 1e-12:
        return None
    particular = n.conj()/np.vdot(n,n).real
    from scipy.linalg import null_space
    kernel = null_space(n.reshape(1,-1), rcond=1e-12)
    if kernel.shape[1]:
        coeff = particular+kernel@np.linalg.lstsq(a@kernel, -a@particular, rcond=1e-12)[0]
    else:
        coeff = particular
    return embedding@coeff


def coefficient_record(c: np.ndarray) -> list:
    return [[float(v.real),float(v.imag)] for v in c]


def numerical_trial(a, held, embed, norm_idx, c) -> dict:
    out = column_metrics(a@embed)
    if c is None:
        out["normalization_admissible"] = False
        return out
    out.update(normalization_admissible=True, coefficient_vector=coefficient_record(c),
               coefficient_vector_hash=canonical_hash(coefficient_record(c)),
               train_physical_residual_rms=float(np.linalg.norm(a@c)/math.sqrt(a.shape[0])),
               heldout_physical_residual_rms=float(np.linalg.norm(held@c)/math.sqrt(held.shape[0])),
               coefficient_norm=float(np.linalg.norm(c)))
    return out


def jet_matrix(cutoff: int, r: str | float, vertex) -> np.ndarray:
    """24×n values and component-major (dx1,dx2,dt) jets, t=y²."""
    mi, mz = cusp_modes(cutoff)
    n, ni = len(mi)+len(mz), len(mi)
    x, z, t = map(Fraction, vertex)
    ya = (arb(t.numerator)/t.denominator).sqrt()
    y = float(ya.mid())
    backend = ValidatedWhittaker(arb(str(r)))
    jet = np.zeros((24,n),dtype=complex)
    for copy in range(6):
        modes = mi if copy == 0 else mz
        u,v = PhysicalOperator.frequencies(copy,modes)
        off = 0 if copy == 0 else ni
        for k, mode in enumerate(modes):
            mag = (arb(mode[2])/(1 if copy == 0 else 5)).sqrt()
            radial,_ = backend.radial(mag,ya)
            dy,_ = backend.radial_y_derivative(mag,ya)
            ph = np.exp(2j*np.pi*(u[k]*float(x+(0 if copy==0 else copy-1))+v[k]*float(z)))
            val = float(radial.real.mid())*ph
            j = off+k
            jet[copy,j] = val
            jet[6+3*copy,j] = 2j*np.pi*u[k]*val
            jet[6+3*copy+1,j] = 2j*np.pi*v[k]*val
            jet[6+3*copy+2,j] = float(dy.real.mid())*ph/(2*y)
    return jet


@lru_cache(None)
def exact_zero_jet_columns(cutoff: int) -> tuple[int,...]:
    """Prove zero vertex jets in Q[zeta_10], before column equilibration.

Each orbit has a common radial factor. Both its phase sum and its two
frequency-weighted sums must vanish in Q[X]/Phi_10(X).
"""
    import global_quarter_turn as gq
    import ford_vertex_equivariance as geom
    mi,mz=cusp_modes(cutoff); ni=len(mi)
    # X^4=X^3-X^2+X-1 in the tenth cyclotomic field.
    powers=[(Fraction(1),Fraction(0),Fraction(0),Fraction(0))]
    for k in range(1,10):
        a,b,c,d=powers[-1]
        powers.append((-d,a+d,b-d,c+d))
    zeros=[]
    for j,col in enumerate(gq.sector_basis(cutoff,-1)):
        infinity=min(col)<ni
        all_zero=True
        for x,z,t in geom.VERTICES.values():
            for copy in ((0,) if infinity else range(1,6)):
                sums=[[Fraction(0)]*4 for _ in range(3)]
                for idx,sign in col.items():
                    a,b,_=(mi[idx] if infinity else mz[idx-ni])
                    u,v=(Fraction(a),Fraction(b)) if infinity else (Fraction(2*a-b,5),Fraction(a+2*b,5))
                    angle=u*(x+(0 if copy==0 else copy-1))+v*z
                    exponent=10*angle
                    assert exponent.denominator==1
                    root=powers[int(exponent)%10]
                    for k,weight in enumerate((Fraction(sign),sign*u,sign*v)):
                        sums[k]=[aa+weight*bb for aa,bb in zip(sums[k],root)]
                if any(any(row) for row in sums):
                    all_zero=False
        if all_zero:
            zeros.append(j)
    return tuple(zeros)


@lru_cache(None)
def reduced_jet_constraints():
    """Exactly remove all equations implied by parabolic symmetry and oddness.

The nullspace has canonical free coordinates, so the resulting 14 equations
use those coordinates alone. No tiny implied numerical row is equilibrated.
"""
    import ford_vertex_equivariance as geom
    ambient_rows=geom.constraint_rows(parabolic_only=True)
    piv=geom.rref_sparse(ambient_rows)
    free=[j for j in range(144) if j not in piv]
    ns=geom.nullspace(ambient_rows,144)
    full=geom.constraint_rows()
    restricted=[{j:sum((v*ns[j][k] for k,v in row.items()),Fraction(0))
                 for j in range(len(ns))} for row in full]
    rr=geom.rref_sparse(restricted)
    operator=np.zeros((len(rr),144))
    for i,row in enumerate(rr.values()):
        for k,v in row.items():
            operator[i,free[k]]=float(v)
    assert len(free)==20 and len(rr)==14
    return operator,free


def vertex_constraint_embedding(cutoff: int, r: str | float, odd: np.ndarray):
    """Hooked to exact geometric projectors; midpoint nullspace is diagnostic."""
    import ford_vertex_equivariance as geom
    operator,free=reduced_jet_constraints()
    jets=np.vstack([jet_matrix(cutoff,r,point) for point in geom.VERTICES.values()])
    c = operator@jets@odd
    # Identically zero columns cannot be divided by their roundoff norm.
    zero_cols=exact_zero_jet_columns(cutoff)
    c[:,list(zero_cols)]=0
    nonzero=[j for j in range(odd.shape[1]) if j not in zero_cols]
    # Equilibrate constraint columns and rows before numerical rank detection.
    cn = np.linalg.norm(c,axis=0)
    cn[cn==0] = 1
    scaled = c[:,nonzero]/cn[nonzero]
    rn = np.linalg.norm(scaled,axis=1)
    active = rn > 1e-13
    scaled = scaled[active]/rn[active,None]
    _u,s,vh = svd(scaled,full_matrices=True)
    rank = int(np.sum(s > 1e-10*(s[0] if len(s) else 1)))
    kernel = vh[rank:].conj().T/cn[nonzero,None]
    basis = np.zeros((odd.shape[1],kernel.shape[1]+len(zero_cols)),dtype=complex)
    basis[np.array(nonzero),:kernel.shape[1]]=kernel
    for k,j in enumerate(zero_cols):
        basis[j,kernel.shape[1]+k]=1
    # QR affects only coordinates, not the allowed field subspace.
    if basis.shape[1]:
        basis = np.linalg.qr(basis)[0]
    embedding = odd@basis
    return embedding, {"raw_constraint_rows": 144, "exact_independent_geometric_rows":c.shape[0],
                       "exact_parabolic_odd_ambient_jet_dimension":len(free),
                       "exact_zero_jet_column_indices":list(zero_cols),
                       "constraint_rank_exact_upper":min(14,odd.shape[1]-len(zero_cols)), "constraint_rank_numerical": rank,
                       "remaining_unknowns": embedding.shape[1], "constraint_rank_certified": False,
                       "constraint_zero_row_count_numerical": int(np.sum(~active)),
                       "tolerance": 1e-10, "constraint_sigma": [float(v) for v in s],
                       "nullspace_defect": float(np.linalg.norm(c@basis)) if basis.size else 0,
                       "duplicate_constraint_rows_used_as_weights": False}


def symmetry_defect(cutoff: int, coeff: np.ndarray) -> float:
    import global_quarter_turn as gq
    offset = 0
    worst = 0.0
    for modes in cusp_modes(cutoff):
        idx = {(a,b):k for k,(a,b,_) in enumerate(modes)}
        for k,(a,b,_) in enumerate(modes):
            worst = max(worst, abs(coeff[offset+k]+coeff[offset+idx[gq.rot_i((a,b))]]))
        offset += len(modes)
    return float(worst)


def run(cutoffs: list[int], count: int, held_count: int, grid_count: int) -> dict:
    ctx.prec = BITS
    train = sample_points(count,173)
    held = sample_points(held_count,947)
    assert not ({(p['x1'],p['x2'],p['y']) for p in train}&{(p['x1'],p['x2'],p['y']) for p in held})
    out = {"schema":"ford-spectral-stability/v1", "tier":"EXPERIMENTAL",
           "source_commit":subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
           "bits":BITS, "cutoffs":cutoffs, "normalization":"a_infinity,(1,0)=1",
           "physical_row_weighting":"equal physical row weights; RMS divided by sqrt(row count)",
           "spectral_search_objective":"smallest physical singular value on an orthonormal full-coefficient basis; no column equilibration in scan objective",
           "seeds":{"training":173,"heldout":947},"train_points":train,"heldout_points":held,
           "train_points_hash":canonical_hash(train),"heldout_points_hash":canonical_hash(held),
           "grid_count":grid_count,"fixed_parameter_comparison":[],"spectral_scans":[],
           "global_hejhal_defect_certified":False,"hejhal_existence_certified":False,
           "rung4_certified":False,"continuum_certification_attempted":False}
    for M in cutoffs:
        t0=time.time()
        op=PhysicalOperator(train,M); hp=PhysicalOperator(held,M)
        ni=next(i for i,m in enumerate(op.mi) if m[:2]==(1,0))
        odd=sector_matrix(M,-1,op.n)
        for r in FIXED_R:
            a,work=op.assemble(r); h,_=hp.assemble(r)
            adapted,constraints=vertex_constraint_embedding(M,r,odd)
            rec={"M":M,"r":r,"work":work,"constraints":constraints,"systems":{}}
            for name,embed in (("baseline",np.eye(op.n)),("odd",odd),("odd_vertex_jets",adapted)):
                c=normalized_lstsq(a@embed,embed,ni)
                rec['systems'][name]=numerical_trial(a,h,embed,ni,c)
                if c is not None:
                    rec['systems'][name]['odd_symmetry_defect']=symmetry_defect(M,c)
            out['fixed_parameter_comparison'].append(rec)
        for region,(lo,hi) in REGIONS.items():
            cache={}
            def evaluate(r):
                key=float(r)
                if key not in cache:
                    a,work=op.assemble(format(key,'.17g'))
                    adapted,con=vertex_constraint_embedding(M,key,odd)
                    cache[key]=(a,work,adapted,con)
                return cache[key]
            grid=np.linspace(float(lo),float(hi),grid_count)
            for name in ('baseline','odd','odd_vertex_jets'):
                def objective(r):
                    a,_w,ad,_c=evaluate(r)
                    embed=np.eye(op.n) if name=='baseline' else odd if name=='odd' else ad
                    if embed.shape[1]==0 or np.linalg.norm(embed[ni])<1e-12:
                        return math.inf
                    return column_metrics(a@embed)['sigma_min_physical']/np.linalg.norm(embed[:,0])
                values=[objective(r) for r in grid]
                if all(not math.isfinite(v) for v in values):
                    out['spectral_scans'].append({'M':M,'region':region,'system':name,
                        'normalized_trial_exists':False,'reason':'prescribed infinity coefficient normalization unavailable in permitted space'})
                    continue
                candidates=[]
                for k in range(1,len(grid)-1):
                    if values[k]<=values[k-1] and values[k]<=values[k+1]:
                        opt=minimize_scalar(objective,bounds=(grid[k-1],grid[k+1]),method='bounded',
                                            options={'xatol':1e-8,'maxiter':32})
                        candidates.append((float(opt.x),float(opt.fun),False))
                candidates.extend((float(grid[k]),float(values[k]),True) for k in (0,len(grid)-1))
                r,obj,boundary=min(candidates,key=lambda c:c[1])
                a,work,ad,con=evaluate(r)
                embed=np.eye(op.n) if name=='baseline' else odd if name=='odd' else ad
                c=normalize_singular_vector(a@embed,embed,ni)
                h,_=hp.assemble(format(r,'.17g'))
                rec={'M':M,'region':region,'system':name,'r_M':format(r,'.17g'),
                     'boundary_minimum':boundary,'candidate_minima':[
                         {'r':format(x,'.17g'),'sigma_min':v,'boundary':b} for x,v,b in candidates],
                     'grid_r':[format(float(v),'.17g') for v in grid], 'grid_sigma_min':values,
                     'work':work,'constraints':con if name=='odd_vertex_jets' else None,
                     'evaluation_count_shared':len(cache),**numerical_trial(a,h,embed,ni,c)}
                if c is not None:
                    rec['odd_symmetry_defect']=symmetry_defect(M,c)
                out['spectral_scans'].append(rec)
        print(f'M={M} completed in {time.time()-t0:.1f}s',flush=True)
    out['operator_independent_checks']=[PhysicalOperator(train,cutoffs[-1]).independently_check(r) for r in FIXED_R]
    out['cutoff_drift']={}
    for r in FIXED_R:
        for name in ('baseline','odd','odd_vertex_jets'):
            records=[q for q in out['fixed_parameter_comparison'] if q['r']==r]
            out['cutoff_drift'][f'fixed:{r}:{name}']=drift_records(records,lambda q:q['systems'][name])
    for region in REGIONS:
        for name in ('baseline','odd','odd_vertex_jets'):
            records=[q for q in out['spectral_scans'] if q['region']==region and q['system']==name]
            out['cutoff_drift'][f'scan:{region}:{name}']=drift_records(records,lambda q:q)
    # A numerical diagnostic threshold is explicit; it never certifies convergence.
    out['continuum_gate']={'passed':False,'policy':'require last two deltas<0.1, last two r changes<1e-4, '
                           'interior minima, full observable rank, heldout RMS<1e-6 before Arb trial work',
                           'reason':'independent analytic residual and same-trial mass certificates are absent'}
    return out


def drift_records(records, unwrap):
    ds=[]
    for lo,hi in zip(records,records[1:]):
        a,b=unwrap(lo),unwrap(hi)
        if 'coefficient_vector' not in a or 'coefficient_vector' not in b:
            ds.append({'from_M':lo['M'],'to_M':hi['M'],'delta_l2':None,'reason':'no normalized admissible trial'})
            continue
        decode=lambda q:np.array([complex(*v) for v in q['coefficient_vector']])
        d=drift(decode(a),decode(b),lo['M'],hi['M'])
        if 'r_M' in lo and 'r_M' in hi:
            d['r_change']=abs(float(hi['r_M'])-float(lo['r_M']))
        ds.append(d)
    return ds


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cutoffs',default='2,3,4,5,6')
    parser.add_argument('--points',type=int,default=32)
    parser.add_argument('--heldout-points',type=int,default=16)
    parser.add_argument('--grid-count',type=int,default=17)
    parser.add_argument('--json-out',type=Path,default=ROOT/'ford_spectral_stability.json')
    ns=parser.parse_args()
    out=run([int(v) for v in ns.cutoffs.split(',')],ns.points,ns.heldout_points,ns.grid_count)
    out['code_sha256_policy']='sha256-utf8-lf/v1'
    out['code_sha256']={p:text_hash((ROOT/p).read_bytes()) for p in
                        ['ford_equivariant_collocation.py','ford_vertex_equivariance.py','global_quarter_turn.py',
                         'track_b_two_cusp_data.py','track_b_two_cusp_hejhal.py']}
    out['result_definition_hash']=canonical_hash(out)
    ns.json_out.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
