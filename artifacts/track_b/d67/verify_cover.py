"""Independent exact d67 cusp and finite Ford cutoff replay."""
from fractions import Fraction as F
from math import gcd, isqrt
from itertools import product
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
DESIGN=HERE.parents[2]/'docs/D67_TRACK_B_DESIGN.md'

def need(value,why):
    if not value:raise ArithmeticError(why)

def times(x,y):
    a,b=x;c,d=y
    return (a*c-17*b*d,a*d+b*c+b*d)

def norm(x):
    a,b=x;return a*a+a*b+17*b*b

def determinant(x,y):return x[0]*y[1]-x[1]*y[0]

def unit_ideal(c,l):
    v=(c,times(c,(0,1)),l,times(l,(0,1)))
    return gcd(*(abs(determinant(v[i],v[j]))
                 for i in range(4) for j in range(i)))==1

def height(c,l,a,b):
    N=norm(c)
    p=times(l,(c[0]+c[1],-c[1]))
    u,v=F(p[0],N),F(p[1],N)
    x,y=a-u,b-v
    return F(1,N)-x*x-x*y-17*y*y

def reduced_forms():
    out=[]
    for a in range(1,isqrt(67//3)+1):
        for b in range(-a,a+1):
            if (b*b+67)%(4*a):continue
            c=(b*b+67)//(4*a)
            if a>c or gcd(gcd(a,b),c)!=1:continue
            if (abs(b)==a or a==c) and b<0:continue
            out.append((a,b,c))
    return out

def verify(data=None):
    if data is None:data=json.loads((HERE/'cover_witnesses.json').read_text())
    need(reduced_forms()==[(1,1,17)],'class-number enumeration')
    need(len(data)==712,'cover leaf count')
    remaining=set();lowest=F(1)
    for i,row in enumerate(data):
        a,b,w=map(F,(row['a'],row['b'],row['w']))
        c=tuple(row['c']);l=tuple(row['l'])
        need(0<norm(c)<=23 and unit_ideal(c,l),f'primitive row {i}')
        need(w>0 and w.denominator&(w.denominator-1)==0,'non-dyadic leaf')
        corners=((a,b),(a+w,b),(a,b+w),(a+w,b+w))
        v=min(height(c,l,*p) for p in corners)
        need(v>=F(1,34),f'floor bound {i}')
        lowest=min(lowest,v)
        need((a,b,w) not in remaining,'duplicate leaf')
        remaining.add((a,b,w))
    def recurse(a,b,w):
        if (a,b,w) in remaining:
            remaining.remove((a,b,w));return
        need(w>=F(1,512),'uncovered square')
        h=w/2
        for i,j in ((0,0),(1,0),(0,1),(1,1)):
            recurse(a+i*h,b+j*h,h)
    recurse(F(-1,2),F(-1,2),F(1))
    need(not remaining,'overlapping or stray leaves')
    # Classification of all possible denominator norms below 36.
    norms=sorted({norm((a,b)) for a,b in product(range(-7,8),range(-2,3))
                  if 0<norm((a,b))<36})
    need(norms==[1,4,9,16,17,19,23,25,29],'norm gap failed')
    need(F(1,36)<F(1,34),'omitted sphere reaches floor')
    need(F(4,67)<1 and F(64*9,67)>1,'dual mode/Y bound')
    return dict(schema='d67-ford-floor-cover/v1',d=67,D=-67,
                group='PSL2(O_-67),full,level-one',
                design_sha256=hashlib.sha256(DESIGN.read_bytes()).hexdigest(),
                cover_sha256=hashlib.sha256((HERE/'cover_witnesses.json').read_bytes()).hexdigest(),
                leaves=len(data),min_square_witness=str(lowest),
                floor_squared_lower='1/34',candidate_c_norm_max=29,
                candidate_l_norm_strict_upper=169,
                first_omitted_c_norm=36,forms=[list(x) for x in reduced_forms()],
                cover_verified=True,faces_verified=False,stabilizers_verified=False,
                spectral_exclusion_verified=False)

if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
