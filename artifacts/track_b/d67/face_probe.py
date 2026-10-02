import time,json
from fractions import Fraction as F
from pathlib import Path
from itertools import product
from verify_cover import norm,unit_ideal,times,height

def center(c,l):
 N=norm(c);p=times(l,(c[0]+c[1],-c[1]));return F(p[0],N),F(p[1],N)
rows=[(c,l) for c in product(range(-6,7),range(-1,2)) if 0<norm(c)<=29 and c>(0,0) for l in product(range(-14,15),range(-3,4)) if norm(l)<169 and unit_ideal(c,l)]
print('rows',len(rows),flush=True)
def coeff(c,l):
 u,v=center(c,l);return 2*u+v,u+34*v,F(1,norm(c))-u*u-u*v-17*v*v
Cs=[coeff(*r) for r in rows]
P=((F(-1,2),F(-1,2)),(F(1,2),F(-1,2)),(F(1,2),F(1,2)),(F(-1,2),F(1,2)))
def clip(poly,A,B,C):
 if not poly:return ()
 out=[]
 for s,t in zip(poly,poly[1:]+poly[:1]):
  x=A*s[0]+B*s[1]+C;y=A*t[0]+B*t[1]+C
  if x>=0:out.append(s)
  if (x<0<y) or (y<0<x):
   v=x/(x-y);out.append((s[0]+v*(t[0]-s[0]),s[1]+v*(t[1]-s[1])))
 return tuple(dict.fromkeys(out))
def area(poly):return sum(s[0]*t[1]-s[1]*t[0] for s,t in zip(poly,poly[1:]+poly[:1]))/2
start=time.monotonic();faces=[]
for i,(c,l) in enumerate(rows):
 ai,bi,ci=Cs[i];poly=P
 for j,(aj,bj,cj) in enumerate(Cs):
  if i==j:continue
  poly=clip(poly,ai-aj,bi-bj,ci-cj)
  if len(poly)<3:break
 if len(poly)>=3 and area(poly)>0:faces.append((c,l,poly))
 if i%500==0:print(i,len(faces),time.monotonic()-start,flush=True)
print('faces',len(faces),'area',sum((area(p) for c,l,p in faces),F()),'min vertex',min(height(c,l,*z) for c,l,p in faces for z in p),flush=True)
print('norms',sorted({norm(c) for c,l,p in faces}),flush=True)

Path(__file__).with_name('face_probe_result.json').write_text(json.dumps(dict(schema='d67-face-probe/v1',candidate_rows=len(rows),positive_area_faces=len(faces),projected_area=str(sum((area(p) for c,l,p in faces),F())),lowest_vertex_squared=str(min(height(c,l,*z) for c,l,p in faces for z in p)),active_denominator_norms=sorted({norm(c) for c,l,p in faces}),face_pairings_verified=False,stabilizers_verified=False,geometry_gate_certified=False),indent=2)+'\n')
