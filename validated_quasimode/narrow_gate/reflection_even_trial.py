#!/usr/bin/env python3
"""Exact rational real + t-reflection projection of Picard Q-odd Fourier seed.

The involution T:(x,t,y)->(x,-t,y) normalizes PSL2(Z[i]) and conjugates
Q to Q^-1. Hence Re/T-even projection preserves automorphy and Q oddness.
It is essential for using half-square Fourier Parseval, which does not hold
coefficientwise for a generic 40-orbit Q-odd seed.
"""
from fractions import Fraction as F
import json,hashlib
from pathlib import Path
from level_one_quasimode_analytic import modes_and_orbits
ROOT=Path(__file__).resolve().parent
src=ROOT/'analytic_quasimode_coefficients_50.json'


def real_reflection_coefficients():
    source=json.loads(src.read_text())
    modes,orbits,S=modes_and_orbits(50)
    amap={};mode_index={(a,b):i for i,(a,b,n) in enumerate(modes)}
    for row,o in enumerate(source['orbits']):
        x=F(o['coefficient_real'])
        for j,ix in enumerate(orbits[row]):
            a,b,n=modes[ix]
            assert tuple(o['modes'][j])==(a,b)
            amap[(a,b)]=F(int(S[row,ix]))*x
    assert len(amap)==len(modes)
    avg={(a,b):(amap[(a,b)]+amap[(a,-b)])/2 for a,b,n in modes}
    assert all(avg[(a,b)]==avg[(a,-b)] for a,b,n in modes)
    assert all(avg[(-b,a)]==-avg[(a,b)] for a,b,n in modes)
    assert avg[(1,0)]==1
    result=[]
    for row,orb in enumerate(orbits):
        a,b,n=modes[orb[0]]
        value=avg[(a,b)]*int(S[row,orb[0]])
        assert all(avg[modes[ix][:2]]==int(S[row,ix])*value for ix in orb)
        result.append(value)
    return source,modes,avg,result


def as_decimal(x,decimal_places=12):
    assert (x*10**decimal_places).denominator==1
    z=x*10**decimal_places
    assert z.denominator==1
    v=z.numerator
    return ('-' if v<0 else '')+str(abs(v)//10**decimal_places)+'.'+str(abs(v)%10**decimal_places).zfill(decimal_places)


def main():
    source,modes,avg,coeff=real_reflection_coefficients()
    out={'M':50,'r':source['r'],'lambda':source['lambda'], 'normalization':'a(1,0) = 1',
         'projection':'exact Re and t-even projection of frozen complex 40-orbit seed',
         'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
         'orbits':[{'modes':o['modes'], 'coefficient_real':as_decimal(v), 'coefficient_imag':'0.000000000000'} for o,v in zip(source['orbits'],coeff)]}
    f=ROOT/'analytic_quasimode_coefficients_50_reflection_even.json'
    f.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    assert len(coeff)==40
    print('PASS: exact rational reflection-even and real Q-odd projection')
    print('coefficients:',len(coeff),'original SHA',out['source_sha256'],'projected SHA',hashlib.sha256(f.read_bytes()).hexdigest())
if __name__=='__main__':main()
