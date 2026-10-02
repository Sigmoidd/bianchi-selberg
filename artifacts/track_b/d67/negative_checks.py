"""Fail-closed checks for the d67 cover, without modifying stored artifacts."""
from copy import deepcopy
from pathlib import Path
import json
from verify_cover import verify

HERE=Path(__file__).resolve().parent
original=json.loads((HERE/'cover_witnesses.json').read_text())
checks=[]

def rejected(label,modify):
    payload=deepcopy(original);modify(payload)
    try:verify(payload)
    except (ArithmeticError,KeyError,TypeError,ValueError) as error:
        checks.append(dict(name=label,rejected=True,reason=str(error)))
        return
    raise ArithmeticError('negative control accepted: '+label)

if __name__=='__main__':
    verify(original)
    rejected('missing cover square',lambda p:p.pop())
    rejected('nonprimitive denominator',lambda p:p[0].update(c=[0,0]))
    rejected('shifted witness square',lambda p:p[0].update(a='0'))
    rejected('wrong field witness',lambda p:p[0].update(c=[1,2]))
    print(json.dumps(dict(baseline=True,checks=checks),indent=2))
