"""Sensitivity checks for the new exact d67 geometry gate."""
from copy import deepcopy
import hashlib,json
from pathlib import Path
from verify_geometry import verify
ROOT=Path(__file__).resolve().parent

def run():
    original=json.loads((ROOT/'geometry_result.json').read_text())
    cover=json.loads((ROOT/'cover_witnesses.json').read_text())
    verify(original,cover);checks=[]
    changes=[('wrong field',lambda p:p.update(d=43)),
             ('wrong level',lambda p:p.update(level_hnf=[2,0,1])),
             ('incomplete denominator cutoff',lambda p:p.update(c_norm_max=23)),
             ('missing face',lambda p:p['faces'].pop()),
             ('nonunimodular face map',lambda p:p['faces'][0]['matrix'][0].__setitem__(0,99)),
             ('missing vertex isotropy',lambda p:next(x for x in p['vertex_stabilizers'] if x['matrices'])['matrices'].pop()),
             ('missing edge',lambda p:p['edges'].pop()),
             ('altered minimum height',lambda p:p.update(lowest_vertex_height_squared='1/34'))]
    for label,mutate in changes:
        payload=deepcopy(original);mutate(payload)
        try:verify(payload,cover)
        except (ArithmeticError,KeyError,ValueError,TypeError) as e:
            checks.append(dict(name=label,rejected=True,reason=str(e)))
        else:raise ArithmeticError('negative control accepted: '+label)
    return dict(baseline=True,checks=checks)

if __name__=='__main__':print(json.dumps(run(),indent=2))
