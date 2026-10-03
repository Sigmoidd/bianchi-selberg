"""Field-specific fail-closed sensitivity checks; input copies only."""
from copy import deepcopy
from pathlib import Path
import json
import hashlib
from verify_geometry import verify

ROOT=Path(__file__).resolve().parent
base=json.loads((ROOT/'geometry_result.json').read_text())
cover=json.loads((ROOT/'cover_witnesses.json').read_text())
checklist=[]

def negative(label,change_geometry=None,change_cover=None):
    g=deepcopy(base);w=deepcopy(cover)
    if change_geometry:change_geometry(g)
    if change_cover:
        change_cover(w)
        # Give the bad payload a matching hash, so the semantic guard is
        # exercised rather than stopping at provenance validation.
        g['cover_sha256']=hashlib.sha256(
            (json.dumps(w,indent=2)+'\n').encode()).hexdigest()
    try:verify(g,w)
    except (ArithmeticError,ValueError,TypeError,KeyError,IndexError) as e:
        checklist.append(dict(name=label,rejected=True,reason=str(e)))
        return
    raise ArithmeticError(f'NEGATIVE CONTROL WAS ACCEPTED: {label}')

if __name__=='__main__':
    verify(base,cover)
    negative('wrong field',lambda g:g.update(d=67))
    negative('missing lower face',lambda g:g['faces'].pop())
    negative('corrupt pairing matrix',lambda g:g['faces'][3]['matrix'][0].__setitem__(0,99))
    negative('missing vertex stabilizer',lambda g:g['vertex_stabilizers'].pop())
    negative('corrupt edge stabilizer',lambda g:g['edges'][0].update(pointwise_stabilizer=[[[1,0],[0,0],[0,0],[1,0]]]))
    negative('cover gap',change_cover=lambda w:w.pop())
    negative('nonprimitive cover row',change_cover=lambda w:w[0].update(c=[0,0]))
    negative('claimed taller minimum',lambda g:g.update(lowest_vertex_height_squared='1/20'))
    print(json.dumps(dict(geometry_baseline=True,negative_checks=checklist),indent=2))
