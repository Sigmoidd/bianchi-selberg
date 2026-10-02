"""Sensitivity checks for the new mapped geometry and claim gates."""
from copy import deepcopy
import json
from reference_mesh import ROOT
from verify_reference_mesh import verify_mesh


def run():
    base=json.loads((ROOT/'reference_mesh.json').read_text())
    verify_mesh(base)
    results=[]
    def reject(name,mutate):
        data=deepcopy(base);mutate(data)
        try:verify_mesh(data,replay_geometry=False)
        except (ArithmeticError,KeyError,IndexError,ValueError) as error:
            results.append(dict(name=name,rejected=True,reason=str(error)));return
        raise ArithmeticError('mutation was accepted: '+name)
    reject('wrong_field',lambda x:x.update(d=43))
    reject('wrong_group',lambda x:x.update(full_level_one=False))
    reject('unearned_spectral_claim',lambda x:x.update(spectral_exclusion_certified=True))
    reject('geometry_hash',lambda x:x.update(geometry_sha256='0'*64))
    reject('missing_floor_triangle',lambda x:x['triangles'].pop())
    reject('wrong_pairing',lambda x:x['triangles'][0].update(paired_triangle=1))
    reject('missing_periodic_edge',lambda x:x['side_pairs'].pop())
    reject('floor_vertex_changed',lambda x:x['triangles'][0]['vertices'][0].__setitem__(0,'-49/100'))
    result=dict(schema='d67-reference-mesh-negative-checks/v1',all_rejected=all(x['rejected'] for x in results),checks=results,spectral_exclusion_certified=False)
    (ROOT/'negative_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':print(json.dumps(run(),indent=2))
