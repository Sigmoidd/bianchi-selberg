"""Reject malformed topology and a changed exact prolongation weight.

Mutations are applied to temporary copies, never the replay inputs.
"""
import argparse
import json
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile


def run(directory):
    rejected=[]
    with tempfile.TemporaryDirectory(dir=directory) as tmp:
        prefix=Path(tmp)/'mutated'
        for name, suffix, offset, content, expected in [
            ('wrong_schema', '.topology.bin', 0, b'BADMAGIC', 'wrong topology schema'),
            ('moved_initial_vertex', '.topology.bin', 48,
             struct.pack('<q', 0), 'initial point mismatch'),
            ('changed_face_weight', '.rows.bin', 16,
             struct.pack('<Q', 17), 'row weight mismatch')]:
            for ext in ('.topology.bin', '.rows.bin'):
                dest=prefix.with_suffix(ext)
                if dest.exists():dest.unlink()
                source=directory/('full'+ext)
                if ext==suffix:shutil.copyfile(source,dest)
                else:dest.symlink_to(source)
            with prefix.with_suffix(suffix).open('r+b') as dest:
                dest.seek(offset);dest.write(content)
            result=subprocess.run([str(directory/'verify_moment_map'),
                str(directory/'moment_input.txt'),str(prefix)],capture_output=True,text=True)
            if result.returncode==0 or expected not in result.stderr:
                raise ArithmeticError(f'{name} was not rejected for {expected}: {result.stderr[-500:]}')
            rejected.append(dict(mutation=name,rejected=True,reason=expected))
    return dict(schema='d67-face-moment-negative-checks/v1',checks=rejected,
                spectral_exclusion_certified=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--output',type=Path);args=p.parse_args()
    data=json.dumps(run(args.directory),indent=2)+'\n'
    if args.output:args.output.write_text(data)
    print(data,end='')
