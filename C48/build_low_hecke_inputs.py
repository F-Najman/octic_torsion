"""Serialize exactly the admitted low-support source divisors, preserving IDs."""
from pathlib import Path
import hashlib
import runpy

P=Path(__file__).resolve().parent/'data'
src=P/'clean_run/lower_support_arithmetic.py'
a=runpy.run_path(str(src))
s=runpy.run_path(str(P/'source_partition_inputs_complete.py'))
assert a['complete'] is True
out=[]
for row,(degs,indices,vector,admitted) in enumerate(a['tests'],1):
    if not admitted:
        continue
    assert sum(degs)==8 and len(degs)==2
    pairs=[s['representative_pairs'][degs[0]][indices[0]-1],a['all_source_pairs'][degs[1]][indices[1]-1]]
    out.append([row,degs,indices,pairs])
assert len(out)==18
def magma(x):
    return str(x)
text='/* Exact residual inputs from the completed 2076-test low filter. */\n'
text+=f'C48LowArithmeticSHA := "{hashlib.sha256(src.read_bytes()).hexdigest()}";\n'
text+=f'C48LowRowIDs := {magma([z[0] for z in out])};\n'
text+=f'C48LowSupports := {magma([z[1] for z in out])};\n'
text+=f'C48LowIndices := {magma([z[2] for z in out])};\n'
text+=f'C48LowPairs := {magma([z[3] for z in out])};\n'
(P/'low_hecke_inputs.m').write_text(text)
print('C48_LOW_HECKE_INPUTS_EXPORTED',len(out))
