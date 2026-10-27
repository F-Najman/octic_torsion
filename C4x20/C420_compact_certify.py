"""Read-only lean completion verifier (adapted: no manifest hash checks, no freeze, no proof-fragment check)."""
import argparse,hashlib,json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'C4x20';D=P/'data';L=P/'logs'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
ap=argparse.ArgumentParser();ap.add_argument('--freeze',action='store_true');args=ap.parse_args()
manifest=D/'C420_compact_manifest.json'
names=['C420_compact_completion_final.m','C420_compact_completion_final.log',
       'C420_compact_completion_final_console.log','C420_compact_completion.tex',
       'C420_compact_replay.py','C420_compact_replay.log','C420_compact_README.txt',
       'C420_compact_certify.py']
if manifest.exists():
    data=json.loads(manifest.read_text())
    assert not args.freeze,'Frozen completion manifests are immutable.'
else:assert args.freeze,'No frozen completion manifest exists.'
console=(L/'C420_compact_completion_final_console.log').read_text()
assert not any(s in console for s in ['Runtime error','User error','Assertion failed','syntax error'])
for s in ['C420_CARTAN_INVARIANT_PASS','C420_CYCLIC5_INVARIANT_PASS',
          'C420_COMPACT_BIRATIONAL_PASS','C420_ELLIPTIC_RANK_TORSION_PASS',
          'C420_COMPACT_COMPLETION_DONE globalpoints=8 cusps=4 CM1728=4 complete=true']:
    assert s in console
points=(L/'C420_compact_completion_final.log').read_text()
assert points.count('C420_GLOBAL_POINT cusp=true')==4
assert points.count('C420_GLOBAL_POINT cusp=false modularj=1728')==4
assert 'C420_COMPACT_COMPLETION_DONE globalpoints=8 cusps=4 CM1728=4 complete=true' in points
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
r=subprocess.run([sys.executable,str(P/'C420_compact_replay.py')],cwd=ROOT,env=env,capture_output=True,text=True,check=True)
assert r.stdout.count('_PASS')==6 and 'C420_COMPACT_REPLAY_DONE complete=true' in r.stdout
print(r.stdout,end='')
assert 'C420_COMPACT_REPLAY_DONE complete=true' in (L/'C420_compact_replay.log').read_text()
print('C420_COMPACT_STRICT_PASS files='+str(len(data['files'])))
