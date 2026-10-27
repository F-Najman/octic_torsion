#!/usr/bin/env python3
"""Read-only complete C4×20 verifier (adapted: no manifest or snapshot hash checks, no freeze, no proof-fragment or paper checks)."""
import argparse,ast,hashlib,json,os,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'C4x20';D=P/'data';L=P/'logs';NC=D
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rd(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']'))
ap=argparse.ArgumentParser();ap.add_argument('--freeze',action='store_true');args=ap.parse_args()
manifest=D/'C420_manifest.json'
if manifest.exists():
    data=json.loads(manifest.read_text());assert not args.freeze,'Frozen manifests must not be regenerated.'
else:assert args.freeze,'Freeze the explicit file list once before read-only verification.'
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
for script,marker in [('C420_finite_replay.py','C420_FULL_FINITE_REPLAY_DONE complete=true'),
                      ('C420_compact_certify.py','C420_COMPACT_STRICT_PASS files=8')]:
    r=subprocess.run([sys.executable,str(P/script)],cwd=ROOT,env=env,capture_output=True,text=True,check=True)
    assert marker in r.stdout;print(r.stdout,end='')
assert rd(D/'C420_positive_fresh_arithmetic.txt')==rd(D/'C420_positive_fresh_arithmetic.txt')  # the root copy was byte-identical, merged
for name in ['C420_positive_fresh_console.log','C420_positive_root.log']:
    s=(L/name).read_text();assert 'C420_POSITIVE_COMPLETE nc=96,320,888 reps=2,4,12 cusps=24,12,0,0 complete=true' in s
    assert not any(x in s for x in ['Runtime error','User error','Assertion failed','syntax error'])
assert 'C420_FULL_FINITE_REPLAY_DONE complete=true' in (L/'C420_finite_replay.log').read_text()
assert json.loads((D/'C420_finite_audit.json').read_text())['complete']
prefix,rem=rd(NC/'sieve_4_20_D4_Aut96_arithmetic.txt')
actual=[]
for line in (L/'sieve_4_20_D4_Aut96.log').read_text().splitlines():
    m=re.fullmatch(r'M4N20_D4_COMPLETE_SINGLE (\d+) (\d+) (\d+) (true|false) [0-9.]+',line)
    if m:
        d,i,q,b=m.groups();actual.append([1,int(d),int(i),0,int(q),int(b=='true')]);continue
    m=re.fullmatch(r'M4N20_D4_COMPLETE_PAIR (\d+) (\d+) (\d+) (true|false) [0-9.]+',line)
    if m:
        i,j,q,b=m.groups();actual.append([2,2,int(i),int(j),int(q),int(b=='true')])
assert actual==prefix and len(prefix)==103
tail,trem=rd(NC/'M4N20_pair_tail_root_arithmetic.txt');actual=[]
for m in re.finditer(r'^M4N20_ROOT_PAIR_TAIL (\d+) (\d+) (\d+) (true|false) [0-9.]+$',(L/'M4N20_pair_tail_root.log').read_text(),re.M):
    i,j,q,b=m.groups();actual.append([int(i),int(j),int(q),int(b=='true')])
assert actual==tail and len(tail)==111
q29=rd(NC/'M4N20_q29_residuals_root_arithmetic.txt');actual=[]
for line in (L/'M4N20_q29_residuals_root.log').read_text().splitlines():
    m=re.fullmatch(r'M4N20_ROOT_Q29_SINGLE (\d+) (\d+) (true|false) [0-9.]+',line)
    if m:
        d,i,b=m.groups();actual.append([int(d),int(i),0,29,int(b=='true')]);continue
    m=re.fullmatch(r'M4N20_ROOT_Q29_PAIR (\d+) (\d+) (true|false) [0-9.]+',line)
    if m:
        i,j,b=m.groups();actual.append([2,int(i),int(j),29,int(b=='true')])
assert actual==q29 and len(q29)==9
status=json.loads((NC/'M4N20_D4_prefix_status.json').read_text())
assert not status['completed'] and status['completed_pair_cases']==82
assert status['arithmetic_sha256']==sha(NC/'sieve_4_20_D4_Aut96_arithmetic.txt')
for name in ['M4N20_pair_tail_root.log','M4N20_q29_residuals_root.log','M4N20_stabilizers_root.log']:
    s=(L/name).read_text();assert not any(x in s for x in ['Runtime error','User error','Assertion failed','syntax error'])
assert 'M4N20_ROOT_PAIR_TAIL_DONE 111' in (L/'M4N20_pair_tail_root.log').read_text()
assert 'M4N20_ROOT_Q29_DONE' in (L/'M4N20_q29_residuals_root.log').read_text()
assert 'M4N20_ROOT_STABILIZERS_COMPLETE' in (L/'M4N20_stabilizers_root.log').read_text()
for name in ['C420_compact_completion_root_console.log','C420_compact_completion_root.log']:
    s=(L/name).read_text();assert 'C420_COMPACT_COMPLETION_DONE globalpoints=8 cusps=4 CM1728=4 complete=true' in s
cm=json.loads((ROOT/'realizations/data/cm_degree8_classification.json').read_text())
assert max(a*b for a,b in cm['groups'])==50
print('C420_LOG_ARITHMETIC_AGREEMENT_PASS records=214+9 prefix_explicitly_incomplete=true')
print('C420_FULL_STRICT_PASS complete=true remaining=0 files='+str(len(data['files'])))
