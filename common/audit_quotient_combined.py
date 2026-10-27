#!/usr/bin/env python3
"""Exact nonempty-support coverage for the new degree8 quotient sieves."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[1]
n=int(sys.argv[1]);assert n in (46,50)
expected={46:{(4,):5,(5,):4,(6,):5,(7,):24,(8,):71,(4,4):275},50:{(6,):6,(7,):20,(8,):107}}[n]
files=([f'N46_p3_single_completed.log',f'N46_p3_mixed_quotient.log'] if n==46 else ['N50_p3_quotient.log'])
seen={};allcounts=[];details=[]
for fn in files:
 s=(ROOT/f'C{n}'/'logs'/fn).read_text();assert not re.search(r'Runtime error|User error|System error|Assertion failed|syntax error|failed memory request|^ERROR',s,re.I|re.M)
 if 'single_completed' in fn:
  assert s.rstrip().endswith('SINGLE_BRANCH_CHECKPOINT_DONE level=46 p=3 D=8 completed_degrees=4,5,6,7,8')
 else:
  assert f'SUMMARY N={n} M=1 p=3 D=8 shard=0/1 total_survivors=0' in s
  assert s.rstrip().endswith(f'SIEVE_RESUME_DONE FromDegree={8 if n==46 else 6} MaxPart={4 if n==46 else 8}')
 c=list(map(int,re.search(r'NONCUSP counts deg 1\.\.8: \[([^\]]+)\]',s)[1].replace(',',' ').split()));allcounts.append(c)
 for m in re.finditer(r'RESULT type=\[([^\]]+)\] total=(\d+) survivors=(\d+)',s):
  part=tuple(map(int,m[1].replace(',',' ').split()));total=int(m[2]);assert m[3]=='0' and part not in seen
  assert part in expected and total==expected[part];seen[part]=total
  if 'single_completed' in fn:
   ii=[int(m2[1]) for m2 in re.finditer(r'SMART_CASE_RESULT type=\[\s*'+str(part[0])+r'\s*\] index=(\d+) killed=true',s)]
   assert sorted(ii)==list(range(1,total+1))
 rejects=[]
 for m in re.finditer(r'QUOTIENT_REJECT type=\[([^\]]+)\] index=(\d+) base=(\d+)',s):
  part=tuple(map(int,m[1].replace(',',' ').split()));i=int(m[2]);assert part in expected and 1<=i<=expected[part]
  assert int(m[3])==(23 if n==46 else 25);rejects.append((part,i))
 assert len(rejects)==len(set(rejects))
 details.append({'log':fn,'sha256':hashlib.sha256((ROOT/f'C{n}'/'logs'/fn).read_bytes()).hexdigest(),'quotient_rejections':len(rejects)})
assert seen==expected and all(c==allcounts[0] for c in allcounts)
c=allcounts[0]
ind=(ROOT/'common'/'logs'/'independent_counts.log').read_text();m=re.search(rf'DEGREE8_FAST_COUNT m=1 n={n} p=3 .*?closed_place_counts=\[([^\]]+)\]',ind);assert m
assert c==list(map(int,m[1].replace(',',' ').split()))
def parts(n,b=None):
 if n==0:yield ();return
 for a in range(min(n,b or n),0,-1):
  for rest in parts(n-a,a):yield (a,)+rest
possible={p for d in range(1,9) for p in parts(d) if all(c[a-1]>0 for a in p)}
assert possible==set(expected)
report={'degree':8,'level':n,'prime':3,'complete':True,'survivors':0,'support_types':[{'type':list(p),'total':expected[p]} for p in sorted(expected)],'tested_divisors':sum(expected.values()),'closed_place_counts':c,'independent_counts_agree':True,'partition_types_accounted_for':66,'details':details}
(ROOT/f'C{n}'/'data'/f'N{n}_combined_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print('DEGREE8_QUOTIENT_COMBINED_AUDIT_DONE')
