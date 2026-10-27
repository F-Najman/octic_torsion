#!/usr/bin/env python3
"""Third independent GF4 rank implementation using bit-plane columns."""
import ast,hashlib,json
from collections import Counter
from itertools import combinations_with_replacement
from pathlib import Path
P=Path(__file__).resolve().parent

def mul(v,a):
 l,h=v
 return [(0,0),(l,h),(h,l^h),(l^h,l)][a]
def gf4rank(cols):
 basis={}
 for v in cols:
  l,h=v
  while l|h:
   k=(l|h).bit_length()-1;a=((l>>k)&1)+2*((h>>k)&1)
   if k not in basis:
    basis[k]=mul((l,h),[0,1,3,2][a]);break
   x,y=mul(basis[k],a);l^=x;h^=y
 return len(basis)
reports=[]
for name in ['double_cusp_3_21_D4_arithmetic.txt','double_cusp_3_21_D4_dt_arithmetic.txt']:
 f=P/'data'/name;s=f.read_text();N,D,g,ds,keys,As,bad,tested,total=ast.literal_eval(s.replace('[*','[').replace('*]',']'))
 assert (N,D,g)==(21,4,25) and bad==[] and tested==total==17742
 assert Counter(ds)=={1:24,3:8} and len(keys)==len(As)==104
 jets={}
 for (i,m),A in zip(keys,As):
  assert len(A)==25 and all(len(r)==2*m*ds[i-1] for r in A)
  assert all(0<=a<=3 for r in A for a in r)
  jets[i,m]=[(sum((A[r][c]&1)<<r for r in range(25)),sum((A[r][c]>>1)<<r for r in range(25))) for c in range(len(A[0]))]
 rats=[i for i,d in enumerate(ds,1) if d==1];cubs=[i for i,d in enumerate(ds,1) if d==3]
 counts=Counter()
 for ms in combinations_with_replacement(rats,4):
  cols=sum((jets[i,m] for i,m in Counter(ms).items()),[]);assert len(cols)==8
  counts[gf4rank(cols)]+=1
 for i in rats:
  for j in cubs:
   cols=jets[i,1]+jets[j,1];assert len(cols)==8;counts[gf4rank(cols)]+=1
 assert counts=={8:17742}
 reports.append({'arithmetic':name,'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'ranks':dict(counts),'local_blocks':104})
(P/'data'/'root_replay_3_21.json').write_text(json.dumps({'degree':8,'group':[3,21],'independent_GF4_bitplane_replay':True,'datasets':reports},indent=2)+'\n')
print(json.dumps(reports,indent=2));print('ROOT_M3N21_TWO_EXPORTS_BITPLANE_REPLAY_DONE')
