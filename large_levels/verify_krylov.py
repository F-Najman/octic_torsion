#!/usr/bin/env python3
"""Independent modular determinant replay of a char-zero winding certificate."""
import ast,json,re
from pathlib import Path
H=Path(__file__).resolve().parent
s=(H/'data/winding529_krylov_matrix_v2.m').read_text()
p=int(re.search(r'\bp:=(\d+);',s)[1]);d=int(re.search(r'\bd:=(\d+);',s)[1])
v=ast.literal_eval(re.search(r'\bv:=(\[.*?\]);',s,re.S)[1]);a=ast.literal_eval(re.search(r'\bA:=(\[.*?\]);',s,re.S)[1])
assert len(v)==d and len(a)==d*d
A=[a[i*d:(i+1)*d] for i in range(d)]
vs=[];u=v[:]
for k in range(d):
 vs.append(u);u=[sum(u[i]*A[i][j] for i in range(d))%p for j in range(d)]
B=[r[:] for r in vs];rk=0
for j in range(d):
 q=next((i for i in range(rk,d) if B[i][j]%p),None)
 if q is None:continue
 B[rk],B[q]=B[q],B[rk];inv=pow(B[rk][j],-1,p);B[rk]=[(x*inv)%p for x in B[rk]]
 for i in range(rk+1,d):
  c=B[i][j];B[i]=[(x-c*y)%p for x,y in zip(B[i],B[rk])]
 rk+=1
assert rk==d==45
out={'N':529,'character_order':253,'dimension':d,'splitting_prime':p,'Krylov_rank':rk,'verified':True}
(H/'data/winding529_krylov_replay.json').write_text(json.dumps(out,indent=2)+'\n')
print('WINDING529_KRYLOV_INDEPENDENT_REPLAY_DONE',out)
