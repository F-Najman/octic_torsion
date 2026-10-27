#!/usr/bin/env python3
"""Independent integer-lattice replay of all C57 Hecke lift choices."""
import ast
from collections import Counter
from itertools import product
from pathlib import Path
import sys
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
p=Path(sys.argv[1]);raw=p.read_text().replace('[*','[').replace('*]',']').replace('<','(').replace('>',')')
invariants,relations,records=ast.literal_eval(raw);assert invariants==[105,210,0]
h=hermite_normal_form(Matrix(relations).T);assert h.shape==(3,2)
pivots=[]
for j in range(h.cols):
 nz=[i for i in range(h.rows) if h[i,j]];k=max(nz);assert int(h[k,j])>0;pivots.append(k)
assert pivots==[0,1]
def canonical(v):
 v=list(v)
 for j in reversed(range(h.cols)):
  k=pivots[j];q=v[k]//int(h[k,j]);v=[v[i]-q*int(h[i,j]) for i in range(3)]
 return tuple(v)
assert Counter(r[1] for r in records)=={7:2,8:1}
assert {r[0] for r in records}=={1,2,3}
keys=[];total=nonzero=0;bad=[]
for idx,e,key,p,sp,targets,observed in records:
 assert key[0]==e and key not in keys;keys.append(key)
 assert p[2]==sp[2]==e
 assert sum(d*m for d,m,f,c,sc in targets)==6*e
 options=[]
 for d,m,f,c,sc in targets:
  assert f in (1,2) and c[2]==sc[2]==f*d
  if f==1:
   options.append([[ (m-2*a)*(c[j]-sc[j]) for j in range(3)] for a in range(m+1)])
  else:
   assert m%2==0 and c==sc;options.append([[0,0,0]])
 expected=[]
 for choices in product(*options):
  v=[4*(sp[j]-p[j])+sum(q[j] for q in choices) for j in range(3)]
  assert v[2]==0;expected.append(canonical(v))
 assert len(expected)==len(observed)
 assert Counter(expected)==Counter(map(canonical,observed)),(idx,'class multisets differ')
 zeros=expected.count((0,0,0));total+=len(expected);nonzero+=len(expected)-zeros
 if zeros:bad.append({'index':idx,'degree':e,'survivors':zeros})
 # Simultaneous interchange of input and output deck lifts negates every class.
 assert Counter(canonical([-a for a in v]) for v in expected)==Counter(canonical([-a for a in v]) for v in map(canonical,observed))
 print('C57_INTEGER_REPLAY_CASE',{'index':idx,'degree':e,'lift_choices':len(expected),'zero_classes':zeros})
print('C57_GOOD78_INTEGER_REPLAY_DONE',{'source_cases':len(records),'operator_classes':total,'nonzero_classes':nonzero,
      'remaining_cases':bad,'all_split_inert_lift_choices_covered':True,'independent_integer_HNF_replay':True})
