"""Independent binary replay of all C75 cusp kernels and131790 ranks."""
import ast,re,json
from pathlib import Path
from collections import Counter,deque
P=Path(__file__).resolve().parent
source=(P/'data'/'N75_allbad_orbit_matrices_root.m').read_text()+(P/'data'/'N75_allbad_orbit_kernels_root.m').read_text()
def expand_range(match):
    a,b=int(match[1]),int(match[2])
    step=int(match[3]) if match[3] else 1
    assert step>0
    return repr(list(range(a,b+1,step)))
source=re.sub(r'\[\s*(-?\d+)\s*\.\.\s*(-?\d+)(?:\s+by\s+(-?\d+))?\s*\]',expand_range,source)
tate=json.loads((P/'data'/'audit75_tate_parities.json').read_text())
def value(name):
    matches=re.findall(rf'\b{name}:=(.*?);',source,re.S)
    assert len(matches)==1,name
    return ast.literal_eval(matches[0])

matrices={}
for name,n,m,seq in re.findall(r'(\w+(?:\[<\d+,\d+>\])?):=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',source,re.S):
    n,m=int(n),int(m)
    seq=ast.literal_eval(seq)
    assert len(seq)==n*m and all(x in (0,1) for x in seq) and name not in matrices
    matrices[name]=(m,[sum(seq[i*m+j]<<j for j in range(m)) for i in range(n)])
def get(name):
    return matrices[name][1]
def rank(rows):
    piv={}
    for x in rows:
        while x:
            k=x.bit_length()-1
            if k not in piv:
                piv[k]=x
                break
            x^=piv[k]
    return len(piv)
def multiply(A,B):
    out=[]
    for row in A:
        x=0
        while row:
            bit=row&-row
            x^=B[bit.bit_length()-1]
            row^=bit
        out.append(x)
    return out
def plus(A,B):
    assert len(A)==len(B)
    return [a^b for a,b in zip(A,B)]
def power(A,n):
    out=[1<<i for i in range(len(A))]
    while n:
        if n&1:
            out=multiply(out,A)
        A=multiply(A,A)
        n>>=1
    return out
def space_equal(A,B):
    return rank(A)==rank(B)==rank(A+B)


g=73;I=[1<<i for i in range(g)];C=get('CartierMatrix');M=get('Diamond2Matrix');BM=get('ImageBasis')
assert len(BM)==rank(BM)==71 and rank(C)==69
assert power(M,10)==I and multiply(M,C)==multiply(C,M)
assert space_equal(plus(C,I),BM) and space_equal(multiply(BM,M),BM)
ds=value('ds');widths=value('CuspWidths');perm=value('CuspPermutation2');steps=value('MultiplicitySteps')
assert len(ds)==24 and Counter(zip(ds,widths,steps))=={(1,25,2):10,(1,75,1):10,(4,5,1):2,(4,15,1):2}
assert sorted(perm)==list(range(1,25))
# Match the entire native abstract cusp action to the independent Tate action.
tp={int(k):v for k,v in tate['diamond2_closed_cusp_permutation'].items()};places=tate['target_closed_cusps']
for width in [25,75,5,15]:
 actual=[i for i,w in enumerate(widths) if w==width]
 expected=[i for i,p in enumerate(places) if p['width']==width]
 mapping={actual[0]:expected[0]};todo=deque(actual[:1])
 while todo:
  i=todo.popleft();j=mapping[i];assert (ds[i],steps[i])==(places[j]['degree'],places[j]['source_push_step'])
  ni,nj=perm[i]-1,tp[j]
  if ni in mapping:assert mapping[ni]==nj
  else:mapping[ni]=nj;todo.append(ni)
 assert set(mapping)==set(actual) and set(mapping.values())==set(expected)
orbits=[];left=set(range(1,25))
while left:
 rep=min(left);orb=[];i=rep
 while i not in orb:orb.append(i);i=perm[i-1]
 assert i==rep;left-=set(orb);orbits.append(orb)
assert list(map(len,orbits))==[10,10,2,2]
source_count=local_count=0
for orb in orbits:
 rep=orb[0]
 for m in range(steps[rep-1],8//ds[rep-1]+1,steps[rep-1]):
  SW=get(f'sourcekernels[<{rep},{m}>]');SA=get(f'basisjets[<{rep},{m}>]')
  assert rank(SW)==len(SW) and rank(SW)+rank(SA)==g and multiply(SW,SA)==[0]*len(SW)
  source_count+=1
  for k,i in enumerate(orb):
   DM=power(M,k);Di=power(M,(-k)%10)
   W=get(f'kernelspaces[<{i},{m}>]');A=get(f'basisjets[<{i},{m}>]');J=get(f'localjets[<{i},{m}>]')
   assert space_equal(W,multiply(SW,Di)) and A==multiply(DM,SA) and J==multiply(BM,A)
   assert rank(W)+rank(A)==g and multiply(W,A)==[0]*len(W)
   assert space_equal(SW,multiply(SW,power(M,len(orb))))
   local_count+=1
assert source_count==16 and local_count==128
columns={}
for i in range(1,25):
 for m in range(steps[i-1],8//ds[i-1]+1,steps[i-1]):
  nc,rows=matrices[f'localjets[<{i},{m}>]'];assert nc==m*ds[i-1] and len(rows)==71
  columns[i,m]=[sum(((row>>j)&1)<<r for r,row in enumerate(rows)) for j in range(nc)]
tested=bad=0
def check(i,remaining,cols):
 global tested,bad
 if i==24:
  if remaining==0:
   assert len(cols)==8;tested+=1;bad+=rank(cols)!=8
  return
 for m in range(0,remaining//ds[i]+1,steps[i]):
  check(i+1,remaining-m*ds[i],cols+(columns[i+1,m] if m else []))
check(0,8,[]);assert tested==131790 and bad==0
out={'complete':True,'genus':73,'Cartier_rank':69,'formal_image_dimension':71,'source_kernels':source_count,'transported_kernels':local_count,'admissible_cuspidal_divisors':tested,'exceptions':bad,'Tate_cusp_action_independently_matched':True}
(P/'data'/'C75_allbad_independent_root.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out));print('C75_ALLBAD_INDEPENDENT_REPLAY_PASS')
