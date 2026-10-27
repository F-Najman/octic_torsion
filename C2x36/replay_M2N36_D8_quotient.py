#!/usr/bin/env python3
"""Replay all source coverage and quotient membership with actual relations.

The finite Jacobian divisor-to-class map remains the exact Magma computation.
Its presentation, subgroup arithmetic, residue-degree multipliers, support
coverage, pair additions, and final decision assembly are independently checked.
"""
import ast, contextlib, hashlib, io, json, runpy, time
from collections import Counter, deque
from pathlib import Path
P=Path(__file__).resolve().parent/'data';t0=time.monotonic();bound={}
def read(n):
 raw=(P/n).read_bytes();bound[n]=hashlib.sha256(raw).hexdigest()
 return ast.literal_eval(raw.decode().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
with contextlib.redirect_stdout(io.StringIO()):
 env=runpy.run_path(str(P.parent/'replay_M2N36_D8_direct.py'))
FF,EC,pack=env['FF'],env['EC'],env['pack']
gd,fields,singles,small,pairs,survivors=read('filter_2_36_D8_quotient_arithmetic.txt')
inv,relations,cuspdegrees,hgens,hinv=gd
assert inv==[6,84,252,0] and relations==[[6,0,0,0],[0,84,0,0],[0,0,252,0]]
assert hinv==[2,42,126] and Counter(cuspdegrees)=={1:9,2:3,3:3}
assert all(v[-1]==0 for v in hgens)
def red(v):assert v[-1]==0;return tuple(v[i]%inv[i] for i in range(3))
def add(a,b):return tuple((x+y)%inv[i] for i,(x,y) in enumerate(zip(a,b)))
hh={red([0]*4)};todo=deque(hh);gens=list(map(red,hgens))
while todo:
 a=todo.popleft()
 for g in gens:
  b=add(a,g)
  if b not in hh:hh.add(b);todo.append(b)
assert len(hh)==2*42*126==10584
assert 6*84*252//len(hh)==12
fielddict={d:f for d,f in fields};expected={};own={};Fdict={}
for D in range(3,9):
 a=read(f'M2N36_D8_direct_full{D}_arithmetic.txt')
 f,recs=(a[0],a[7]) if D==8 else (a[1],a[8])
 assert fielddict[D]==f;Fdict[D]=FF(f)
 expected[D]={i:ky for i,(d,ky) in enumerate(recs,1) if d==D};own[D]=[]
def check_image(D,key,e,w):
 F=Fdict[D];b,c=pack(key[:D]),pack(key[D:]);E=EC(F,b,c)
 Q2=E.add((b,b),(b,b));assert Q2 is not None
 targetbc=E.bc((0,0),Q2)
 assert F.degree(targetbc)==e and e*w==D and w in [1,2]
 sourceQ=(b,b);R=None
 for _ in range(18):R=E.add(R,Q2)
 assert R is None
 # Doubling an exact36-point has exact18; its independent2-point persists.
 return targetbc
seen={D:set() for D in range(3,9)};want_survivors=[];totals=Counter()
for row in singles:
 D,i,key,e,w,cl,flag=row
 assert i not in seen[D] and expected[D][i]==key;seen[D].add(i)
 check_image(D,key,e,w);assert flag==int(red(cl) in hh)
 own[D].append(row);totals[D]+=flag
 if flag:want_survivors.append([1,D,i])
for D in range(3,9):assert seen[D]==set(expected[D])
assert [len(own[D]) for D in range(3,9)]==[2,5,30,126,434,2115]
assert [totals[D] for D in range(3,9)]==[0,0,0,18,39,196]
allsmall=read('M2N36_D8_direct_small_all_points.txt');smd={d:rows for d,f,rows in allsmall}
smalls={D:[] for D in [3,4,5]}
for row in small:
 D,i,key,e,w,cl=row;assert i==len(smalls[D])+1
 assert key==smd[D][i-1][0]+smd[D][i-1][1];check_image(D,key,e,w)
 red(cl);smalls[D].append(row)
assert [len(smalls[D]) for D in [3,4,5]]==[48,120,720]
types=[[3,3],[3,4],[3,5],[4,4]];want_indices=[];pc=Counter()
for D,E in types:
 want_indices += [(D,i,E,j) for i in range(1,len(own[D])+1) for j in range(1,len(smalls[E])+1)]
assert len(pairs)==len(want_indices)==2376
for row,idx in zip(pairs,want_indices):
 D,i,E,j,cl,flag=row;assert tuple(row[:4])==idx
 assert red(cl)==add(red(own[D][i-1][5]),red(smalls[E][j-1][5]))
 assert flag==int(red(cl) in hh);pc[D,E]+=flag
 if flag:want_survivors.append([2,D,i,E,j])
assert survivors==want_survivors and len(survivors)==461
assert [pc[tuple(t)] for t in types]==[20,32,84,72]
log=(P.parent/'logs'/'filter_2_36_D8_quotient.log').read_text()
assert 'M2N36_D8_QUOTIENT_FILTER_DONE 2712 2376 461 complete=true' in log
assert 'Runtime error' not in log and 'Assertion failed' not in log
report={'status':'complete necessary-condition filter; not an exclusion',
 'group':[2,36],'prime':5,'genus_source':41,'genus_target':7,
 'map':'(P2,Q36) -> (P2,2Q36)','generic_map_degree':4,
 'actual_Jacobian_relations':relations,'global_target_cusp_image_order':len(hh),
 'target_image_index':12,'all_positive_source_points_verified':True,
 'single_tests':2712,'single_pass_counts':dict(totals),
 'mixed_tests':2376,'mixed_pass_counts':{str(k):v for k,v in pc.items()},
 'all_residue_degree_weights_verified':True,'remaining_test_records':461,
 'verification_boundary':__doc__,'sha256':bound,
 'seconds':round(time.monotonic()-t0,3)}
(P/'replay_M2N36_D8_quotient.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('M2N36_D8_QUOTIENT_REPLAY_DONE')
