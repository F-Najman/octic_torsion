#!/usr/bin/env python3
"""Replay portable-record coverage and every logged N28 Hecke decision.

The nonlinear IsPrincipal calculations are Magma's exact certificates;
this independent checker verifies their exhaustive assembly, rather than
claiming to reimplement divisor arithmetic.
"""
import ast, hashlib, json, re, runpy
from pathlib import Path
P = Path(__file__).resolve().parent / 'data'
runpy.run_path(str(P.parent/'replay_N28_D8_fullaut_coverage.py'))
bound = {}
def content(name):
    data = (P.parent/'logs'/name if name.endswith('.log') else P/name).read_bytes()
    bound[name] = hashlib.sha256(data).hexdigest()
    return data.decode()
def read(name):
    return ast.literal_eval(content(name).replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
records, coverage = read('N28_D8_fullaut_Hecke_input.txt')
sigs, pairs, vs, tests, remaining = read('N28_D8_quintic_correction_arithmetic.txt')
assert remaining == [] and len(sigs) == 4 and len(set(map(tuple,sigs))) == 2
assert sorted(sum(pairs, [])) == [1,2,3,4]
wantvs=[]
for pair in pairs:
    assert len(pair) == 2 and sigs[pair[0]-1] == sigs[pair[1]-1]
    for sign in [-1,1]:
        v=[0]*4;v[pair[0]-1]=-sign;v[pair[1]-1]=sign;wantvs.append(v)
assert vs == wantvs
low={1:[],2:[]}
for n,q,corrections,flags in tests:
    want = [[0]*4]
    if q != 19:
        want += [[v[i]+q*w[i] for i in range(4)] for v in vs for w in vs]
    assert corrections == want and len(set(map(tuple,want))) == len(want)
    assert len(flags) == len(want) and all(x in [0,1] for x in flags)
    low[n].append((q,flags))
assert all(low[n] and all(low[n][-1][1][i] == 0 for i in range(len(low[n][-1][1]))) for n in low)
text = content('N28_D8_quintic_cusp_corrections.log')
assert 'N28_D8_FRESH_T5_GENERIC_REFERENCE_EQUAL 30' in text
assert 'N28_D8_QUINTIC_CORRECTIONS_DONE []' in text
high = {}
for shard in [0,1]:
    text = content(f'N28_D8_fullaut_high_Hecke_shard{shard}.log')
    assert 'Runtime error' not in text
    rx=r'N28_D8_FULLAUT_HIGH_HECKE_CASE (\d+) (\d+) (\d+) (\d+) (\d+) (true|false)'
    seen={}
    for s,n,d,i,q,flag in re.findall(rx,text):
        s,n,d,i,q=map(int,[s,n,d,i,q])
        assert s == shard and (n-1)%2 == shard and records[n-1][:2] == [d,i]
        assert d in [6,7,8] and q in [5,11]
        seen.setdefault(n,[]).append((q,flag == 'true'))
    expected={n for n in range(3,50) if (n-1)%2 == shard}
    assert set(seen) == expected
    for n,flags in seen.items():
        assert flags in [[(5,False)],[(5,True),(11,False)],[(5,True),(11,True)]]
    survivors=[n for n,flags in seen.items() if flags[-1][1]]
    summary=re.search(r'N28_D8_FULLAUT_HIGH_HECKE_SURVIVORS \d+ (\[[^\]]*\])',text)
    assert summary and ast.literal_eval(summary[1]) == survivors
    assert re.search(rf'N28_D8_FULLAUT_HIGH_HECKE_DONE shard={shard} tested={len(expected)} remaining={len(survivors)} complete=true',text)
    high.update(seen)
pending=[n for n,flags in sorted(high.items()) if flags[-1][1]]
third=[]
for n in pending:
    name=f'N28_D8_fullaut_third_{n}_arithmetic.txt'
    if not (P/name).exists():
        continue
    rows=read(name)
    assert len(rows)==1 and rows[0][:4]==[n,8,records[n-1][1],19]
    assert isinstance(rows[0][4],bool)
    third += rows
survive=[n for n in pending if not any(r[0]==n and not r[4] for r in third)]
report={'group':[2,28],'degree':8,'records':49,'degrees':[2,5,13,29],
        'quintic_representatives_killed':2,'quintic_corrections_per_q5':17,
        'higher_representatives':47,'two_prime_survivors':pending,
        'third_prime_decisions':third,'unresolved_representatives':survive,
        'all_modular_points_excluded':not survive,'checksums':bound}
(P/'replay_N28_D8_Hecke.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
print('N28_D8_HECKE_ASSEMBLY_REPLAY_DONE')
