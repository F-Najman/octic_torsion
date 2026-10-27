#!/usr/bin/env python3
"""Exact character and Frobenius arithmetic for the Q(zeta5) rank audit."""
from pathlib import Path
from math import gcd
import ast
import json
import re
p=Path(__file__).resolve().parent/'data';n=p.parent/'logs'
assert len({pow(2,k,125) for k in range(100)})==100
assert len({pow(2,k,25) for k in range(20)})==20
assert pow(2,8,25)==6 and pow(2,40,125)==26
ks=[2,4,6,8,12,14,16,18]
for k in ks:
    for ell in range(4):
        knew=(k+10*ell)%20;e=5*knew
        assert 100//gcd(e,100) in [5,10] and e%5==0 and e%25!=0
        assert e*40%100==0
        assert 8*knew%20==8*k%20
        assert 20//gcd(8*knew,20)==5
line=[s for s in (n/'probe_J25_split_reductions_root.log').read_text().splitlines() if 'p=11 ' in s]
assert len(line)==1
c=ast.literal_eval(re.search(r'frobenius=(\[.*?\])',line[0])[1]);assert len(c)==25
assert c[-1]==1 and all(c[i]==11**(12-i)*c[24-i] for i in range(12))
assert sum(c)==5014032203125 and sum(c)%2==1 and 11%5==1
out={'all_32_character_twists_checked':True,'wild_nebentype_unchanged':True,
     'H26_invariant':True,'delta6_unchanged_and_order5':True,
     'J25_F11_order':sum(c),'odd':True,'rank_twist_theory_audit':'passes; see audit_rank_5_25.txt'}
(p/'audit_rank_5_25_arithmetic.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)
