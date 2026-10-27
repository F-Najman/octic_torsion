#!/usr/bin/env python3
"""Read-only mathematical audit; writes only this audit's own JSON."""
from pathlib import Path
from fractions import Fraction
from math import gcd
import hashlib,json,re
H=Path(__file__).resolve().parent
P=H.parent
cache=json.loads((P/'common/data/lmfdb_newforms_cache.json').read_text())
def ds(n):return [d for d in range(1,n+1) if n%d==0]
def phi(n):return sum(gcd(a,n)==1 for a in range(1,n+1))
out={}
for n,g in [(45,41),(55,81),(57,85),(63,97),(65,121),(75,145)]:
    forms=[];total=0
    for m in ds(n):
        for f in cache[str(m)]:
            assert f['analytic_rank_proved'] is True
            assert f['analytic_rank'] is not None
            mult=len(ds(n//m));total+=mult*f['dim']
            if f['analytic_rank']>0:forms.append({'oldform_multiplicity':mult,**f})
    assert total==g
    I=Fraction(n*phi(n),2)
    for p in ds(n):
        if p>1 and all(p%d for d in range(2,int(p**.5)+1)):I*=Fraction(p+1,p)
    out[n]={'genus':g,'weighted_newform_dimension':total,'positive_rank_forms':forms,
            'rank_zero':not forms,'Abramovich_lower_bound':str(Fraction(325,32768)*I)}
assert out[45]['rank_zero'] and out[55]['rank_zero'] and out[75]['rank_zero']
assert [f['label'] for f in out[57]['positive_rank_forms']]==['57.2.a.a']
assert [f['label'] for f in out[65]['positive_rank_forms']]==['65.2.a.a']
assert [f['label'] for f in out[63]['positive_rank_forms']]==['63.2.s.a']
clog=(H/'logs/check_characters.log').read_text()
assert 'indices=[ 47, 59 ] common_kernel={ 1, 5, 25, 38, 58, 62 }' in clog
assert 'OCTIC_REMAINING_RANK_CHARACTER_DONE' in clog
s=(H/'logs/winding69_octic.log').read_text()
rows=re.findall(r'WINDING N=69 order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+) time=',s)
assert len(rows)==4 and all(d==w for o,c,d,w,k in rows)
assert sum(int(d)*phi(int(o)) for o,c,d,w,k in rows)==133
assert 'bad_orders=[]' in s and 'WINDING_CHARACTER_PROBE_DONE' in s
out[69]={'genus':133,'weighted_winding_dimension':133,'rank_zero':True,
         'Abramovich_lower_bound':str(Fraction(325*2112,32768))}
assert Fraction(out[69]['Abramovich_lower_bound'])>16
assert Fraction(out[75]['Abramovich_lower_bound'])>16
out['all_cuspidal2']={
 '75':'excluded by rankzero and gonality>16; Parent implies2(D-C) principal, then base change gives h0(2D)>1 overQ',
 '69':'excluded by the same rankzero argument and gonality>16',
 '45':'not excluded by this bound: lower bound<16',
 '55':'not excluded by this bound: lower bound<16'}
files=['audit_remaining.py','check_characters.m','logs/check_characters.log','logs/winding69_octic.log','../common/winding_probe.m','../common/data/lmfdb_newforms_cache.json']
out['files']=[{'path':f,'sha256':hashlib.sha256((H/f).read_bytes()).hexdigest()} for f in files]
(H/'data/remaining_rank_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('OCTIC_REMAINING_RANK_AUDIT_DONE',','.join(map(str,[45,55,57,63,65,69,75])))
