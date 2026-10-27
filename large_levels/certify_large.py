#!/usr/bin/env python3
"""Exact octic inequalities, complete rank provenance, and CM minima.

The global-q3 argument is independent of dyadic good-reduction points.
The extra local case95 requires the separately completed degree8 F2 log.
"""
import ast
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
INPUTS = HERE / 'inputs'  # the copied inputs now live in common/data, logs/ and data/
COMMON = HERE.parent / 'common' / 'data'
CACHE = json.loads((COMMON / 'lmfdb_newforms_cache.json').read_text())
CM_TEXT = (COMMON / 'dcm_list_all_min_orders_X1_10k.m').read_text()
CM = {int(n): int(d) for n, d, _ in re.findall(
    r'\[\*\s*(\d+),\s*(\d+),\s*(\d+),', CM_TEXT)}
assert sorted(CM) == list(range(1, 10001))
GLOBAL = [(119,3),(121,2),(125,3),(133,64),(143,67),(169,3),
          (187,122),(209,34),(221,3),(247,144),(289,3),
          (323,101),(361,3)]
WINDING = [115,161,253,299,391,437,529]
TARGETS = sorted(n for n, _ in GLOBAL) + WINDING + [95]


def divisors(n):
    return [i for i in range(1,n+1) if n % i == 0]


def phi(n):
    return sum(math.gcd(i,n) == 1 for i in range(1,n+1))


def index(n):
    ans = Fraction(n*n, 2)
    for p in divisors(n):
        if p > 1 and len(divisors(p)) == 2:
            ans *= Fraction(p*p-1,p*p)
    assert ans.denominator == 1
    return int(ans)


def diamond_order(n,a):
    assert math.gcd(n,a) == 1
    power, order = a%n, 1
    while power not in (1,n-1):
        power, order = power*a%n, order+1
    return order


def rank_input(n):
    total, positive = 0, []
    for m in divisors(n):
        assert str(m) in CACHE
        mult = len(divisors(n//m))
        for f in CACHE[str(m)]:
            assert f['analytic_rank_proved'] is True
            assert f['analytic_rank'] is not None
            total += mult*f['dim']
            if f['analytic_rank'] > 0:
                positive.append({'level':m,'multiplicity':mult,**f})
    return total, positive


def winding_input(n):
    name = 'winding529_krylov_modp_v2.log' if n == 529 else f'winding{n}.log'
    source = (HERE/'logs'/name).read_text()
    assert 'WINDING_CHARACTER_PROBE_DONE' in source
    assert not re.search(r'Runtime error|Assertion failed|SystemError|User error',source)
    rows = re.findall(r'WINDING N=(\d+) order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+) time=([0-9.]+)',source)
    assert rows
    chars = [dict(zip(['N','order','conductor','dimension','winding_dimension',
                      'coefficient_degree'],map(int,r[:6]))) for r in rows]
    assert all(r['N'] == n and (r['order'] == 1 or
               r['dimension'] == r['winding_dimension']) for r in chars)
    final = re.findall(r'WINDING_FINAL N=(\d+) genus=(\d+) total_dimension=(\d+) bad_orders=(\[.*?\]) admissible=(\[.*?\])',source)
    assert len(final) == 1
    nn,g,total = map(int,final[0][:3])
    bad = ast.literal_eval(final[0][3])
    assert nn == n and g == total and bad in ([],[1])
    assert 3 in ast.literal_eval(final[0][4])
    assert sum(phi(c['order'])*c['dimension'] for c in chars) == g
    if n == 529:
        replay = json.loads((HERE/'data'/'winding529_krylov_replay.json').read_text())
        assert replay['verified'] and replay['Krylov_rank'] == replay['dimension'] == 45
        assert replay['splitting_prime'] == 2531
        assert 'WINDING_GOOD_PRIME p=2531 all_denominators_invertible=true' in source
        assert 'WINDING_KRYLOV N=529 order=253 span_rank=45 dimension=45 verified=true' in source
    return g, chars, name


rows, character_cases = [], []
for n,a,method in [(n,a,'global3') for n,a in GLOBAL] + [
        (n,3,'global3') for n in WINDING] + [(95,11,'local2')]:
    idx = index(n)
    genus = 1 + Fraction(idx,12) - Fraction(
        sum(phi(d)*phi(n//d) for d in divisors(n)),4)
    assert genus.denominator == 1
    gon = Fraction(325*idx,32768)
    order = diamond_order(n,a)
    assert 8 % order and CM[n] > 8
    row = {'N':n,'a':a,'method':method,'index':idx,'genus':int(genus),
           'gonality_lower_bound':str(gon),'diamond_order':order,
           'minimum_CM_degree':CM[n]}
    if n in WINDING:
        g,chars,name = winding_input(n)
        assert g == genus
        row.update(rank_source='exact winding submodules',winding_log='inputs/'+name,
                   character_orbits=chars)
    else:
        total,positive = rank_input(n)
        assert total == genus
        nontriv = sorted({(f['level'],f['char_orbit_label']) for f in positive
                          if f['char_order'] > 1})
        assert all(m == n for m,_ in nontriv)
        character_cases.append((n,a,[lab for _,lab in nontriv]))
        row.update(rank_source='proved archived newforms',
                   newform_dimension_sum=total,positive_rank_forms=positive,
                   nontrivial_positive_characters=nontriv)
    if method == 'global3':
        k = 7 if a in (3,pow(3,-1,n)) else 8
        assert n%3 and 8*k < gon
        row.update(q=3,kq=k,principal_function_degree_bound=8*k)
    else:
        ff = (HERE/'logs'/'finite_fields_2_octic.log').read_text()
        assert ff.rstrip().endswith('DONE')
        assert re.search(rf'TORSION n={n} p=2 : k with Z/{n} over F_p\^k = \[\s*\]',ff)
        assert gon > 32
        row.update(p=2,principal_function_degree_bound=32,
                   good_reduction_degrees_through8=[])
    rows.append(row)

script = ['SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);',
          'AttachSpec("mdmagma/Magma/magma.spec");']
for n,a,labs in character_cases:
    script.append(f'N:={n};a:={a};')
    for lab in labs:
        script += [f'idx:=[j:j in [1..N-1]|GCD(j,N) eq 1 and ConreyCharacterOrbitLabel(N,j) eq Sprintf("%o.%o",N,"{lab}")];',
                   'assert #idx gt 0 and &and[ConreyCharacterValue(N,j,a) eq 1:j in idx];',
                   f'printf "OCTIC_CHARACTER_KERNEL N=%o a=%o orbit=%o indices=%o verified=true\\n",N,a,"{lab}",idx;']
    script.append(f'printf "OCTIC_FINITE_IMAGE_INPUT N=%o a=%o nontrivial_positive_characters=%o verified=true\\n",N,a,{len(labs)};')
script += ['print "OCTIC_LARGE_CHARACTER_CERTIFICATE_DONE";','quit;']
(HERE/'check_large_characters.m').write_text('\n'.join(script)+'\n')
out = {'degree':8,'excluded_cyclic_orders':sorted(r['N'] for r in rows),
       'rows':sorted(rows,key=lambda r:r['N']),
       'unresolved_original_large_frontier':[77,81,85,91],
       'status':'all octic numerical, complete rank and CM inputs verified; exact character log required',
       'rank_input_sha256':hashlib.sha256((COMMON/'lmfdb_newforms_cache.json').read_bytes()).hexdigest(),
       'CM_input_sha256':hashlib.sha256((COMMON/'dcm_list_all_min_orders_X1_10k.m').read_bytes()).hexdigest()}
check_log = HERE/'logs'/'check_large_characters.log'
if check_log.exists():
    text = check_log.read_text()
    assert text.rstrip().endswith('OCTIC_LARGE_CHARACTER_CERTIFICATE_DONE')
    assert not re.search(r'Runtime error|Assertion failed|SystemError|User error',text)
    for n,a,labs in character_cases:
        assert f'OCTIC_FINITE_IMAGE_INPUT N={n} a={a} nontrivial_positive_characters={len(labs)} verified=true' in text
        for lab in labs:
            assert f'OCTIC_CHARACTER_KERNEL N={n} a={a} orbit={lab} ' in text
    out['status'] = 'complete exact degree8 numerical, rank, character, CM and finite-field certificate'
    out['exact_character_log_sha256'] = hashlib.sha256(check_log.read_bytes()).hexdigest()
(HERE/'data'/'large_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('OCTIC_LARGE_NUMERICAL_DONE',len(rows),out['excluded_cyclic_orders'])
for r in out['rows']:
    print(r['N'],r['a'],r['method'],r['diamond_order'],r['gonality_lower_bound'],r['minimum_CM_degree'])
