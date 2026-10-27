#!/usr/bin/env python3
"""Exact minimal-obstruction ledger for Phi(8).

The classification is complete exactly when every minimal obstruction
has an audited exclusion record. Completion is recomputed from evidence.
Pairs (m,n), with m|n, denote C_m x C_n. The Weil pairing forces phi(m)|8, and Khawaja gives S(8)=Primes(23).
"""
from pathlib import Path
import hashlib
import json
from sympy import divisors, primefactors, totient

HERE = Path(__file__).resolve().parent/'data'
ROOT = Path(__file__).resolve().parents[1]
PRIMES = {2, 3, 5, 7, 11, 13, 17, 19, 23}
# phi(m)|8 forces p-1|8 for p|m, so p is 2,3 or5. Odd prime
# exponents are at most1, and the 2-exponent is at most4.
FIRST_FACTORS = {2**a*3**b*5**c for a in range(5) for b in range(2)
                 for c in range(2) if 8 % int(totient(2**a*3**b*5**c)) == 0}
assert FIRST_FACTORS == {1,2,3,4,5,6,8,10,12,15,16,20,24,30}
INFINITE = ({(1,n) for n in range(1,29)} | {(1,n) for n in (30,32,36)}
            | {(2,2*n) for n in range(1,13)}
            | {(3,3*n) for n in range(1,5)}
            | {(4,4*n) for n in range(1,4)}
            | {(5,5),(6,6)})
assert len(INFINITE) == 52



def subgroup(a, b):
    return b[0] % a[0] == 0 and b[1] % a[1] == 0


def subgroups(g):
    return {(int(a), int(b)) for a in divisors(g[0])
            for b in divisors(g[1]) if b % a == 0}


def frontier(base):
    assert all(subgroups(g) <= base for g in base)
    cyclic = {n for m, n in base if m == 1}
    # For a minimal cyclic obstruction n, either n is prime or every n/p
    # belongs to the cyclic base. Thus n <= max(base cyclic)*max(PRIMES).
    bound = max(cyclic)*max(PRIMES)
    answer = []
    for m in sorted(FIRST_FACTORS):
        ns = range(1, bound+1) if m == 1 else sorted(cyclic)
        for n in ns:
            g = (m, n)
            if n % m or g in base or not set(primefactors(n)) <= PRIMES:
                continue
            if subgroups(g)-{g} <= base:
                answer.append(g)
    return sorted(answer)


def label(g):
    return f'C{g[1]}' if g[0] == 1 else f'C{g[0]} x C{g[1]}'


def read_certificates():
    # The groups of the accepted proof ledger, with the folder of this repository that holds
    # each proof (common/data/ledger.json; no folder for a group excluded by a published
    # result); see README. The SHA-256 of the evidence files are checked by the tests.
    path = HERE/'ledger.json'
    data = json.loads(path.read_text())
    for kind in ('realizations', 'exclusions'):
        for record in data[kind]:
            assert record['audited'] is True
            assert record['degree'] == 8
            assert record['folder'] is None or (ROOT/record['folder']).is_dir(), record['folder']
    return data


def main():
    records = read_certificates()
    known = INFINITE | {tuple(r['group']) for r in records['realizations']}
    excluded = {tuple(r['group']) for r in records['exclusions']}
    assert not any(subgroup(e, g) for e in excluded for g in known)
    minimal = frontier(known)
    eliminated = [g for g in minimal if any(subgroup(e, g) for e in excluded)]
    remaining = [g for g in minimal if g not in eliminated]
    data = {
        'degree': 8,
        'status': 'COMPLETE' if not remaining else 'INCOMPLETE RESEARCH',
        'classification_complete': not remaining,
        'prime_divisors': sorted(PRIMES),
        'first_invariant_factors': sorted(FIRST_FACTORS),
        'phi_infinite': sorted(INFINITE),
        'known_lower_bound': sorted(known),
        'classification': sorted(known) if not remaining else None,
        'minimal_outside_infinite': frontier(INFINITE),
        'minimal_outside_known': minimal,
        'certified_exclusions': sorted(excluded),
        'eliminated_minimal_obstructions': eliminated,
        'remaining_minimal_obstructions': remaining,
        'sources': {
            'infinite_families': 'Najman--Varivoda, arXiv:2602.03513v4, Theorem 1.1(2)',
            'prime_bound': 'Khawaja, arXiv:2304.14284v3, Theorem 1.2',
            'warning': 'No exclusion in another degree is automatically transferred to degree 8.'
        }
    }
    (HERE/'candidates.json').write_text(json.dumps(data, indent=2)+'\n')
    text = [f'Phi(8): {data["status"]}',
            f'Published infinite-family groups: {len(INFINITE)}',
            f'Certified occurring groups: {len(known)}',
            f'Certified exclusion records: {len(excluded)}',
            f'Remaining minimal obstruction groups: {len(remaining)}', '',
            'Known occurring groups:', ', '.join(map(label, sorted(known))), '',
            'Remaining minimal obstruction groups:', ', '.join(map(label, remaining)), '',
            'A minimal obstruction is a subgroup which must be excluded or realized.',
            'It is not an assertion that this subgroup actually occurs.',
            'Realizing a new group can expose additional larger minimal obstructions.',
            'CANDIDATES_DEGREE8_DONE']
    (HERE/'candidates.txt').write_text('\n'.join(text)+'\n')
    print('\n'.join(text))


if __name__ == '__main__':
    main()
