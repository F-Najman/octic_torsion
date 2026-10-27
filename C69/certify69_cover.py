#!/usr/bin/env python3
"""Audit the completed octic69 cover computation and bind its evidence."""
from pathlib import Path
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
BASE = HERE.parent


def main():
    log = (HERE / 'logs' / 'N69_cover_good_D8.log').read_text()
    assert not re.search(r'Runtime error|Assertion failed|System Error|User error', log, re.I)
    assert re.search(r'C69_D8_GOOD_DONE nc7_cases=3 nc7_remaining=0 nc8_cases=2 nc8_tests=0 nc8_remaining=0 complete=true time=[0-9.]+\s*$', log)
    assert 'COVER_INDEPENDENT_FIELD_VERIFIED 69 13' in log
    assert set(map(int, re.findall(r'COVER_GENERATOR_FROBENIUS_VERIFIED (\d+)', log))) == {1, 2, 3}
    assert 'C69_COVER_CUSP_CHARACTERS 2' in log
    base7 = re.findall(r'C69_BASE_A5_TEST (\d+).*?principal (true|false)', log, re.S)
    assert len(base7) == 3 and sorted(x[1] for x in base7) == ['false', 'false', 'true']
    lifted = re.findall(r'C69_COVER_CASE_TEST (\d+) (\d+) principal (true|false)', log)
    assert len(lifted) == 8 and {x[1] for x in lifted} == set(map(str, range(1, 9)))
    assert all(x[2] == 'false' for x in lifted)
    assert 'C69_D8_BASE_IMAGES curves=16 images=2 marked_points=352' in log
    base8 = re.findall(r'C69_D8_BASE_TEST index=(\d+) base_degree=(\d+) principal=(true|false)', log)
    assert base8 == [('1', '8', 'false'), ('2', '8', 'false')]
    source = (HERE / 'N69_cover_good_D8.m').read_text()
    assert ':i in [0..Degree(Parent(r[1]))-1]' in source
    assert 'RelativeDecomposition' in source and 'Valuation(F!u,Q) eq 1' in source
    rank_manifest = HERE / 'data' / 'remaining_rank_certificate.json'
    rank_data = json.loads(rank_manifest.read_text())
    assert rank_data['69']['rank_zero'] and rank_data['69']['weighted_winding_dimension'] == 133
    winding = (HERE / 'logs' / 'winding69_octic.log').read_text()
    rows = re.findall(r'WINDING N=69 order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+)', winding)
    assert rows and sum(int(d) * int(c) for _, _, d, _, c in rows) == 133
    assert all(d == w for _, _, d, w, _ in rows)
    assert winding.count('WINDING_CHARACTER_PROBE_DONE') == 1
    assert not re.search(r'Runtime error|Assertion failed|System Error|User error', winding, re.I)
    cusp = json.loads((HERE / 'data' / 'N69_cusp_source_audit.json').read_text())
    assert cusp['all_F2_rational_source_cusps_Q_rational']
    assert cusp['diamond5_trivial_on_quotient']
    files = {
        BASE / name for name in ('C69/N69_cover_good_D8.m', 'C69/logs/N69_cover_good_D8.log',
        'C69/N69_cusp_source_audit.py', 'C69/data/N69_cusp_source_audit.json',
        'C69/logs/N69_cusp_source_audit.log', 'C69/certify69_cover.py', 'common/independent_counts.m',
        'common/logs/independent_counts.log', 'C69/audit_remaining.py', 'C69/check_characters.m',
        'C69/logs/check_characters.log', 'C69/logs/winding69_octic.log', 'common/winding_probe.m',
        'common/data/lmfdb_newforms_cache.json', 'C69/audit69_finite.m', 'C69/logs/audit69_finite.log',
        'C69/audit69_cusps.py', 'C69/data/audit69_cusps.json', 'C69/logs/audit69_cusps.log')
    }
    files.add(rank_manifest)
    entries = []
    for p in sorted({p.resolve() for p in files}):
        entries.append({'path': str(p.relative_to(BASE)),
                        'sha256': hashlib.sha256(p.read_bytes()).hexdigest()})
    result = {'degree': 8, 'excluded_group': [1, 69], 'status': 'complete',
              'source_genus': 133, 'quotient_genus': 13, 'base_genus': 7,
              'closed_source_counts_7_8': [66, 44], 'base_images_7_8': [3, 2],
              'good7_lift_tests': 8, 'good8_base_tests': 2,
              'all_bad_branch': 'rank zero and rational gonality greater than16',
              'source_cusp_complement': 'Every F2-rational source cusp is Q-rational.',
              'path_base': 'work/degree8', 'files': entries}
    (HERE / 'data' / 'N69_cover_manifest.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'files'}, indent=2))
    print('OCTIC69_COMPLETE_CERTIFICATE_DONE')


if __name__ == '__main__':
    main()
