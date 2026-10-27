#!/usr/bin/env python3
"""Independent group-law coverage check for all lower C40 supports."""
import json
from audit_fixed_octic_inputs import Field,TateCurve,literals,HERE

data = literals(HERE/'fixed_lower_inputs.py')
expected = {4:(4,32),5:(9,72),6:(10,80),7:(51,408)}
result = []
assert [row[0] for row in data['inputs']] == [4,5,6,7]
for degree,polynomial,representatives in data['inputs']:
    number,total = expected[degree]
    assert len(representatives) == number
    F = Field(polynomial)
    assert F.d == degree
    keys = set()
    sizes = []
    for raw in representatives:
        b,c = [F.encode(v) for v in raw]
        E = TateCurve(F,b,c)
        assert E.multiple(40) is None and E.multiple(20) is not None and E.multiple(8) is not None
        assert E.tate_pair(E.P) == (b,c)
        orbit = {F.key(*E.tate_pair(E.multiple(a))) for a in [1,3,7,9,11,13,17,19]}
        assert keys.isdisjoint(orbit)
        keys.update(orbit)
        sizes.append(len(orbit))
    assert len(keys) == total
    result.append({'degree':degree,'representatives':number,'closed_places':total,'orbit_sizes':sizes})
report = {'level':40,'degree':8,'supports':result,'exact_order40':True,
          'complete_lower_diamond_partitions':True}
(HERE/'fixed_lower_inputs_replay.json').write_text(json.dumps(report,indent=2)+'\n')
print('C40_LOWER_INPUT_INDEPENDENT_REPLAY_DONE representatives=74 closed_places=592 exact_order40=true')
