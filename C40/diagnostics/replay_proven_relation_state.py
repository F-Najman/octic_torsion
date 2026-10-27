#!/usr/bin/env python3
"""Check the saved native proven-kernel data and abstract class map.

This is an integer consistency replay of the recovered native data, not a
replacement for Magma's proved factor-base generation/relation search.
"""
import ast
from hashlib import sha256
import json
from math import prod
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]/'data'/'diagnostics'
path=HERE/'proven_relation_state.py'
data={st.targets[0].id:ast.literal_eval(st.value) for st in ast.parse(path.read_text()).body}
assert data['proven'] is True and data['order']==31579088486400
moduli=[2,2,2,2,2,16,68880,895440]
assert prod(moduli)==31579088486400
assert data['abstract_relations']==[[m if j==i else 0 for j in range(9)] for i,m in enumerate(moduli)]
degrees=data['factorbase_degrees'];images=data['factorbase_abstract_classes'];relations=data['raw_relations']
assert len(images)==len(degrees)==616 and len(relations)==766
assert all(len(v)==9 for v in images)
assert all(v[-1]==d for v,d in zip(images,reversed(degrees)))
for r in relations:
    assert len(r)==616 and sum(a*d for a,d in zip(r,degrees))==0
    image=[sum(a*v[j] for a,v in zip(reversed(r),images)) for j in range(9)]
    assert image[-1]==0 and all(x%m==0 for x,m in zip(image,moduli))
report={'level':40,'degree':8,'native_relation_proof_preserved':True,
        'raw_relation_vectors':766,'factorbase_places':616,
        'all_relations_map_to_zero':True,'degree_coordinate_matches_reversed_factorbase':True,
        'abstract_class_number':prod(moduli),'sha256':sha256(path.read_bytes()).hexdigest(),
        'scope':'integer consistency; native proven generation and principal relations remain Magma evidence'}
(HERE/'proven_relation_state_replay.json').write_text(json.dumps(report,indent=2)+'\n')
print('C40_RECOVERED_RELATION_REPLAY_DONE factorbase=616 raw_relations=766 all_degree0=true all_kernel_images0=true exact_order=31579088486400')
