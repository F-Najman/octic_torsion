#!/usr/bin/env python3
"""Check the proved native presentation and its exact factor-base map."""
import ast
import json
from pathlib import Path
from sympy import Matrix,eye
from sympy.matrices.normalforms import hermite_normal_form

HERE = Path(__file__).resolve().parent/'data'
def literals(path):
    return {s.targets[0].id:ast.literal_eval(s.value) for s in ast.parse(path.read_text()).body}
raw = literals(HERE/'diagnostics/proven_relation_state.py')
geometry = literals(HERE/'factorbase_geometry.py')
assert raw['proven'] and raw['order'] == 31579088486400
relations = raw['abstract_relations']
images = raw['factorbase_abstract_classes']
degrees = list(reversed(raw['factorbase_degrees']))
assert geometry['relations'] == relations
assert geometry['factorbase_classes'] == images
assert geometry['factorbase_degrees'] == degrees
assert len(images) == 616 and len(relations) == 8
assert all(len(v) == 9 and v[-1] == d for v,d in zip(images,degrees))
# Prove that the saved factor-base images generate the full abstract Picard
# presentation, including its free degree-one direction.
full_lattice = hermite_normal_form(Matrix(relations+images).T)
assert full_lattice == eye(9)
R = Matrix([v[:-1] for v in relations])
H = hermite_normal_form(R.T)
assert abs(int(H.det())) == raw['order']
def zero(vector):
    assert len(vector) == 9 and vector[-1] == 0
    v = list(vector[:-1])
    for j in range(7,-1,-1):
        q = v[j]//int(H[j,j])
        for i in range(j+1):
            v[i] -= q*int(H[i,j])
    return not any(v)
assert len(raw['raw_relations']) == 766
for row in raw['raw_relations']:
    assert len(row) == 616
    # Native raw relation order is the unreversed factor-base order.
    actual = [sum(c*images[615-j][i] for j,c in enumerate(row)) for i in range(9)]
    assert zero(actual)
report = {'level':40,'degree':8,'native_proven_flag':True,
          'finite_order':raw['order'],'factorbase_primes':616,
          'raw_relations_checked':766,'factorbase_map_surjective':True,
          'degree_coordinates_exact':True,'actual_relation_lattice':True}
(HERE/'proved_picard_presentation_replay.json').write_text(json.dumps(report,indent=2)+'\n')
print('C40_PROVED_PICARD_PRESENTATION_REPLAY_DONE raw_relations=766 factorbase=616 surjective=true degree_exact=true')
