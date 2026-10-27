"""Check both proved Picard presentations without identifying their bases."""
import ast
import hashlib
import json
from pathlib import Path
from sympy import Matrix,eye
from sympy.matrices.normalforms import hermite_normal_form

HERE=Path(__file__).resolve().parent/'data'
def literals(path):
    return {s.targets[0].id:ast.literal_eval(s.value) for s in ast.parse(path.read_text()).body}

for dirname in ('clean_run','clean_run_v2'):
    P=HERE/dirname
    raw=literals(P/'proven_relation_state.py')
    geometry=literals(HERE/'clean_run'/'factorbase_geometry.py')  # the clean_run_v2 copy was identical, merged
    assert raw['proven'] is True and raw['order']==4246732800
    relations=raw['abstract_relations'];images=raw['factorbase_images']
    degrees=raw['generator_factorbase_degrees']
    assert degrees==list(reversed(raw['raw_factorbase_degrees']))
    assert geometry['relations']==relations and geometry['factorbase_classes']==images
    assert geometry['factorbase_degrees']==degrees
    n=len(images)
    assert n==810 and len(relations)==8
    assert all(len(v)==9 and v[-1]==d for v,d in zip(images,degrees))
    assert hermite_normal_form(Matrix(relations+images).T)==eye(9)
    H=hermite_normal_form(Matrix([v[:-1] for v in relations]).T)
    assert abs(int(H.det()))==raw['order']
    def zero(vector):
        assert len(vector)==9 and vector[-1]==0
        v=list(vector[:-1])
        for j in range(7,-1,-1):
            q=v[j]//int(H[j,j])
            for i in range(j+1):
                v[i]-=q*int(H[i,j])
        return not any(v)
    for row in raw['raw_relations']:
        assert len(row)==n
        actual=[sum(c*images[n-1-j][i] for j,c in enumerate(row)) for i in range(9)]
        assert zero(actual)
    report={'complete':True,'native_proven_flag':True,'finite_order':raw['order'],
            'factorbase_primes':n,'raw_relations_checked':len(raw['raw_relations']),
            'factorbase_map_surjective':True,'degree_coordinates_exact':True,
            'actual_relation_lattice_checked':True,
            'state_sha256':hashlib.sha256((P/'proven_relation_state.py').read_bytes()).hexdigest(),
            'geometry_sha256':hashlib.sha256((HERE/'clean_run'/'factorbase_geometry.py').read_bytes()).hexdigest()}
    (P/'proved_presentation_replay.json').write_text(json.dumps(report,indent=2)+'\n')
    print('C48_PROVED_PRESENTATION_PASS',dirname,n,len(raw['raw_relations']),flush=True)
