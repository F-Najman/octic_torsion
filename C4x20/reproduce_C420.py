#!/usr/bin/env python3
"""Optional full native reproduction in a NEW directory (about 800CPU sec)."""
import argparse,ast,subprocess
from pathlib import Path
P=Path(__file__).resolve().parent;ROOT=P.parent
ap=argparse.ArgumentParser();ap.add_argument('--directory',type=Path,required=True);args=ap.parse_args()
D=args.directory.resolve();assert not D.exists(),'Use a new directory; frozen outputs must not be overwritten.'
assert D.is_relative_to(ROOT),'Keep reproduction inside the repository.'
D.mkdir(parents=True)
def run(name,assignment,source,marker,timeout):
    wrapper=D/(name+'.m');console=D/(name+'_console.log')
    wrapper.write_text(assignment+'\nload '+repr(str(source.relative_to(ROOT))).replace("'",'"')+';\n')
    with console.open('w') as out:
        result=subprocess.run(['magma','-b',str(wrapper)],cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,timeout=timeout)
    text=console.read_text();assert result.returncode==0 and marker in text,text[-3000:]
    assert not any(s in text for s in ['Runtime error','User error','Assertion failed','syntax error']),text[-3000:]
    print('C420_NATIVE_STAGE_PASS',name,flush=True)
run('positive','C420PositiveStem:="'+str(D/'positive')+'";',P/'C420_positive_reproduce.m','C420_POSITIVE_COMPLETE nc=',180)
run('all_supports','C420FullStem:="'+str(D/'all_supports')+'";',P/'C420_full_reproduce.m','C420_FRESH_COMPLETE singles=18 pairs=192 survivors=4',1800)
run('compact','C420LogPath:="'+str(D/'compact_points.log')+'";',P/'C420_compact_completion_final.m','C420_COMPACT_COMPLETION_DONE globalpoints=8 cusps=4 CM1728=4 complete=true',120)
def rd(path):return ast.literal_eval(path.read_text().replace('[*','[').replace('*]',']'))
assert rd(D/'positive_arithmetic.txt')==rd(P/'data'/'C420_positive_fresh_arithmetic.txt')
tests,remaining=rd(D/'all_supports_arithmetic.txt');assert remaining==[[1,4,4,0],[1,4,12,0],[2,2,1,82],[2,2,2,84]]
cases={tuple(t[:4]) for t in tests};assert len(cases)==210
assert rd(D/'all_supports_stabilizers.txt')==rd(ROOT/'C4x20/data/M4N20_stabilizers_root_arithmetic.txt')
print('C420_FULL_NATIVE_REPRODUCTION_PASS positive_places=1304 support_cases=210 compact_points=8 complete=true')
