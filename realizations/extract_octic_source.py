#!/usr/bin/env python3
"""Derive the exact ten octic entries from the archived primary source."""
import hashlib
import json
from pathlib import Path
import re
import sympy as s

HERE=Path(__file__).resolve().parent
source=HERE/'data'/'LowDegreePlaces'
x,y,t=s.symbols('x y t')
entries=[]
for lineno,line in enumerate(source.read_text().splitlines(),1):
    m=re.match(r'N = (\d+), degv = 8,.*?\[(.*) = 0, (.*) = 0\]',line)
    if not m:
        continue
    n=int(m[1]);f=s.sympify(m[2].replace('^','**'));g=s.sympify(m[3].replace('^','**'))
    if s.degree(f,x)==8:
        roots=s.solve(g,y);assert len(roots)==1
        xp=t;yp=roots[0].subs(x,t);ff=f.subs(x,t)
    else:
        # A primitive absolute generator a=x+y for the quartic+quadratic case.
        G=s.groebner([f,g.subs(y,t-x)],x,t,order='lex')
        assert len(G.polys)==2
        roots=s.solve(G.polys[0].as_expr(),x);assert len(roots)==1
        xp=roots[0];yp=t-xp;ff=G.polys[1].as_expr()
    def coeffs(p):
        return [str(c) for c in reversed(s.Poly(p,t).all_coeffs())]
    entries.append({'index':len(entries)+1,'N':n,'source_line':lineno,
                    'source_text':line,'f':coeffs(ff),'x':coeffs(xp),'y':coeffs(yp)})
assert len(entries)==10
out=HERE/'data'
data={'source':'/home/fnajman/degree5_torsion/refs/vanhoeij/LowDegreePlaces',
      'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'examples':entries}
(out/'vanhoeij_octic_input.json').write_text(json.dumps(data,indent=2)+'\n')
rows=['<'+str(e['N'])+','+', '.join('['+','.join(e[k])+']' for k in ('f','x','y'))+'>' for e in entries]
(out/'vanhoeij_octic_input.m').write_text('OcticExamples:=[*\n'+',\n'.join(rows)+'\n*];\n')
print('OCTIC_PRIMARY_SOURCE_EXTRACT_DONE examples=10 complete=true')
