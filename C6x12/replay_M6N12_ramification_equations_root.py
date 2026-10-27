#!/usr/bin/env python3
"""Verify the24 degree8 ramification primes independently of native factorization."""
from pathlib import Path
import ast,json,re
from functools import lru_cache
from sympy import Add,Mul,Pow,Poly,QQ,Symbol,Rational,sqrt,sympify
from sympy.polys.polyclasses import ANP
P=Path(__file__).resolve().parent
u,v,z=map(Symbol,['u','v','z']);zr=(1+sqrt(-3))/2;K=QQ.algebraic_field(zr)
def parse(s):return sympify(s.replace('^','**'),locals={'u':u,'v':v,'z':z})
raw=(P/'data'/'M6N12_explicit_pencil_witness_arithmetic_root.txt').read_text()
# Sprint renders the two function strings without quotation marks.
raw,nquoted=re.subn(r'\[\*\s*(u\^8[^,\[\]]*),\s*([^,\[\]]*)\s*\*\]',
    lambda m:'[*'+json.dumps(m[1].strip())+','+json.dumps(m[2].strip())+'*]',raw)
assert nquoted==24
native=ast.literal_eval(raw.replace('[*','[').replace('*]',']'))
assert [r[0] for r in native]==list(range(127,133))
model=(P.parent/'mdmagma/models_X1_m_n/X1_6_12.txt').read_text()
curve=parse(re.search(r'Xz\s*:=\s*(.*?);',model,re.S)[1])
irreducible={};summaries=[]
for index,dim,deg,ds,ms,gens in native:
    assert (dim,deg,ds,ms)==(2,8,[8]*4,[1]*4) and len(gens)==4
    pencil=parse((P/'data'/f'M6N12_char0_pencil_{index}_expression_root.txt').read_text())
    distinct=set()
    for g1s,g2s in gens:
        g1,g2=parse(g1s),parse(g2s)
        assert not g1.has(v) and Poly(g1,u).degree()==8
        fg=Poly(g1.subs(z,zr),u,extension=zr).monic()
        fkey=str(fg.as_expr())
        if fkey not in irreducible:
            assert fg.is_irreducible;irreducible[fkey]=True
        mod=[K.from_sympy(c) for c in fg.all_coeffs()]
        def const(a):return ANP([K.convert(a)],mod,K)
        zero,one=const(0),const(1);uu=ANP([K.one,K.zero],mod,K)
        zz=const(zr)
        @lru_cache(None)
        def simple(e):
            if e==u:return uu
            if e==z:return zz
            if e.is_Number:return const(e)
            if e.is_Add:
                out=zero
                for q in e.args:out+=simple(q)
                return out
            if e.is_Mul:
                out=one
                for q in e.args:out*=simple(q)
                return out
            assert e.is_Pow and e.exp.is_Integer
            return simple(e.base)**int(e.exp)
        pg2=Poly(g2,v);assert pg2.degree()==1
        av,bv=simple(pg2.nth(1)),simple(pg2.nth(0));assert not av.is_zero
        vv=-bv/av
        pointkey=(fkey,str(vv.to_list()));assert pointkey not in distinct;distinct.add(pointkey)
        @lru_cache(None)
        def ad(e):
            if e==u:return uu,one,zero
            if e==v:return vv,zero,one
            if e==z:return zz,zero,zero
            if e.is_Number:return const(e),zero,zero
            if e.is_Add:
                out=[zero,zero,zero]
                for q in e.args:
                    a=ad(q);out=[x+y for x,y in zip(out,a)]
                return tuple(out)
            if e.is_Mul:
                a,b,c=one,zero,zero
                for q in e.args:
                    x,y,t=ad(q);a,b,c=a*x,b*x+a*y,c*x+a*t
                return a,b,c
            assert e.is_Pow and e.exp.is_Integer
            x,y,t=ad(e.base);n=int(e.exp)
            if n==0:return one,zero,zero
            a=x**n;factor=const(n)*x**(n-1)
            return a,factor*y,factor*t
        cv,cu,cw=ad(curve);assert cv.is_zero and not cw.is_zero
        fv,fu,fw=ad(pencil)
        assert (fu-fw*cu/cw).is_zero
        print('M6N12_INDEPENDENT_RAMIFICATION_PRIME',index,len(distinct),8,flush=True)
    assert len(distinct)==4
    summaries.append({'index':index,'independent_closed_degrees':[8]*4,
                      'all_on_curve':True,'all_smooth_in_u':True,'all_df_zero':True})
report={'degree':8,'group':[6,12],'pencils':summaries,'verified_primes':24,
 'distinct_irreducible_u_polynomials':len(irreducible),
 'argument':'Each independently irreducible degree8 u polynomial and linear v equation defines a closeddegree8 smooth point. All24 points lie on the exact characteristic-zero model and annul the implicit derivative of the corresponding explicit pencil. Four distinct degree8 points fill the entire ramificationdegree32 for each degree8 map of a genus9 curve, so no further ramification or quartic support is possible.',
 'boundary':'Independent polynomial irreducibility over Qzeta6, quotient-field arithmetic, curve membership and exact derivatives. Native Riemann--Roch dimension and pole-degree/class identities are checked separately.'}
(P/'data'/'replay_M6N12_ramification_equations_root.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('M6N12_RAMIFICATION_EQUATIONS_INDEPENDENT_REPLAY_DONE')
