"""Prepare standalone exact ideal-identity replays, with no group search.

Requires all eight completed arithmetic/witness streams. Rough points are
paired only with their actual successful RR witness; smooth points are
matched directly to the saved actual factorbase ideal.
"""
import ast
from pathlib import Path
import runpy

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'C48'/'data'
D=P/'clean_run_v2'
L=ROOT/'C48'/'logs'/'clean_run_v2'  # identity logs and their witnesses
S=ROOT/'C48'/'clean_run_v2'  # generated replay scripts
def literal(path):
    return {x.targets[0].id:ast.literal_eval(x.value) for x in ast.parse(path.read_text()).body}
def magma(v):
    if isinstance(v,int):
        return str(v)
    assert isinstance(v,list)
    if all(isinstance(x,int) for x in v):
        return '['+','.join(map(str,v))+']'
    return '[*'+','.join(magma(x) for x in v)+'*]'

geom=literal(P/'clean_run/factorbase_geometry.py')  # identical to the clean_run_v2 copy, merged
binding=literal(D/'function_field_bindings.py')
source=runpy.run_path(str(P/'source_partition_inputs_complete.py'))
orbits=runpy.run_path(str(P/'clean_run/high_target_aut_orbits.py'))
keys={r[2]:r[3] for r in orbits['cases']}
records=[]
commonA=None
for shard in range(8):
    ar=literal(D/f'high_shard{shard}_arithmetic.py')
    ws=literal(L/f'high_shard{shard}_identities.log.witnesses.py')['witnesses']
    assert ar['complete'] is True and ar['relations']==geom['relations']
    assert [r[0] for r in ar['tests']]==ar['expected_cluster_ids']
    stream=iter(ws)
    for cid,d,idx,e,v,passed in ar['tests']:
        if e>5:
            w=next(stream)
            assert e==d and w[1]==0 and w[2]==v
            if commonA is None:
                commonA=w[5]
            assert w[5]==commonA
        else:
            w=[]
        records.append([cid,d,e,v,w])
    assert next(stream,None) is None
assert sorted(r[0] for r in records)==list(range(1,3508))
assert commonA is not None

common=[('C48ReplayFieldEquation',geom['field_equation']),
        ('C48ReplayFlags',geom['finite_flags']),('C48ReplayDegrees',geom['factorbase_degrees']),
        ('C48ReplayIdealBases',geom['factorbase_ideal_bases']),('C48ReplayImages',geom['factorbase_classes']),
        ('C48ReplayRelations',geom['relations']),('C48ReplayA',commonA),
        ('C48ReplayCoordinates',binding['target_modular_coordinates']),
        ('C48ReplayBaseSignature',binding['base_cusp_signature']),
        ('C48ReplayFieldPolynomials',[source['field_polynomials'][d] for d in range(3,9)]),
        ('C48ReplayKeys',[keys[cid] for cid in range(1,3508)])]
with (S/'high_witness_replay_v3_common.m').open('w') as f:
    for name,value in common:
        f.write(f'{name} := {magma(value)};\n')

template='''/* Independent exact replay; no ClassGroup or RR search. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "common/pair_af_place_linear.m";
load "C48/clean_run_v2/high_witness_replay_v3_common.m";
load "CASEFILE";
procedure C48ReplayHighWitnesses()
 Y:=MDX11(2,24,GF(5));FC:=FunctionField(Curve(Y));F0,mp:=AlgorithmicFunctionField(FC);F:=RationalExtensionRepresentation(F0);
 B:=BaseField(F);k:=ConstantField(F);PP:=PolynomialRing(k);
 function RatData(a)
  return [[Integers()!x:x in Coefficients(Numerator(a))],[Integers()!x:x in Coefficients(Denominator(a))]];
 end function;
 assert [RatData(a):a in Coefficients(DefiningPolynomial(F))] eq [[[Integers()!u:u in z]:z in v]:v in C48ReplayFieldEquation];
 function GetFunction(data)
  return F![B!(PP!v[1])/B!(PP!v[2]):v in data];
 end function;
 coords:=[F!mp(z):z in Y`_coordinates];assert coords eq [GetFunction(v):v in C48ReplayCoordinates];
 OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[**];
 for j in [1..#C48ReplayDegrees] do
  OO:=C48ReplayFlags[j] eq 1 select OF else OI;
  I:=ideal<OO|[GetFunction(a):a in C48ReplayIdealBases[j]]>;
  assert Degree(Place(I)) eq C48ReplayDegrees[j];Append(~fb,I);
 end for;
 assert #fb eq 810;
 function PairDescriptor(values)
  K:=Parent(values[1]);kk:=PrimeField(K);ee:=LCM([Degree(MinimalPolynomial(v,kk)):v in values]);found:=false;
  for bits in [[1,0],[0,1]] cat [[1,a]:a in [1..4]] do
   gen:=&+[K!bits[i]*values[i]:i in [1,2]];pol:=MinimalPolynomial(gen,kk);
   if Degree(pol) eq ee then chosen:=bits;found:=true;break;end if;
  end for;
  assert found;M:=Matrix(kk,[Eltseq(gen^i,kk):i in [0..ee-1]]);assert Rank(M) eq ee;
  cs:=[Eltseq(Solution(M,Vector(kk,Eltseq(v,kk)))):v in values];
  assert &and[&+[K!cs[j][i+1]*gen^i:i in [0..ee-1]] eq values[j]:j in [1,2]];
  return [ee] cat chosen cat [Integers()!v:v in Coefficients(pol)] cat [Integers()!v:v in &cat cs];
 end function;
 function DivisorIdeals(sparse)
  IF:=ideal<OF|1>;II:=ideal<OI|1>;degree:=0;seen:={Integers()|};
  for row in sparse do
   j,e:=Explode(row);assert j in [1..810] and e gt 0 and j notin seen;Include(~seen,j);
   degree+:=e*C48ReplayDegrees[j];
   if C48ReplayFlags[j] eq 1 then IF*:=fb[j]^e;else II*:=fb[j]^e;end if;
  end for;
  return IF,II,degree;
 end function;
 AF,AI,degree:=DivisorIdeals(C48ReplayA);assert degree eq 1;
 if IsOne(AI) then APlace:=Place(AF);assert IsOne(AI);else assert IsOne(AF);APlace:=Place(AI);end if;
 assert Degree(APlace) eq 1;
 b,c:=Explode(coords);dd:=b*(1-b-c);E:=EllipticCurve([F!0,c,F!0,dd,F!0]);Q:=E![b,b];units:=[b,dd,c^2-4*dd];
 for i in [1..12] do Append(~units,(i*Q)[1]);end for;
 assert [Valuation(u,APlace):u in units] eq C48ReplayBaseSignature;
 avec:=[&+[row[2]*C48ReplayImages[row[1]][i]:row in C48ReplayA]:i in [1..9]];
 inv:=[C48ReplayRelations[i][i]:i in [1..8]];
 assert [[Integers()!u:u in row]:row in C48ReplayRelations] eq [[(i eq j select inv[i] else 0):j in [1..9]]:i in [1..8]];
 fields:=AssociativeArray();for d in [3..8] do fields[d]:=ext<GF(5)|Polynomial(GF(5),C48ReplayFieldPolynomials[d-2])>;end for;
 done:=[];rough:=0;smooth:=0;
 for item in C48ReplayCases do
  cid,d,e,v,w:=Explode(item);key:=C48ReplayKeys[cid];assert key[1] eq d;
  K:=fields[d];bc:=[&+[K!key[2+j*d+i]*K.1^i:i in [0..d-1]]:j in [0,1]];
  if #w gt 0 then
   assert e eq d and w[2] eq 0 and w[3] eq v and w[6] eq C48ReplayA;
   EF,EI,ed:=DivisorIdeals(w[4]);RF,RI,rd:=DivisorIdeals(w[5]);assert ed eq 13 and rd eq 13;
   aa:=GetFunction(w[7]);assert aa ne 0;
   // Reconstruct the actual point ideal FROM the saved E,R,function.
   // Primality, both orders and its intrinsic moduli pair are then checked
   // independently, avoiding the original common-zero inverse entirely.
   ID:=EF/(RF*(aa*OF));JD:=EI/(RI*(aa*OI));IP:=ID*AF^e;JP:=JD*AI^e;
   assert IsIntegral(IP) and IsIntegral(JP);
   factors:=[*<z[1],z[2]>:z in Factorization(IP)*];
   for z in Factorization(JP) do Append(~factors,<z[1],z[2]>);end for;
   assert #factors eq 1 and factors[1][2] eq 1;
   p:=Place(factors[1][1]);assert Degree(p) eq e;
   PF,PI:=Ideals(Divisor(p));assert PF eq IP and PI eq JP;
   actualbc:=[Evaluate(z,p):z in coords];assert PairDescriptor(actualbc) eq PairDescriptor(bc);
   DF,DI:=Ideals(Divisor(p)-e*Divisor(APlace));assert DF eq ID and DI eq JD;rough+:=1;
  else
   p:=C48PairAFPlaceLinear(coords,bc);assert Degree(p) eq e;
   assert e le 5;I:=Ideal(p);hits:=[j:j in [1..810]|C48ReplayDegrees[j] eq e and fb[j] eq I];assert #hits eq 1;j:=hits[1];
   vv:=[(d div e)*(C48ReplayImages[j][i]-e*avec[i]):i in [1..9]];
   assert vv[9] eq v[9] and &and[(vv[i]-v[i]) mod inv[i] eq 0:i in [1..8]];smooth+:=1;
  end if;
  Append(~done,cid);
  printf "C48_WITNESS_REPLAY_CASE cluster=%o both_orders_exact=true time=%o\\n",cid,Cputime(tt);
 end for;
 assert done eq C48ReplayExpectedIDs;
 printf "C48_WITNESS_REPLAY_DONE shard=SHARD tested=%o rough=%o smooth=%o complete=true time=%o\\n",#done,rough,smooth,Cputime(tt);
end procedure;
C48ReplayHighWitnesses();quit;
'''
for shard in range(8):
    selected=sorted((r for r in records if (r[0]-1)%8==shard),key=lambda r:r[0])
    data=D/f'high_witness_replay_v3_shard{shard}_data.m'
    data.write_text('C48ReplayCases := '+magma(selected)+';\nC48ReplayExpectedIDs := '+magma([r[0] for r in selected])+';\n')
    body=template.replace('CASEFILE',str(data.relative_to(ROOT))).replace('SHARD',str(shard))
    (S/f'high_witness_replay_v3_shard{shard}.m').write_text(body)
print('C48_HIGH_WITNESS_REPLAY_V3_PREPARED cases=3507 shards=8')
