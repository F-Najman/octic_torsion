#!/usr/bin/env python3
"""Prepare a standalone geometric replay; no class-group search is needed.

Default: require all four complete shards and export all 96 witnesses.
--indices 1,2 permits an explicitly partial diagnostic while shards run.
"""
import argparse
import ast
from pathlib import Path

HERE = Path(__file__).resolve().parent/'data'

def literals(path,partial=False):
    text = path.read_text()
    if partial and not text.rstrip().endswith('True'):
        text += '\n]\n'
    data = {}
    for statement in ast.parse(text).body:
        assert isinstance(statement,ast.Assign) and len(statement.targets) == 1
        data[statement.targets[0].id] = ast.literal_eval(statement.value)
    return data

def magma(value):
    if isinstance(value,int):
        return str(value)
    assert isinstance(value,list)
    if all(isinstance(v,int) for v in value):
        return '['+','.join(str(v) for v in value)+']'
    return '[*'+','.join(magma(v) for v in value)+'*]'

parser = argparse.ArgumentParser()
parser.add_argument('--indices',help='Explicit diagnostic subset, e.g. 1,2')
parser.add_argument('--stem',help='Separate output name; preserves previous diagnostics')
args = parser.parse_args()
wanted = set(map(int,args.indices.split(','))) if args.indices else set(range(1,97))
partial = bool(args.indices)
geometry = literals(HERE/'factorbase_geometry.py')
master = literals(HERE/'fixed_octic_inputs.py')
bindings = literals(HERE/'function_field_bindings.py')
records = {}
for shard in sorted({(i-1)%4 for i in wanted}):
    arithmetic = literals(HERE/f'fork_shard{shard}_arithmetic.py',partial)
    witness_path = HERE.parent/'logs'/f'fork_shard{shard}_identities.log.witnesses.py'
    # A partial witness stream has no trailing complete assignment.
    text = witness_path.read_text()
    if not text.rstrip().endswith(']'):
        text += '\n]\n'
    witnesses = ast.literal_eval(ast.parse(text).body[0].value)
    if not partial:
        assert arithmetic['complete'] and len(witnesses) == len(arithmetic['tests']) == 24
    for test,witness in zip(arithmetic['tests'],witnesses):
        index,vector,member = test
        assert vector == witness[2] and index not in records
        records[index] = witness
assert wanted <= records.keys()
stem = args.stem or ('witness_replay_partial' if partial else 'witness_replay_all')
assert stem and all(c.isalnum() or c == '_' for c in stem)
data_path = HERE/f'{stem}_data.m'
with data_path.open('w') as output:
    for name,value in [
        ('C40ReplayFieldEquation',geometry['field_equation']),
        ('C40ReplayFiniteFlags',geometry['finite_flags']),
        ('C40ReplayDegrees',geometry['factorbase_degrees']),
        ('C40ReplayIdealBases',geometry['factorbase_ideal_bases']),
        ('C40ReplayFieldPolynomial',master['field_polynomial']),
        ('C40ReplayTateBC',master['representatives_bc']),
        ('C40ReplayTateFunctions',bindings['universal_tate_coordinates']),
        ('C40ReplaySutherlandFunctions',bindings['sutherland_coordinates']),
        ('C40ReplayBaseCuspSignature',bindings['base_cusp_signature']),
        ('C40ReplayCases',[[i,records[i]] for i in sorted(wanted)])]:
        output.write(f'{name} := {magma(value)};\n')
source = f'''/* Standalone direct-principal-identity replay.  Explicit cases {sorted(wanted)}. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");
Attach("common/fast_hecke_11.m");
load "{data_path.relative_to(HERE.parents[1])}";
X:=MDX1(40,GF(3));C:=Curve(X);
F0,coordinateMap:=AlgorithmicFunctionField(FunctionField(C));F:=RationalExtensionRepresentation(F0);
B:=BaseField(F);kk:=ConstantField(F);PP:=PolynomialRing(kk);
function RatData(a)
 return [[Integers()!x:x in Coefficients(Numerator(a))],[Integers()!x:x in Coefficients(Denominator(a))]];
end function;
assert [RatData(a):a in Coefficients(DefiningPolynomial(F))] eq [[[Integers()!u:u in z]:z in v]:v in C40ReplayFieldEquation];
function GetFunction(data)
 return F![B!(PP!v[1])/B!(PP!v[2]):v in data];
end function;
assert [F!coordinateMap(z):z in X`_coordinates] eq [GetFunction(v):v in C40ReplayTateFunctions];
assert [F!coordinateMap(z):z in X`_coordinates_xy] eq [GetFunction(v):v in C40ReplaySutherlandFunctions];
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
fb:=[**];
for j in [1..#C40ReplayIdealBases] do
 OO:=C40ReplayFiniteFlags[j] eq 1 select OF else OI;
 I:=ideal<OO|[GetFunction(a):a in C40ReplayIdealBases[j]]>;
 assert Degree(Place(I)) eq C40ReplayDegrees[j];Append(~fb,I);
end for;
assert #fb eq 616;
function DivisorIdeals(sparse)
 IF:=ideal<OF|1>;II:=ideal<OI|1>;degree:=0;seen:={{Integers()|}};
 for row in sparse do
  j,e:=Explode(row);assert j in [1..616] and e gt 0 and j notin seen;Include(~seen,j);
  degree+:=e*C40ReplayDegrees[j];
  if C40ReplayFiniteFlags[j] eq 1 then IF*:=fb[j]^e;else II*:=fb[j]^e;end if;
 end for;
 return IF,II,degree;
end function;
vA:=C40ReplayCases[1][2][6];AF,AI,degree:=DivisorIdeals(vA);assert degree eq 1;
basecusps:=[];
for cusp in Cusps(X) do
 if Degree(cusp) ne 1 then continue;end if;
 I,J:=Ideals(DivisorGroup(F)!FunctionFieldDivisor(Divisor(cusp)));
 if I eq AF and J eq AI then Append(~basecusps,cusp);end if;
end for;
assert #basecusps eq 1;c0:=basecusps[1];assert CuspSignature(X,c0) eq C40ReplayBaseCuspSignature;
K:=ext<GF(3)|PolynomialRing(GF(3))!C40ReplayFieldPolynomial>;
for item in C40ReplayCases do
 index:=item[1];w:=item[2];assert w[2] eq 0 and w[6] eq vA;
 EF,EI,ed:=DivisorIdeals(w[4]);RF,RI,rd:=DivisorIdeals(w[5]);assert ed eq 25 and rd eq 25;
 b:=GetFunction(w[7]);assert b ne 0;
 bc:=[Evaluate(PolynomialRing(K)!v,K.1):v in C40ReplayTateBC[index]];
 P:=MDPlaceCommonZeros(X`_coordinates,bc);assert Degree(P) eq 8;
 DD:=DivisorGroup(F)!FunctionFieldDivisor(Divisor(P)-8*Divisor(c0));
 ID,JD:=Ideals(DD);
 assert EF/(RF*(b*OF)) eq ID and EI/(RI*(b*OI)) eq JD;
 printf "C40_INDEPENDENT_PRINCIPAL_WITNESS index=%o both_orders_exact=true time=%o\\n",index,Cputime(tt);
end for;
printf "C40_INDEPENDENT_PRINCIPAL_REPLAY_DONE cases=%o classgroup_search=false time=%o\\n",#C40ReplayCases,Cputime(tt);
quit;
'''
(HERE.parent/f'{stem}.m').write_text(source)
print('C40_PRINCIPAL_WITNESS_REPLAYER_PREPARED',stem,'cases',sorted(wanted))
