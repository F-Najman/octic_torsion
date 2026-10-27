SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);
if assigned N then N:=StringToInteger(N);else N:=115;end if;
D:=FullDirichletGroup(N);bad:=[];allchi:=[];total:=0;
for chi in GaloisConjugacyRepresentatives(D) do
 if IsOdd(chi) then continue;end if;
 tt:=Cputime();
 chimin:=MinimalBaseRingCharacter(chi);
 S:=CuspidalSubspace(ModularSymbols(chimin,2,1));
 dd:=Dimension(S);
 printf "WINDING_SPACE_READY N=%o order=%o dimension=%o time=%o\n",N,Order(chi),dd,Cputime(tt);
 if N eq 529 and Order(chi) eq 253 then
  vambient:=Vector(Eltseq(WindingElement(S)))*ProjectionMatrix(S);
  v:=Vector(Coordinates(VectorSpace(S),vambient));
  A:=HeckeOperator(S,2);
  printf "WINDING_T2_READY N=%o order=%o time=%o\n",N,Order(chi),Cputime(tt);
  K:=BaseRing(A);
  coeffs:=&cat[Eltseq(K!x):x in Eltseq(v) cat Eltseq(A)];
  pp:=0;
  for mult in [2..200] do
   cand:=253*mult+1;
   if IsPrime(cand) and &and[Denominator(c) mod cand ne 0:c in coeffs] then pp:=cand;break;end if;
  end for;
  assert pp gt 0;F:=GF(pp);
  printf "WINDING_GOOD_PRIME p=%o all_denominators_invertible=true\n",pp;
  roots:=Roots(PolynomialRing(F)!DefiningPolynomial(K));assert #roots gt 0;r:=roots[1][1];
  red:=func<x|&+[F!c*r^(i-1):i->c in Eltseq(K!x)]>;
  vm:=Vector(F,[red(x):x in Eltseq(v)]);
  Am:=Matrix(F,dd,dd,[red(x):x in Eltseq(A)]);
  vs:=[];u:=vm;
  for k in [1..dd] do Append(~vs,u);u:=u*Am;end for;
  ww:=Rank(Matrix(vs));
  if assigned KrylovFile then
   Write(KrylovFile,Sprintf("p:=%o;d:=%o;v:=%o;A:=%o;",#F,dd,[Integers()!x:x in Eltseq(vm)],[Integers()!x:x in Eltseq(Am)]) : Overwrite:=true);
  end if;
  printf "WINDING_KRYLOV N=%o order=%o span_rank=%o dimension=%o verified=%o time=%o\n",N,Order(chi),ww,dd,ww eq dd,Cputime(tt);
  assert ww eq dd;
 else
  ww:=Dimension(WindingSubmodule(S));
 end if;
 printf "WINDING N=%o order=%o conductor=%o relative_dimension=%o winding_dimension=%o coefficient_degree=%o time=%o\n",N,Order(chi),Conductor(chi),dd,ww,Degree(BaseRing(S)),Cputime(tt);
 total+:=dd*EulerPhi(Order(chi));
 Append(~allchi,chi);if dd ne ww then Append(~bad,chi);end if;
end for;
I:=EulerPhi(N)*N*&*[Rationals()|1+1/p[1]:p in Factorization(N)]/2;
g:=Integers()!(1+I/12-&+[EulerPhi(d)*EulerPhi(N div d):d in Divisors(N)]/4);
assert total eq g;
admissible:=[a:a in [2..N-2]|GCD(a,N) eq 1 and &and[Evaluate(chi,a) eq 1:chi in bad]];
printf "WINDING_FINAL N=%o genus=%o total_dimension=%o bad_orders=%o admissible=%o\n",N,g,total,[Order(chi):chi in bad],admissible;
print "WINDING_CHARACTER_PROBE_DONE";
quit;
