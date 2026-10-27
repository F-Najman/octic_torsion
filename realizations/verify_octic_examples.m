/* Exact full-torsion certificates for all10 archived vanHoeij octic places.
   No earlier-degree conclusion or dependency is used. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
load "realizations/data/vanhoeij_octic_input.m";
out8:=Open("realizations/logs/verify_octic_examples.log","w");
procedure Note8(s)
 fprintf out8,"%o\n",s;Flush(out8);print s;
end procedure;
R<T>:=PolynomialRing(Rationals());
function FromXY(xx,yy)
 r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx);
 s:=(xx*yy-yy+1)/(xx*yy);b:=r*s*(r-1);c:=s*(r-1);
 return EllipticCurve([1-c,-b,-b,0,0]);
end function;
function EvalMod(f,a)
 L:=Parent(a);
 return &+[L!Numerator(Coefficient(f,i))/(L!Denominator(Coefficient(f,i)))*a^i:i in [0..Degree(f)]];
end function;
for idx in [1..#OcticExamples] do
 N,fc,xc,yc:=Explode(OcticExamples[idx]);
 f:=R!fc;xp:=R!xc;yp:=R!yc;
 assert Degree(f) eq 8 and IsIrreducible(f) and IsMonic(f);
 K<a>:=NumberField(f);xx:=Evaluate(xp,a);yy:=Evaluate(yp,a);
 E:=FromXY(xx,yy);assert Discriminant(E) ne 0;P:=E![0,0];
 assert N*P eq E!0 and &and[(N div q)*P ne E!0:q in PrimeDivisors(N)];
 if N eq 33 then
  assert xx^4-3*xx^3+3*xx^2-3*xx+1 eq 0;
  assert yy^2+(xx^3-3*xx^2+xx-3)*yy+2*xx^2-xx+2 eq 0;
 end if;
 Q<U>:=PolynomialRing(K);bb2,bb4,bb6,bb8:=Explode(bInvariants(E));
 rt:=Roots(4*U^3+bb2*U^2+2*bb4*U+bb6);n2:=#rt;assert n2 in [0,1,3];
 if IsEven(N) then
  assert n2 gt 0;lower:=n2 eq 3 select 2*N else N;
  inv:=n2 eq 3 select [2,N] else [N];
 else
  lower:=n2 eq 3 select 4*N else (n2 eq 1 select 2*N else N);
  inv:=n2 eq 3 select [2,2*N] else (n2 eq 1 select [2*N] else [N]);
 end if;
 Note8(Sprintf("OCTIC_EXACT_POINT index=%o N=%o field_degree=8 j_degree=%o rational_two_roots=%o torsion_lower=%o anticipated_invariants=%o time=%o",idx,N,Degree(MinimalPolynomial(jInvariant(E))),n2,lower,inv,Cputime(t0)));
 bound:=0;used:=0;df:=Discriminant(f);
 for p in PrimesUpTo(97) do
  if p eq 2 then continue;end if;
  if Numerator(df) mod p eq 0 or Denominator(df) mod p eq 0 then continue;end if;
  if &or[Denominator(z) mod p eq 0:z in Coefficients(f) cat Coefficients(xp) cat Coefficients(yp)] then continue;end if;
  fp:=Polynomial(GF(p),f);assert IsSquarefree(fp);
  for fg in Factorization(fp) do
   gg:=fg[1];d:=Degree(gg);
   if p^d gt 10^12 then continue;end if;
   if d eq 1 then L:=GF(p);aa:=-Coefficient(gg,0)/Coefficient(gg,1);
   else L<aa>:=ext<GF(p)|gg>;end if;
   try
    xe:=EvalMod(xp,aa);ye:=EvalMod(yp,aa);Ee:=FromXY(xe,ye);
    if Discriminant(Ee) eq 0 then continue;end if;
    PP:=Ee![0,0];assert N*PP eq Ee!0;
    nn:=#Ee;assert nn mod lower eq 0;
    bound:=GCD(bound,nn);used+:=1;
    Note8(Sprintf("OCTIC_GOOD_REDUCTION index=%o p=%o residue_degree=%o order=%o gcd=%o time=%o",idx,p,d,nn,bound,Cputime(t0)));
   catch er
    Note8(Sprintf("OCTIC_DISPLAYED_MODEL_UNAVAILABLE index=%o p=%o residue_degree=%o",idx,p,d));
   end try;
   if bound eq lower then break;end if;
  end for;
  if bound eq lower then break;end if;
 end for;
 Note8(Sprintf("OCTIC_TORSION_ORDER_BOUND index=%o N=%o lower=%o upper=%o primes_used=%o",idx,N,lower,bound,used));
 assert bound eq lower;
 Note8(Sprintf("OCTIC_FULL_TORSION_CERTIFIED index=%o invariants=%o exact=true time=%o",idx,inv,Cputime(t0)));
end for;
Note8(Sprintf("OCTIC_VANHOEIJ_EXAMPLES_DONE count=%o complete=true time=%o",#OcticExamples,Cputime(t0)));
delete out8;quit;
