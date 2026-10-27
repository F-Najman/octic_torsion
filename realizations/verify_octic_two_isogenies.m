/* The same good-reduction bounds from verify_octic_examples prove the
   full target torsion once the displayed subgroup attains orderN. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t0:=Cputime();
load "realizations/data/vanhoeij_octic_input.m";
R<T>:=PolynomialRing(Rationals());tested:=0;
for idx in [1..#OcticExamples] do
 N,fc,xc,yc:=Explode(OcticExamples[idx]);if IsOdd(N) then continue;end if;
 f:=R!fc;K<a>:=NumberField(f);xx:=Evaluate(R!xc,a);yy:=Evaluate(R!yc,a);
 r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx);s:=(xx*yy-yy+1)/(xx*yy);
 b:=r*s*(r-1);c:=s*(r-1);E:=EllipticCurve([1-c,-b,-b,0,0]);P:=E![0,0];
 T2:=(N div 2)*P;assert T2 ne E!0 and 2*T2 eq E!0;
 Et,phi:=IsogenyFromKernel(E,Polynomial([-T2[1],K!1]));
 PP:=phi(P);m:=N div 2;
 assert m*PP eq Et!0 and &and[(m div q)*PP ne Et!0:q in PrimeDivisors(m)];
 Q<U>:=PolynomialRing(K);bb2,bb4,bb6,bb8:=Explode(bInvariants(Et));
 roots:=Roots(4*U^3+bb2*U^2+2*bb4*U+bb6);
 a1,a2,a3,a4,a6:=Explode(aInvariants(Et));
 points:=[Et![z[1],-(a1*z[1]+a3)/2]:z in roots];
 assert &and[z ne Et!0 and 2*z eq Et!0:z in points];
 if IsOdd(m) then
  assert #points eq 1;V:=PP+points[1];
  assert N*V eq Et!0 and &and[(N div q)*V ne Et!0:q in PrimeDivisors(N)];
  inv:=[N];
 else
  assert N eq 40 and #points eq 3;
  V:=[z:z in points|z ne (m div 2)*PP][1];assert V ne (m div 2)*PP;
  inv:=[2,m];
 end if;
 tested+:=1;
 printf "OCTIC_TWO_ISOGENY_CERTIFIED index=%o N=%o image_point_order=%o target_invariants=%o target_order_bound=%o exact=true time=%o\n",idx,N,m,inv,N,Cputime(t0);
end for;
assert tested eq 7;
printf "OCTIC_TWO_ISOGENIES_DONE targets=7 complete=true time=%o\n",Cputime(t0);quit;
