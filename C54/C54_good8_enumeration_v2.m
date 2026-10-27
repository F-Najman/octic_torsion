/* Fresh exhaustive pure-degree8 image orbits on X1(27)/F5 which
   have rational2-torsion. No degree-seven number-field exclusion used. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(10*10^9);te54:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("C54/hasse_bsgs_count.m");
K<a>:=GF(5^8);F:=GF(5);q:=#K;assert Valuation(q-1,3) eq 1;
seen:={K|};orbitkeys:={};records:=[* *];weighted:=0;jcount:=0;curveCount:=0;rank2curves:=0;properjrecords:=0;
countTime:=0;maxTrials:=0;properHistogram:=AssociativeArray();
for d in [1,2,4,8] do properHistogram[d]:=0;end for;
units:=[k:k in [1..13]|GCD(k,27) eq 1];assert #units eq 9;
function PairKey54(xy)
 return Minimum([Eltseq(xy[1]^(5^i),GF(5)) cat Eltseq(xy[2]^(5^i),GF(5)):i in [0..7]]);
end function;
for jj in K do
 if jj in seen then continue;end if;
 for i in [0..7] do Include(~seen,jj^(5^i));end for;
 jd:=Degree(MinimalPolynomial(jj,F));assert jd in [1,2,4,8];jcount+:=1;
 curves:=Twists(EllipticCurveWithjInvariant(jj));
 for E in curves do
  curveCount+:=1;tcount:=Cputime();
  if jd lt 8 then n:=#E;trials:=0;else n,trials:=EllipticCountHasseBSGS(E);end if;
  countTime+:=Cputime(tcount);maxTrials:=Max(maxTrials,trials);
  if curveCount le 100 then assert n eq #E;end if;
  if n mod 54 ne 0 then continue;end if;
  va:=Valuation(n,3);assert va ge 3;
  psi3:=DivisionPolynomial(E,3);roots3:=Roots(psi3);nr3:=#roots3;assert nr3 in [1,4];
  vb:=nr3 eq 1 select va else va-1;if vb lt 3 then continue;end if;
  found:=false;
  for tr in [1..100] do
   Q:=(n div 3^va)*Random(E);assert 3^vb*Q eq E!0;
   if 3^(vb-1)*Q ne E!0 then found:=true;break;end if;
  end for;assert found;
  P:=3^(vb-3)*Q;assert 27*P eq E!0 and 9*P ne E!0;
  generators:=[P];
  if nr3 eq 4 then
   rank2curves+:=1;tors3:=[];
   for r in roots3 do
    ps:=Points(E,r[1]);assert #ps eq 2 and &and[3*T eq E!0:T in ps];tors3 cat:=ps;
   end for;
   indep:=[T:T in tors3|T ne 9*P and T ne -9*P];assert #indep eq 6;R3:=indep[1];
   generators:=[P+k*R3:k in [0..2]];
  end if;
  b2e,b4e,b6e,b8e:=Explode(bInvariants(E));RR<T>:=PolynomialRing(K);r2:=#Roots(4*T^3+b2e*T^2+2*b4e*T+b6e);assert r2 in [1,3];
  for P27 in generators do
   assert 27*P27 eq E!0 and 9*P27 ne E!0;
   xy:=TateNormalForm_xy(E,P27);dd:=LCM([Degree(MinimalPolynomial(z,F)):z in xy]);
   if dd ne 8 then assert dd in [1,2,4];continue;end if;
   assert &or[Degree(MinimalPolynomial(z,F)) eq 8:z in xy];
   key:=Minimum([PairKey54(TateNormalForm_xy(E,k*P27)):k in units]);
   if key in orbitkeys then continue;end if;Include(~orbitkeys,key);
   Append(~records,[*[[Integers()!c:c in Eltseq(z,F)]:z in xy],n,r2,[Integers()!c:c in Eltseq(jj,F)],nr3*]);
   weighted+:=9*r2;properHistogram[jd]+:=1;if jd lt 8 then properjrecords+:=1;end if;
  end for;
 end for;
 if jcount mod 1000 eq 0 then printf "OCTIC54_GOOD8_ENUM_PROGRESS j_orbits=%o curves=%o representatives=%o weight=%o counttime=%o time=%o\n",jcount,curveCount,#records,weighted,countTime,Cputime(te54);end if;
end for;
assert #seen eq q and #records eq #orbitkeys;
assert jcount eq 5+(5^2-5) div 2+(5^4-5^2) div 4+(5^8-5^4) div 8;
/* The other126 degree8 source places are the degree2 lifts of the
   degree4 X1(27) places with a1+2 splitting two-division cubic. */
assert weighted eq 48087-126;
fieldpoly:=[Integers()!v:v in Eltseq(DefiningPolynomial(K,F))];
Write("C54/data/C54_good8_pointdata.txt",Sprint([*fieldpoly,records*]) : Overwrite:=true);
printf "OCTIC54_GOOD8_ENUM_DONE representatives=%o weighted_source_places=%o other_degree4_lifts=126 j_orbits=%o curves=%o rank2_curves=%o proper_j_records=%o jd_counts=%o max_count_trials=%o complete=true time=%o\n",#records,weighted,jcount,curveCount,rank2curves,properjrecords,[<d,properHistogram[d]>:d in [1,2,4,8]],maxTrials,Cputime(te54);
quit;
