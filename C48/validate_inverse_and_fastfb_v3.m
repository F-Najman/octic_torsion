/* Fresh native inverse controls and exact lookup acceleration benchmark.
   This file runs only against the retained clean proved group. */
load "common/class_group_direct_RR_checked_v4_fastfb.m";
procedure C48ValidateInverseAndFastFB()
 log:=Open("C48/logs/clean_run_v2/inverse_and_fastfb_v3.log","w");
 ar:=Open("C48/data/clean_run_v2/inverse_fixtures_v3.py","w");
 fprintf ar,"fixtures = [\n";checked:=[];tt:=Cputime();
 for dd in [3..8] do
  found:=false;
  for idx in [1..Minimum(20,#sourceRepBC48[dd])] do
   b,c:=Explode(sourceRepBC48[dd][idx]);K:=Parent(b);
   EE:=EllipticCurve([1-c,-b,-b,K!0,K!0]);E2,L2:=C48TwoIsogenousModuli(EE,EE![0,0]);
   bc:=DerickxNormalForm_bc(E2,L2`P,L2`Q);
   if Degree(MinimalPolynomial(bc[1],GF(5))) ne dd then continue;end if;
   yy:=C48PairPlaceLinear(Y`_coordinates,bc);zz:=MDPlaceCommonZeros(Y`_coordinates,bc);
   assert yy eq zz and Degree(yy) eq dd;
   Append(~checked,yy);
   fprintf ar,"[%o,%o,%o,%o,%o],\n",dd,idx,[Eltseq(z,GF(5)):z in sourceRepBC48[dd][idx]],[Eltseq(z,GF(5)):z in bc],1;
   fprintf log,"C48_NATIVE_INVERSE_PASS degree=%o source_index=%o weight=1 time=%o\n",dd,idx,Cputime(tt);Flush(log);
   found:=true;break;
  end for;
  assert found;
 end for;
 found:=false;K4:=sourceFields48[4];k:=GF(5);
 for idx in [1..#sourceRepBC48[8]] do
  b,c:=Explode(sourceRepBC48[8][idx]);K:=Parent(b);
  EE:=EllipticCurve([1-c,-b,-b,K!0,K!0]);E2,L2:=C48TwoIsogenousModuli(EE,EE![0,0]);bc:=DerickxNormalForm_bc(E2,L2`P,L2`Q);
  ee:=LCM([Degree(MinimalPolynomial(z,k)):z in bc]);
  if ee ne 4 then continue;end if;
  gen:=K!0;
  for trygen in [bc[1],bc[2]] cat [bc[1]+a*bc[2]:a in [1..4]] do
   pol:=MinimalPolynomial(trygen,k);if Degree(pol) eq 4 then gen:=trygen;break;end if;
  end for;
  assert Degree(pol) eq 4;
  mat:=Matrix(k,[Eltseq(gen^i,k):i in [0..3]]);cc:=[Eltseq(Solution(mat,Vector(k,Eltseq(z,k)))):z in bc];
  roots:=Roots(PolynomialRing(K4)!pol);assert #roots eq 4;rt:=roots[1][1];
  bc4:=[&+[K4!aa[i+1]*rt^i:i in [0..3]]:aa in cc];
  yy:=C48PairPlaceLinear(Y`_coordinates,bc);zz:=MDPlaceCommonZeros(Y`_coordinates,bc4);
  assert yy eq zz and Degree(yy) eq 4 and 8 div Degree(yy) eq 2;
  Append(~checked,yy);
  fprintf ar,"[8,%o,%o,%o,2,%o],\n",idx,[Eltseq(z,k):z in sourceRepBC48[8][idx]],[Eltseq(z,k):z in bc],[Eltseq(z,k):z in bc4];
  fprintf log,"C48_NATIVE_INVERSE_PASS degree=8 source_index=%o target_degree=4 weight=2 time=%o\n",idx,Cputime(tt);Flush(log);
  found:=true;break;
 end for;
 assert found;fprintf ar,"]\ncomplete = True\n";delete ar;
 Gfast,_,phiFast,wFast:=MDClassGroupDirectRRCheckedV4FastFB(F0,"C48/logs/clean_run_v2/fastfb_control_identities.log");
 Gbase,_,phiBase,wBase:=MDClassGroupDirectRRCheckedV4(F0,"C48/logs/clean_run_v2/baseline_control_identities.log");
 assert Gfast cmpeq G and Gbase cmpeq G;
 for j in [1..#checked] do
  DD:=FunctionFieldDivisor(Divisor(checked[j])-Degree(checked[j])*Divisor(c0));
  t:=Cputime();vbase:=phiBase(DD);tb:=Cputime(t);
  t:=Cputime();vfast:=phiFast(DD);tf:=Cputime(t);assert vbase eq vfast;
  fprintf log,"C48_FASTFB_CONTROL_PASS index=%o degree=%o baseline=%o accelerated=%o vector=%o\n",j,Degree(checked[j]),tb,tf,Eltseq(vfast);Flush(log);
 end for;
 fprintf wFast,"]\n";delete wFast;fprintf wBase,"]\n";delete wBase;
 fprintf log,"C48_INVERSE_AND_FASTFB_COMPLETE native=7 group_comparisons=7 exact=true time=%o\n",Cputime(tt);delete log;
 print "C48_INVERSE_AND_FASTFB_COMPLETE";
end procedure;
C48ValidateInverseAndFastFB();
