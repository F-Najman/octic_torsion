/* Exact same source Hecke operator, entirely in the algorithmic function
   field. Avoid expensive plane-curve-place conversions. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");
load "common/pair_af_place_linear.m";
load "C48/data/fixed_source_inputs.m";
load "C48/data/low_hecke_inputs.m";
procedure C48HeckeLowAF()
 tt:=Cputime();X:=MDX1(48,GF(5));CC:=Curve(X);assert Genus(X) eq 37;
 FC:=FunctionField(CC);F,mp:=AlgorithmicFunctionField(FC);coords:=[mp(FC.1),mp(FC.2)];
 fields:=AssociativeArray();
 for dd in [3..5] do fields[dd]:=ext<GF(5)|Polynomial(GF(5),C48FieldPolynomials[dd-2])>;end for;
 log:=Open("C48/logs/low_hecke_af_xy.log","w");
 ar:=Open("C48/data/low_hecke_af_xy_results.py","w");
 fprintf ar,"input_arithmetic_sha = \"%o\"\ntests = [\n",C48LowArithmeticSHA;
 cache:=AssociativeArray();survivors:=[];
 for i in [1..#C48LowRowIDs] do
  degrees:=C48LowSupports[i];
  fprintf log,"C48_LOW_AF_CASE_START row=%o support=%o time=%o\n",C48LowRowIDs[i],degrees,Cputime(tt);Flush(log);
  killed:=false;witness:=0;
  for q in [7,11,13] do
   AA:=DivisorGroup(F)!0;
   for j in [1..#degrees] do
    dd:=degrees[j];raw:=C48LowPairs[i][j];key:=[q,dd] cat &cat raw;
    if not IsDefined(cache,key) then
     K:=fields[dd];bc:=[&+[K!v[k+1]*K.1^k:k in [0..dd-1]]:v in raw];b,c:=Explode(bc);
     EE:=EllipticCurve([1-c,-b,-b,K!0,K!0]);PP:=EE![0,0];
     assert 48*PP eq EE!0 and 24*PP ne EE!0 and 16*PP ne EE!0;
     pp:=C48PairAFPlaceLinear(coords,TateNormalForm_xy(EE,PP));assert Degree(pp) eq dd;
     fprintf log,"C48_LOW_AF_INPUT q=%o degree=%o time=%o\n",q,dd,Cputime(tt);Flush(log);
     images:=MDIsogeniesFast(EE,PP,q);DD:=DivisorGroup(F)!0;
     fprintf log,"C48_LOW_AF_ISOGENIES q=%o degree=%o orbits=%o time=%o\n",q,dd,#images,Cputime(tt);Flush(log);
     for t in images do
      EG,PG,e:=Explode(t);pair:=TateNormalForm_xy(EG,PG);
      yy:=C48PairAFPlaceLinear(coords,pair);mult:=e*dd/Degree(yy);
      assert Denominator(mult) eq 1;DD+:=(Integers()!mult)*Divisor(yy);
      fprintf log,"C48_LOW_AF_IMAGE q=%o source_degree=%o kernel_orbit=%o image_degree=%o time=%o\n",q,dd,e,Degree(yy),Cputime(tt);Flush(log);
     end for;
     assert Degree(DD) eq (q+1)*dd;
     dp:=C48PairAFPlaceLinear(coords,TateNormalForm_xy(EE,q*PP));assert Degree(dp) eq dd;
     cache[key]:=DD-q*Divisor(dp)-Divisor(pp);assert Degree(cache[key]) eq 0;
    end if;
    AA+:=cache[key];
   end for;
   principal:=IsPrincipal(AA);
   fprintf log,"C48_LOW_AF_TEST row=%o support=%o indices=%o q=%o principal=%o time=%o\n",C48LowRowIDs[i],degrees,C48LowIndices[i],q,principal,Cputime(tt);Flush(log);
   if not principal then killed:=true;witness:=q;break;end if;
  end for;
  if not killed then Append(~survivors,C48LowRowIDs[i]);end if;
  fprintf ar,"[%o,%o,%o,%o],\n",C48LowRowIDs[i],degrees,C48LowIndices[i],witness;Flush(ar);
 end for;
 fprintf ar,"]\nsurvivors = %o\ncomplete = True\n",survivors;delete ar;
 fprintf log,"C48_LOW_AF_DONE tested=%o remaining=%o complete=true time=%o\n",#C48LowRowIDs,#survivors,Cputime(tt);delete log;
 print "C48_LOW_AF_DONE",#C48LowRowIDs,#survivors,Cputime(tt);
end procedure;
C48HeckeLowAF();quit;
