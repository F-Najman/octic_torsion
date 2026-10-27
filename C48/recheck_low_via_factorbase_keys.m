/* No RR descent and no repeated point inversion: evaluate all actual target
   factorbase places once, then match every low source image by exact Tate
   moduli/Frobenius keys. Rechecks the old low vectors in the fresh basis. */
SetQuitOnError(false);
load "C48/data/low_expected_pairs_vectors.m";
function C48NormalizedPairFB(E,P,Q)
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));assert a1 eq 0 and a3 eq 0 and a6 eq 0;
 r:=P[1];x:=Q[1]-r;y:=Q[2];assert P[2] eq 0 and x ne 0 and y ne 0;
 b:=x^3/y^2;c:=(a2+3*r)*x^2/y^2;
 assert (a4+2*a2*r+3*r^2)*x^4/y^4 eq b*(1-b-c);
 return [b,c];
end function;
procedure C48RecheckLowFBKeys()
 assert FR`cld`proven;tt:=Cputime();
 Fcheck,mp:=AlgorithmicFunctionField(FunctionField(Curve(Y)));assert RationalExtensionRepresentation(Fcheck) cmpeq FR;
 coords:=[FR!mp(z):z in Y`_coordinates];
 fb:=Reverse([v`place:v in FR`cld`fb]);ff:=FR`cld`ZnToCl;zn:=Domain(ff);
 A:=DivisorGroup(FR)!FunctionFieldDivisor(Divisor(c0));ap,am:=Support(A);
 assert #ap eq 1 and am eq [1] and Degree(ap[1]) eq 1;
 ja:=Index(fb,ap[1]);assert ja gt 0;ga:=ff(zn.ja);
 table:=AssociativeArray();counts:=[0,0,0,0,0];
 log:=Open("C48/logs/clean_run_v2/low_factorbase_key_replay.log","w");
 ar:=Open("C48/data/clean_run_v2/low_factorbase_key_replay.py","w");
 fprintf ar,"base_factorbase_index = %o\ntarget_table = [\n",ja;
 for j in [1..#fb] do
  pp:=fb[j];dd:=Degree(pp);if dd lt 3 or dd gt 5 then continue;end if;
  if &or[Valuation(z,pp) lt 0:z in coords] then continue;end if;
  bc:=[Evaluate(z,pp):z in coords];b,c:=Explode(bc);d:=b*(1-b-c);
  if d eq 0 or c^2-4*d eq 0 then continue;end if;
  key:=[dd] cat MDM2FrobeniusKey(bc,sourceFields48[dd],5,dd);
  assert not IsDefined(table,key);v:=ff(zn.j)-dd*ga;table[key]:=<j,v>;
  counts[dd]+:=1;
  fprintf ar,"[%o,%o,%o],\n",j,key,Eltseq(v);
 end for;
 assert counts eq [0,0,16,96,672];fprintf ar,"]\nsource_images = [\n";Flush(ar);
 fprintf log,"C48_LOW_FB_TARGET_TABLE_COMPLETE counts=%o places=%o time=%o\n",counts,#Keys(table),Cputime(tt);Flush(log);
 tested:=0;
 for dd in [3..5] do
  K:=sourceFields48[dd];pairs:=C48LowExpectedAllPairs[dd-2];expected:=C48LowExpectedAllClasses[dd-2];
  assert #pairs eq #expected;
  for i in [1..#pairs] do
   raw:=pairs[i];bc:=[&+[K!v[k+1]*K.1^k:k in [0..dd-1]]:v in raw];b,c:=Explode(bc);
   EE:=EllipticCurve([1-c,-b,-b,K!0,K!0]);E2,L2:=C48TwoIsogenousModuli(EE,EE![0,0]);
   targetbc:=C48NormalizedPairFB(E2,L2`P,L2`Q);
   assert targetbc eq DerickxNormalForm_bc(E2,L2`P,L2`Q);
   key:=[dd] cat MDM2FrobeniusKey(targetbc,K,5,dd);
   assert IsDefined(table,key);j,v:=Explode(table[key]);assert Eltseq(v) eq expected[i];tested+:=1;
   fprintf ar,"[%o,%o,%o,%o],\n",dd,i,j,Eltseq(v);
  end for;
  fprintf log,"C48_LOW_FB_SOURCE_DEGREE_PASS degree=%o tested=%o time=%o\n",dd,#pairs,Cputime(tt);Flush(log);
 end for;
 assert tested eq 736;fprintf ar,"]\ncomplete = True\n";delete ar;
 fprintf log,"C48_LOW_FB_KEY_REPLAY_DONE target_places=784 source_images=736 actual_vectors_exact=true complete=true time=%o\n",Cputime(tt);delete log;
 print "C48_LOW_FB_KEY_REPLAY_DONE",tested,Cputime(tt);
end procedure;
C48RecheckLowFBKeys();
