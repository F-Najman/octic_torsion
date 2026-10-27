/* Complete degree7 prefix: source indices348,1102. Every source F5
   rational cusp has global field Q or Q(i), by DEHMZ cusp labels:
   d>2 and 5=1 mod(48/d) forces48/d dividing4; d<=2 gives no fixed
   cusp. Thus13=1 mod4 fixes every degree1 complement and A13 cancels it. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");
load "common/pair_af_place_linear.m";
load "C48/data/fixed_source_inputs.m";
X48d6:=MDX1(48,GF(5));FC48d6:=FunctionField(Curve(X48d6));
F48d6,mp48d6:=AlgorithmicFunctionField(FC48d6);
coords48d6:=[mp48d6(FC48d6.1),mp48d6(FC48d6.2)];
assert Genus(F48d6) eq 37 and 5^3 mod 48 eq 29;
procedure C48HeckeD7Shards()
 child:=false;shard:=1;
 for j in [0..0] do pid:=Fork();if pid eq 0 then child:=true;shard:=j;break;end if;assert pid gt 0;end for;
 idx:=[348,1102][shard+1];dd:=7;q:=13;tt:=Cputime();
 log:=Open(Sprintf("C48/logs/high_degree7_shard%o.log",shard),"w");
 K:=ext<GF(5)|Polynomial(GF(5),C48FieldPolynomials[5])>;
 raw:=C48RepresentativePairs[5][idx];bc:=[&+[K!v[k+1]*K.1^k:k in [0..dd-1]]:v in raw];b,c:=Explode(bc);
 E:=EllipticCurve([1-c,-b,-b,K!0,K!0]);P:=E![0,0];
 assert 48*P eq E!0 and 24*P ne E!0 and 16*P ne E!0;
 fprintf log,"C48_D7_HEC_START source_index=%o q=13 master_sha=%o\n",idx,C48MasterSHA;Flush(log);
 pp:=C48PairAFPlaceLinear(coords48d6,TateNormalForm_xy(E,P));assert Degree(pp) eq dd;
 images:=MDIsogeniesFast(E,P,q);DD:=DivisorGroup(F48d6)!0;
 fprintf log,"C48_D7_ISOGENIES source_index=%o orbits=%o time=%o\n",idx,#images,Cputime(tt);Flush(log);
 for t in images do
  EG,PG,e:=Explode(t);pair:=TateNormalForm_xy(EG,PG);
  yy:=C48PairAFPlaceLinear(coords48d6,pair);mult:=e*dd/Degree(yy);
  assert Denominator(mult) eq 1;DD+:=(Integers()!mult)*Divisor(yy);
  fprintf log,"C48_D7_IMAGE source_index=%o kernel_orbit=%o image_degree=%o time=%o\n",idx,e,Degree(yy),Cputime(tt);Flush(log);
 end for;
 assert Degree(DD) eq (q+1)*dd;
 dp:=C48PairAFPlaceLinear(coords48d6,TateNormalForm_xy(E,q*P));assert Degree(dp) eq dd;
 AA:=DD-q*Divisor(dp)-Divisor(pp);assert Degree(AA) eq 0;
 fprintf log,"C48_D7_PRINCIPALITY_START source_index=%o time=%o\n",idx,Cputime(tt);Flush(log);
 principal:=IsPrincipal(AA);
 fprintf log,"C48_D7_HEC_DONE source_index=%o q=13 principal=%o complete=true time=%o\n",idx,principal,Cputime(tt);delete log;
 print "C48_D7_HEC_DONE",idx,principal,Cputime(tt);
 if child then quit;end if;
end procedure;
C48HeckeD7Shards();quit;
