/* Fresh complete D8 good-reduction branch. Main class-group state is separate. */
SetColumns(0);SetMemoryLimit(8*10^9);SetQuitOnError(true);tgood:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");
Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_diamond_hash_cyclic_bc.m");
outgood:=Open("C39/logs/root_fresh_good_39.log","w");
procedure NoteGood(s)
 fprintf outgood,"%o\n",s;Flush(outgood);print s;
end procedure;
X:=MDX1(39,GF(2));C:=Curve(X);cs:=Cusps(X);osraw:=CuspOrbitsQ(X);
os:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in osraw];
assert Genus(X) eq 33 and {*Degree(c):c in cs*} eq {*1^^12,2^^6,12^^2*};
assert {*Degree(o):o in os*} eq {*1^^12,2^^6,12^^2*};
assert &and[#Support(o) eq 1 and Degree(Support(o)[1]) eq Degree(o):o in os];
places:=[NoncuspidalPlaces(X,d):d in [1..8]];
assert [#p:p in places] eq [0,0,0,0,0,12,0,49];
fp:=FrobeniusPolynomial(JOne(39),2);traces:=[];closed:=[];
for n in [1..8] do
 sn:=-n*Coefficient(fp,66-n);
 for j in [1..n-1] do sn-:=Coefficient(fp,66-j)*traces[n-j];end for;
 Append(~traces,Integers()!sn);count:=2^n+1-Integers()!sn;
 for d in Divisors(n) do if d lt n then count-:=d*closed[d];end if;end for;
 assert count mod n eq 0;Append(~closed,count div n);
 assert closed[n] eq #places[n]+#[c:c in cs|Degree(c) eq n];
end for;
NoteGood(Sprintf("C39_D8_GOOD_SETUP genus=33 all_counts=%o noncuspidal_counts=%o source_cusps_all_Q_lift=true time=%o",closed,[#p:p in places],Cputime(tgood)));
load "realizations/data/vanhoeij_octic_input.m";
n8,fc,xc,yc:=Explode(OcticExamples[6]);assert n8 eq 39;
f8:=Polynomial(GF(2),fc);assert IsIrreducible(f8) and Degree(f8) eq 8;
K8<a8>:=ext<GF(2)|f8>;
pp:=Places(C(K8)![Evaluate(Polynomial(GF(2),xc),a8),Evaluate(Polynomial(GF(2),yc),a8)]);
assert #pp eq 1 and Degree(pp[1]) eq 8 and not IsCusp(X,pp[1]);
known8:=Seqset(DiamondOrbit(X,pp[1]));
NoteGood(Sprintf("C39_D8_GOOD_KNOWN_REDUCTION full_diamond_orbit=%o time=%o",#known8,Cputime(tgood)));
tested:=0;known:=0;unknown:=[];repcounts:=[];
for d in [6,8] do
 reps:=PlacesUpToDiamondCyclicHashBC(X,places[d]);Append(~repcounts,#reps);
 NoteGood(Sprintf("C39_D8_GOOD_PARTITION degree=%o places=%o representatives=%o exact_cover=true time=%o",d,#places[d],#reps,Cputime(tgood)));
 for i in [1..#reps] do
  P:=reps[i];
  if d eq 8 and P in known8 then
   known+:=1;NoteGood(Sprintf("C39_D8_GOOD_CASE degree=8 index=%o/%o known_exact_reduction=true time=%o",i,#reps,Cputime(tgood)));continue;
  end if;
  killed:=false;
  for q in [5,7] do
   NoteGood(Sprintf("C39_D8_GOOD_TEST_START degree=%o index=%o/%o q=%o time=%o",d,i,#reps,q,Cputime(tgood)));
   DD:=Divisor(P);AA:=HeckeOperatorCommonZeros(X,q,DD)-q*DiamondOperator(X,q,DD)-DD;
   principal:=IsPrincipal(AA);
   NoteGood(Sprintf("C39_D8_GOOD_WITNESS degree=%o index=%o/%o q=%o principal=%o time=%o",d,i,#reps,q,principal,Cputime(tgood)));
   if not principal then killed:=true;break;end if;
  end for;
  tested+:=1;
  if not killed then Append(~unknown,<d,i,P>);end if;
 end for;
end for;
complete:=#unknown eq 0;
NoteGood(Sprintf("C39_D8_GOOD_DONE representatives=%o nonknown_tested=%o known_representatives=%o unknown=%o complete=%o time=%o",repcounts,tested,known,#unknown,complete,Cputime(tgood)));
if complete then delete outgood;quit;end if;
