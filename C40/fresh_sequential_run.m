/* Fresh automatic complete C40 run; no retained state or manual interruption. */
/* Fresh degree8 X1(40)/F3 proven setup without inverse-generator flattening.
   All finite cusp/Q-orbit data are recorded, including non-Q lifts. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);t40:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "common/class_group_unflattened.m";
load "common/native_group_generators.m";
out40:=Open("C40/logs/fresh_sequential/classify_octic_40_setup_unflattened.log","w");
procedure Note40(s)
 fprintf out40,"%o\n",s;Flush(out40);print s;
end procedure;
fp:=FrobeniusPolynomial(JOne(40),3);jf3:=Integers()!Evaluate(fp,1);
X:=MDX1(40,GF(3));C:=Curve(X);cs:=Cusps(X);osraw:=CuspOrbitsQ(X);
os:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in osraw];
assert Genus(X) eq 25 and jf3 eq 31579088486400;
assert #os eq 21 and #Seqset([CuspSignature(X,Support(o)[1]):o in os]) eq 21;
assert &+[Degree(c):c in cs] eq 48;
c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
Note40(Sprintf("C40_D8_MODEL_READY genus=25 finite_cusp_degrees=%o rational_orbit_degrees=%o JF3_order=%o time=%o",[Degree(c):c in cs],[Degree(o):o in os],jf3,Cputime(t40)));
SetVerbose("ClassGroup",1);
F40raw:=AlgorithmicFunctionField(FunctionField(C));FR40raw:=RationalExtensionRepresentation(F40raw);
iv40:=ClassGroupAbelianInvariants(F40raw:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
assert FR40raw`cld`proven;
MDProvenClassGroupGenerators(FR40raw);
G,mraw40,phi40:=MDClassGroupUnflattened(F40raw);
assert Invariants(G) eq iv40;
m2:=map<DivisorGroup(C)->G|D:->phi40(FunctionFieldDivisor(D)),g:->CurveDivisor(C,g@@phi40)>;
m1:=Inverse(m2);
Note40(Sprintf("C40_D8_PROVEN_UNFLATTENED_MAP_DONE generator_flattening_omitted=true time=%o",Cputime(t40)));
assert &*[n:n in Invariants(G)|n ne 0] eq jf3;
qc:=[m2(o-Degree(o)*Divisor(c0)):o in os];HQ:=sub<G|qc>;
assert Invariants(HQ) eq [2,2,2,8,120,895440] and 4592*#HQ eq jf3;
raw40:=Open("C40/data/fresh_sequential/proven_relation_state.py","w");
fprintf raw40,"proven = True\norder = 31579088486400\n";
fprintf raw40,"abstract_relations = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G))];
fprintf raw40,"raw_relations = %o\n",[Eltseq(r):r in FR40raw`cld`rels];
fprintf raw40,"factorbase_degrees = %o\n",[Degree(v`place):v in FR40raw`cld`fb];
fprintf raw40,"factorbase_abstract_classes = %o\n",[Eltseq(FR40raw`cld`ZnToCl(Domain(FR40raw`cld`ZnToCl).i)):i in [1..Ngens(Domain(FR40raw`cld`ZnToCl))]];
delete raw40;
Note40(Sprintf("C40_D8_GROUP_READY invariants=%o rational_image=%o time=%o",Invariants(G),Invariants(HQ),Cputime(t40)));
places:=[NoncuspidalPlaces(X,d):d in [1..8]];
assert [#p:p in places] eq [0,0,0,32,72,80,408,708];
Note40(Sprintf("C40_D8_SETUP_DONE noncuspidal_counts=%o continuation_pending=true time=%o",[#p:p in places],Cputime(t40)));

/* Fresh complete restart; both slower prefixes retained as incomplete diagnostics. */
/* Exact degree-g direct RR map with both finite and infinite ideal identities. */
load "common/class_group_direct_RR_checked_v4.m";
Gold40:=G;
G,mraw40RR,phi40RR,freshWitnesses40:=MDClassGroupDirectRRCheckedV4(AlgorithmicFunctionField(FunctionField(C)),"C40/logs/fresh_sequential/direct_RR_principal_identities.log");
assert G cmpeq Gold40;
m2:=map<DivisorGroup(C)->G|DD:->phi40RR(FunctionFieldDivisor(DD)),g:->CurveDivisor(C,g@@phi40RR)>;
m1:=Inverse(m2);
assert [m2(o-Degree(o)*Divisor(c0)):o in os] eq qc;

assert assigned G and assigned HQ and assigned places and assigned qc;
assert Genus(X) eq 25 and [#p:p in places] eq [0,0,0,32,72,80,408,708];
SetVerbose("ClassGroup",0);
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_diamond_hash_cyclic_bc.m");
out40c:=Open("C40/logs/fresh_sequential/classify_octic_40_continuation_direct_RR.log","w");
procedure Note40C(s)
 fprintf out40c,"%o\n",s;Flush(out40c);print s;
end procedure;
load "C40/known_40_reductions.m";
known40classes:={m2(dd40-8*Divisor(c0)):dd40 in known40divisors};
assert #known40classes eq 4 and &and[a40 in HQ:a40 in known40classes];
Note40C(Sprintf("C40_D8_KNOWN_DIVISORS classes=4 weighted_reduction_patterns=%o time=%o",known40sources,Cputime(t40)));

function EffectiveClasses40(weights,gens,limit,G)
 cc:=AssociativeArray();for j in [0..limit] do cc[j]:={G|};end for;cc[0]:={G!0};
 for i in [1..#weights] do
  w:=weights[i];a:=gens[i];
  if w gt limit then continue;end if;
  for j in [w..limit] do cc[j] join:={z+a:z in cc[j-w]};end for;
 end for;
 return cc;
end function;
cs40classes:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];
cc40:=EffectiveClasses40([Degree(c):c in cs],cs40classes,8,G);
qq40:=EffectiveClasses40([Degree(o):o in os],qc,8,G);
assert &and[a in HQ:a in qq40[8]];
allcusp40members:={a:a in cc40[8]|a in HQ};
unknowncusp40:=allcusp40members diff qq40[8];
Note40C(Sprintf("C40_D8_CUSP_CLASSES finite_counts=%o rational_counts=%o finite_degree8_in_H=%o unknown_cuspal_classes=%o time=%o",[#cc40[j]:j in [0..8]],[#qq40[j]:j in [0..8]],#allcusp40members,#unknowncusp40,Cputime(t40)));

ar40:=Open("C40/data/fresh_sequential/class_arithmetic_direct_RR.py","w");
fprintf ar40,"relations = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G))];
fprintf ar40,"rational_cusp_generators = %o\n",[Eltseq(a):a in qc];
fprintf ar40,"rational_cusp_degrees = %o\n",[Degree(o):o in os];
fprintf ar40,"finite_cusp_classes = %o\n",[Eltseq(a):a in cs40classes];
fprintf ar40,"finite_cusp_degrees = %o\n",[Degree(c):c in cs];
fprintf ar40,"finite_effective_class_counts = %o\n",[#cc40[j]:j in [0..8]];
fprintf ar40,"rational_effective_class_counts = %o\n",[#qq40[j]:j in [0..8]];
fprintf ar40,"finite_cuspal_degree8_members = %o\n",#allcusp40members;
fprintf ar40,"unknown_cuspal_classes = %o\n",#unknowncusp40;
fprintf ar40,"known_octic_classes = %o\n",[Eltseq(a):a in known40classes];
fprintf ar40,"known_reduction_prime_patterns = %o\n",known40sources;
fprintf ar40,"frobenius_coefficients = %o\n",Coefficients(fp);

reps40:=AssociativeArray();
for d40 in [4..8] do
 reps40[d40]:=PlacesUpToDiamondCyclicHashBC(X,places[d40]);
 Note40C(Sprintf("C40_D8_EXACT_DIAMOND_PARTITION degree=%o places=%o representatives=%o time=%o",d40,#places[d40],#reps40[d40],Cputime(t40)));
end for;
fprintf ar40,"representative_counts = %o\n",[#reps40[d40]:d40 in [4..8]];
// All quartic classes are cached: pairs reduce only their first support.
quartic40classes:=[];
for i40 in [1..32] do
 Append(~quartic40classes,m2(Divisor(places[4][i40])-4*Divisor(c0)));
 Note40C(Sprintf("C40_D8_QUARTIC_CLASS_CACHE index=%o/32 time=%o",i40,Cputime(t40)));
end for;
fprintf ar40,"quartic_classes = %o\n",[Eltseq(a):a in quartic40classes];
fprintf ar40,"quartic_representative_indices = %o\n",[Index(places[4],P):P in reps40[4]];
fprintf ar40,"tests = [\n";
candidates40:={G|};singlecounts40:=[];singlepasses40:=[];
for d40 in [4..8] do
 passes40:=0;
 for i40 in [1..#reps40[d40]] do
  P40:=reps40[d40][i40];
  a40:=d40 eq 4 select quartic40classes[Index(places[4],P40)] else m2(Divisor(P40)-d40*Divisor(c0));
  found40:={a40+z40:z40 in cc40[8-d40]|a40+z40 in HQ};
  candidates40 join:=found40;passes40+:=#found40;
  fprintf ar40,"[%o,%o,%o,%o],\n",d40,i40,Eltseq(a40),#found40;
  Note40C(Sprintf("C40_D8_SINGLE_CLASS_TEST degree=%o index=%o/%o admitted_total_classes=%o time=%o",d40,i40,#reps40[d40],#found40,Cputime(t40)));
 end for;
 Append(~singlecounts40,#reps40[d40]);Append(~singlepasses40,passes40);
end for;
fprintf ar40,"]\nquartic_pair_tests = [\n";
paircount40:=0;pairpasses40:=0;
for P40 in reps40[4] do
 i40:=Index(places[4],P40);
 for j40 in [1..32] do
  a40:=quartic40classes[i40]+quartic40classes[j40];member40:=a40 in HQ;
  paircount40+:=1;if member40 then pairpasses40+:=1;Include(~candidates40,a40);end if;
  fprintf ar40,"[%o,%o,%o,%o],\n",i40,j40,Eltseq(a40),member40 select 1 else 0;
  Note40C(Sprintf("C40_D8_QUARTIC_PAIR first=%o second=%o in_rational_image=%o time=%o",i40,j40,member40,Cputime(t40)));
 end for;
end for;
assert paircount40 eq #reps40[4]*32;
unknown40:=candidates40 diff known40classes diff qq40[8];
complete40:=#unknown40 eq 0 and #unknowncusp40 eq 0;
fprintf ar40,"]\nunknown_noncuspidal_classes = %o\n",#unknown40;
fprintf ar40,"complete = %o\n",complete40 select "True" else "False";delete ar40;
Note40C(Sprintf("C40_OCTIC_CLASSIFICATION_DONE single_representatives=%o single_admitted_classes=%o quartic_pairs=%o pair_members=%o candidate_classes=%o known_octic_classes=4 unknown_noncuspidal_classes=%o unknown_cuspal_classes=%o complete=%o time=%o",singlecounts40,singlepasses40,paircount40,pairpasses40,#candidates40,#unknown40,#unknowncusp40,complete40,Cputime(t40)));
fprintf freshWitnesses40,"]\n";delete freshWitnesses40;
if not complete40 then
 // Exact short Q-cuspidal representatives permit independent-prime RR refinement.
 // The kernel relations are globally principal because the rational reduction is injective.
 FF40:=FreeAbelianGroup(#qc);hh40:=hom<FF40->G|qc>;
 KK40:=Matrix(Integers(),[Eltseq(FF40!v):v in Generators(Kernel(hh40))]);
 LL40:=LatticeWithBasis(LLL(KK40));short40:=[];
 targets40:=Sort([Eltseq(a):a in (unknown40 join unknowncusp40)]);
 for tv40 in targets40 do
  a40:=G!tv40;v40:=Eltseq(a40@@hh40);w40:=ClosestVector(LL40,Vector(Rationals(),v40));
  sv40:=[v40[j]-Integers()!w40[j]:j in [1..#qc]];assert hh40(FF40!sv40) eq a40;
  Append(~short40,sv40);
 end for;
 shortout40:=Open("C40/data/fresh_sequential/residual_short_cusp_lifts.txt","w");
 fprintf shortout40,"%o\n",[*40,8,[CuspSignature(X,Support(o)[1]):o in os],[Degree(o):o in os],CuspSignature(X,c0),short40,targets40*];delete shortout40;
 Note40C(Sprintf("C40_D8_RESIDUAL_SHORT_LIFTS count=%o independent_prime_RR_pending=true time=%o",#short40,Cputime(t40)));
end if;
if complete40 then delete out40c;delete out40;quit;end if;
