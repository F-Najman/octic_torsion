/* Clean post-export group rebuild. Old diagnostic basis is not reused. */
/* Fresh degree8 C48 diagnostic: proven Picard group of X1(2,24)/F5.
   No earlier-degree point classification is used. Retain the live state. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);t48:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "common/class_group_unflattened.m";
load "common/native_group_generators.m";
out48:=Open("C48/logs/clean_run/target_2_24_setup.log","w");
procedure Note48(s)
 fprintf out48,"%o\n",s;Flush(out48);print s;
end procedure;
function ClosedCounts48(fp,p,limit)
 g:=Degree(fp) div 2;assert Degree(fp) eq 2*g;
 aa:=Reverse(Coefficients(fp));assert aa[1] eq 1;
 ss:=[];bb:=[];
 for m in [1..limit] do
  sm:=-m*aa[m+1]-&+[Integers()|aa[j+1]*ss[m-j]:j in [1..m-1]];
  Append(~ss,sm);nm:=p^m+1-sm;
  bm:=nm-&+[Integers()|d*bb[d]:d in Divisors(m)|d lt m];
  assert bm mod m eq 0;Append(~bb,bm div m);
 end for;
 return bb;
end function;
fp48:=FrobeniusPolynomial(JOne(48),5);assert Degree(fp48) eq 74;
all48:=ClosedCounts48(fp48,5,8);cusp48:=[16,12,0,4,0,0,0,0];
nc48:=[all48[j]-cusp48[j]:j in [1..8]];
Note48(Sprintf("C48_D8_SOURCE_NEWTON all=%o cusps=%o noncuspidal=%o JF5=%o time=%o",all48,cusp48,nc48,Evaluate(fp48,1),Cputime(t48)));
Y:=MDX11(2,24,GF(5));C:=Curve(Y);FF:=FunctionField(C);
assert Genus(Y) eq 13 and Y`_E[1] eq 0 and Y`_E[3] eq 0;
b,c:=Explode([FF!f:f in Y`_coordinates]);d:=b*(1-b-c);
E:=EllipticCurve([0,c,0,d,0]);Q:=E![b,b];fs:=[b,d,c^2-4*d];
for k in [1..12] do A:=k*Q;assert A[3] ne 0 and A[1] ne 0;Append(~fs,A[1]);end for;
cs:=Cusps(Y);groups:=AssociativeArray();
for cc in cs do
 sg:=[Valuation(f,cc):f in fs];
 if IsDefined(groups,sg) then Append(~groups[sg],cc);else groups[sg]:=[cc];end if;
end for;
sigs:=Sort(Setseq(Keys(groups)));
blocks:=[&+[Divisor(cc):cc in groups[sg]]:sg in sigs];
degrees:=[Degree(o):o in blocks];
assert Sort(degrees) eq [1:i in [1..8]] cat [2:i in [1..6]] cat [4:i in [1..5]];
baseBlock:=Index(degrees,1);c0:=Support(blocks[baseBlock])[1];
A24:=JH(48,[25]);assert Dimension(A24) eq 13;
fp24:=FrobeniusPolynomial(A24,5);jf24:=Integers()!Evaluate(fp24,1);
all24:=ClosedCounts48(fp24,5,8);
nc24:=[all24[j]-#[cc:cc in cs|Degree(cc) eq j]:j in [1..8]];
Note48(Sprintf("C48_D8_TARGET_MODEL genus=13 cusps=%o rational_block_degrees=%o noncuspidal=%o JF5=%o prospective_index=%o time=%o",[Degree(cc):cc in cs],degrees,nc24,jf24,jf24 div 921600,Cputime(t48)));
SetVerbose("ClassGroup",1);
F0:=AlgorithmicFunctionField(FunctionField(C));FR:=RationalExtensionRepresentation(F0);
iv:=ClassGroupAbelianInvariants(F0:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
assert FR`cld`proven;MDProvenClassGroupGenerators(FR);
G,_,phi:=MDClassGroupUnflattened(F0);assert Invariants(G) eq iv;
assert &*[n:n in Invariants(G)|n ne 0] eq jf24;
m2:=map<DivisorGroup(C)->G|DD:->phi(FunctionFieldDivisor(DD)),g:->CurveDivisor(C,g@@phi)>;
gens:=[m2(o-Degree(o)*Divisor(c0)):o in blocks];
HQ:=sub<G|gens>;assert Invariants(HQ) eq [2,4,4,120,240];
GG,pi:=quo<G|HQ>;
Note48(Sprintf("C48_D8_TARGET_GROUP invariants=%o rational_image=%o quotient=%o time=%o",Invariants(G),Invariants(HQ),Invariants(GG),Cputime(t48)));
cuspclasses:=[m2(Divisor(cc)-Degree(cc)*Divisor(c0)):cc in cs];
comps48:=AssociativeArray();for j in [0..8] do comps48[j]:={GG|};end for;comps48[0]:={GG!0};
for i in [1..#cs] do
 w:=Degree(cs[i]);a:=pi(cuspclasses[i]);
 for j in [w..8] do comps48[j] join:={z+a:z in comps48[j-w]};end for;
end for;
ar48:=Open("C48/data/clean_run/target_2_24_arithmetic.py","w");
fprintf ar48,"proven = True\nsource_frobenius_coefficients = %o\nsource_closed_counts = %o\nsource_noncuspidal_counts = %o\n",Coefficients(fp48),all48,nc48;
fprintf ar48,"target_frobenius_coefficients = %o\ntarget_closed_counts = %o\ntarget_noncuspidal_counts = %o\n",Coefficients(fp24),all24,nc24;
fprintf ar48,"relations = %o\nblock_signatures = %o\nblock_degrees = %o\nbase_block = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G))],sigs,degrees,baseBlock;
fprintf ar48,"rational_cusp_generators = %o\nfinite_cusp_classes = %o\nfinite_cusp_degrees = %o\nquotient_invariants = %o\ncomplement_counts = %o\n",[Eltseq(a):a in gens],[Eltseq(a):a in cuspclasses],[Degree(cc):cc in cs],Invariants(GG),[#comps48[j]:j in [0..8]];
delete ar48;
SetVerbose("ClassGroup",0);
Note48(Sprintf("C48_D8_TARGET_SETUP_DONE complete=true exclusion_claimed=false complement_counts=%o time=%o",[#comps48[j]:j in [0..8]],Cputime(t48)));
