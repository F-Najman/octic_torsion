/* Fresh degree4 necessary finite-image test, used only as a lemma forD8.
   No earlier-degree exclusion is transferred. */
if not assigned C230Out then C230Out:="C2x30/logs/certified_run_v2";end if;
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "common/class_group_unflattened.m";
load "common/native_group_generators.m";
X:=MDX1(30,GF(7));C:=Curve(X);assert Genus(X) eq 9;
os:=CuspOrbitsQ(X);cs:=Cusps(X);c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
F0,mp:=AlgorithmicFunctionField(FunctionField(C));F:=RationalExtensionRepresentation(F0);
iv:=ClassGroupAbelianInvariants(F0:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
assert F`cld`proven;MDProvenClassGroupGenerators(F);G,_,phi:=MDClassGroupUnflattened(F0);
gens:=[phi(FunctionFieldDivisor(o-Degree(o)*Divisor(c0))):o in os];H:=sub<G|gens>;assert Invariants(H) eq [4,8160];GG,pi:=quo<G|H>;
cuspclasses:=[phi(FunctionFieldDivisor(Divisor(c)-Degree(c)*Divisor(c0))):c in cs];
cc:=AssociativeArray();for d in [0..4] do cc[d]:={GG|};end for;cc[0]:={GG!0};
for i in [1..#cs] do w:=Degree(cs[i]);for d in [w..4] do cc[d] join:={z+pi(cuspclasses[i]):z in cc[d-w]};end for;end for;
qq:=AssociativeArray();for d in [0..4] do qq[d]:={G|};end for;qq[0]:={G!0};
for i in [1..#os] do w:=Degree(os[i]);for d in [w..4] do qq[d] join:={z+gens[i]:z in qq[d-w]};end for;end for;
log:=Open(C230Out cat "/quartic_filter.log","w");
ar:=Open(C230Out cat "/quartic_filter.py","w");
fprintf ar,"relations = %o\nrational_cusp_generators = %o\nfinite_cusp_degrees = %o\nfinite_cusp_classes = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(g):g in gens],[Degree(c):c in cs],[Eltseq(c):c in cuspclasses];
fprintf log,"X130_D4_GROUP G=%o H=%o quotient=%o rational_cusp_degree4_classes=%o time=%o\n",Invariants(G),Invariants(H),Invariants(GG),#qq[4],Cputime(tt);Flush(log);
jfun:=F!mp(jInvariant(EllipticCurve(X`_E)));classes:=AssociativeArray();total:=0;passed:=[];
fprintf ar,"place_classes = {}\ntests = [\n";
for d in [2..4] do
 ps:=[p:p in Places(F,d)|Valuation(jfun,p) ge 0];classes[d]:=[phi(Divisor(p)-d*FunctionFieldDivisor(Divisor(c0))):p in ps];
 fprintf log,"X130_D4_POINTS degree=%o noncuspidal=%o time=%o\n",d,#ps,Cputime(tt);Flush(log);
 for i in [1..#ps] do
  a:=classes[d][i];pass:=-pi(a) in cc[4-d];known:=false;
  if pass and d eq 4 then known:=a in qq[4];end if;
  if pass and not known then Append(~passed,<d,i,Eltseq(a)>);end if;
  fprintf ar,"[%o,%o,%o,%o,%o],\n",[d],[i],Eltseq(a),pass select 1 else 0,known select 1 else 0;
  total+:=1;
 end for;
end for;
for i in [1..#classes[2]] do for j in [i..#classes[2]] do
 a:=classes[2][i]+classes[2][j];pass:=pi(a) eq GG!0;known:=pass and a in qq[4];
 if pass and not known then Append(~passed,<22,i*100000+j,Eltseq(a)>);end if;
 fprintf ar,"[%o,%o,%o,%o,%o],\n",[2,2],[i,j],Eltseq(a),pass select 1 else 0,known select 1 else 0;total+:=1;
end for;end for;
fprintf ar,"]\ncomplete = True\n";delete ar;
fprintf log,"X130_D4_FILTER_DONE tested=%o remaining=%o data=%o complete=true time=%o\n",total,#passed,passed,Cputime(tt);delete log;
print "X130_D4_FILTER_DONE",total,#passed,Cputime(tt);quit;
