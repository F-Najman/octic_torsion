/* Exact refinement and short global F-cusp lifts. No new class-group
   computation and no full-generation assumption. */
SetColumns(0);SetQuitOnError(true);t0:=Cputime();
a:=eval Read("C3x15/data/diagonal_3_15_D4_arithmetic.txt");
G4:=AbelianGroup(a[3]);G7:=AbelianGroup(a[5]);Free:=FreeAbelianGroup(20);
assert [Eltseq(r):r in Rows(RelationMatrix(G4))] eq a[4];
assert [Eltseq(r):r in Rows(RelationMatrix(G7))] eq a[6];
h4:=hom<Free->G4|[G4!v:v in a[7]]>;h7:=hom<Free->G7|[G7!v:v in a[8]]>;
H4:=Image(h4);H7:=Image(h7);assert Invariants(H4) eq [3,3,24,96,96];
assert Invariants(H7) eq Invariants(H4);
for v in Generators(Kernel(h4)) do assert h7(Free!v) eq G7!0;end for;
for v in Generators(Kernel(h7)) do assert h4(Free!v) eq G4!0;end for;
exactpairs:=[];targets:={G7|};
for pair in a[15] do
 i,j:=Explode(pair);cl4:=G4!a[13][i][3];cl7:=G7!a[14][j][4];
 assert cl4 in H4 and cl7 in H7;
 lift7:=h7(Free!(cl4@@h4));
 if lift7 eq cl7 then Append(~exactpairs,pair);Include(~targets,cl7);end if;
end for;
assert #exactpairs eq 12 and #targets eq 12;
KB:=Matrix(Integers(),[Eltseq(Free!v):v in Generators(Kernel(h7))]);
LL:=LatticeWithBasis(LLL(KB));vectors:=[];targs:=Sort([Eltseq(t):t in targets]);
for tv in targs do
 cl:=G7!tv;v:=Eltseq(cl@@h7);w:=ClosestVector(LL,Vector(Rationals(),v));
 short:=[v[i]-Integers()!w[i]:i in [1..20]];assert h7(Free!short) eq cl;
 Append(~vectors,short);print "M3N15_SHORT_CUSP_LIFT",short,tv;
end for;
Write("C3x15/data/M3N15_D4_exact_cusp_lifts.txt",
 Sprint([*3,15,4,a[1],a[2],Index(a[2],1),vectors,targs,exactpairs*]):Overwrite:=true);
print "M3N15_EXACT_CUSP_LIFTS_DONE",#vectors,Cputime(t0);quit;
