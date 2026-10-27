/* Necessary three-multiple diagonal test. ALL F4 effective degree4 divisors
   are included. AtF7 the entirely cuspidal branch has a global F-cusp lift
   and is excluded separately using odd-prime injectivity and gon_F>=6.
   Noncuspidal F7 patterns are single2/3/4 with every cusp complement,
   and2+2 including repetition. No global full-generation claim is used. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");
Free:=FreeAbelianGroup(20);Xs:=[* *];Gs:=[* *];Hs:=[* *];hs:=[* *];ms:=[* *];
css:=[* *];basecs:=[* *];signatures:=[];degrees:=[];
function WeightedVectors(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else [];end if;
 result:=[];
 for k in [0..d div ds[1]] do
  for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do Append(~result,[k] cat v);end for;
 end for;
 return result;
end function;
for p in [4,7] do
 X:=MDX11(3,15,GF(p));C:=Curve(X);FF:=FunctionField(C);
 E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
 ai:=aInvariants(E);assert ai[2] eq 0 and ai[4] eq 0 and ai[5] eq 0;
 assert P eq E![0,0];fs:=[ai[3],ai[1]^3-27*ai[3]];
 for k in [1..7] do
  for j in [0..2] do
   A:=k*Q+j*P;assert A[3] ne 0 and A[1] ne 0 and A[2] ne 0;
   Append(~fs,A[1]);Append(~fs,A[2]);
  end for;
  A:=k*Q;B:=A+P;assert A ne B and A ne -B;Append(~fs,A[1]-B[1]);
 end for;
 groups:=AssociativeArray();cs:=Cusps(X);
 for cusp in cs do
  sig:=[Valuation(f,cusp):f in fs];
  if IsDefined(groups,sig) then Append(~groups[sig],cusp);else groups[sig]:=[cusp];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 20;
 blocks:=[&+[Divisor(c):c in groups[sig]]:sig in sigs];
 ds:=[Degree(o):o in blocks];assert {*d:d in ds*} eq {*1^^16,4^^4*};
 if p eq 4 then signatures:=sigs;degrees:=ds;
 else assert sigs eq signatures and ds eq degrees;end if;
 c0:=Support(blocks[Index(ds,1)])[1];
 print "M3N15_DIAGONAL_CUSP_BLOCKS",p,#sigs,ds,Cputime(t0);
 G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
 gs:=[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..20]];
 h:=hom<Free->G|gs>;H:=Image(h);
 assert Invariants(H) eq [3,3,24,96,96];
 J:=JH(45,[16]);assert Dimension(J) eq 9;
 fp:=FrobeniusPolynomial(J,p eq 4 select 2 else 7);
 expected:=(p eq 4) select Evaluate(fp,1)*Evaluate(fp,-1) else Evaluate(fp,1);
 assert &*[n:n in Invariants(G)|n ne 0] eq Integers()!expected;
 Append(~Xs,X);Append(~Gs,G);Append(~Hs,H);Append(~hs,h);Append(~ms,m2);
 Append(~css,cs);Append(~basecs,c0);
 print "M3N15_DIAGONAL_GROUP",p,Invariants(G),Invariants(H),Cputime(t0);
end for;
G4:=Gs[1];G7:=Gs[2];H4:=Hs[1];H7:=Hs[2];h4:=hs[1];h7:=hs[2];
for v in Generators(Kernel(h4)) do assert h7(Free!v) eq G7!0;end for;
for v in Generators(Kernel(h7)) do assert h4(Free!v) eq G4!0;end for;
images:=[h7(Free!((G4!H4.i)@@h4)):i in [1..Ngens(H4)]];
transport:=hom<H4->G7|images>;
for i in [1..20] do assert transport(H4!h4(Free.i)) eq h7(Free.i);end for;
print "M3N15_DIAGONAL_CUSP_KERNELS_MATCH",#H4,#H7,Cputime(t0);
S:=AssociativeArray();states4:=[* *];X:=Xs[1];m2:=ms[1];c0:=basecs[1];cs:=css[1];
sc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];
vs:=WeightedVectors([Degree(c):c in cs],4);assert #vs eq 5000;
for i in [1..#vs] do
 v:=vs[i];cl:=&+[v[j]*sc[j]:j in [1..#cs]];a:=H4!(3*cl);pred:=transport(a);
 Append(~states4,[*0,i,Eltseq(cl),Eltseq(a),Eltseq(pred)*]);key:=Eltseq(pred);
 if IsDefined(S,key) then Append(~S[key],#states4);else S[key]:=[#states4];end if;
end for;
xs:=NoncuspidalPlaces(X,4);assert #xs eq 102;
assert [#NoncuspidalPlaces(X,d):d in [1..3]] eq [0,0,0];
for i in [1..#xs] do
 cl:=m2(Divisor(xs[i])-4*Divisor(c0));a:=H4!(3*cl);pred:=transport(a);
 Append(~states4,[*4,i,Eltseq(cl),Eltseq(a),Eltseq(pred)*]);key:=Eltseq(pred);
 if IsDefined(S,key) then Append(~S[key],#states4);else S[key]:=[#states4];end if;
end for;
assert #states4 eq 5102;
print "M3N15_DIAGONAL_STATES4",#states4,#Keys(S),Cputime(t0);
X:=Xs[2];C:=Curve(X);m2:=ms[2];c0:=basecs[2];cs:=css[2];
rat:=[c:c in cs|Degree(c) eq 1];assert #rat eq 16;
ratcl:=[m2(Divisor(c)-Divisor(c0)):c in rat];
ps:=[NoncuspidalPlaces(X,d):d in [1..4]];assert #ps[1] eq 0;
print "M3N15_DIAGONAL_COUNTS7",[#ss:ss in ps],Cputime(t0);
pcs:=[* [G7|] *];
for d in [2..4] do Append(~pcs,[m2(Divisor(x)-d*Divisor(c0)):x in ps[d]]);end for;
states7:=[* *];pairs:=[* *];counts:=[];
for d in [2..4] do
 vs:=WeightedVectors([1:i in [1..16]],4-d);subcount:=0;
 for i in [1..#ps[d]] do
  for j in [1..#vs] do
   v:=vs[j];cl:=pcs[d][i]+&+[v[k]*ratcl[k]:k in [1..16]];
   Append(~states7,[*d,i,j,Eltseq(cl)*]);key:=Eltseq(3*cl);subcount+:=1;
   if IsDefined(S,key) then
    for i4 in S[key] do Append(~pairs,[*i4,#states7*]);end for;
   end if;
  end for;
 end for;
 Append(~counts,subcount);
end for;
subcount:=0;
for i in [1..#ps[2]] do
 for j in [i..#ps[2]] do
  cl:=pcs[2][i]+pcs[2][j];Append(~states7,[*22,i,j,Eltseq(cl)*]);key:=Eltseq(3*cl);subcount+:=1;
  if IsDefined(S,key) then for i4 in S[key] do Append(~pairs,[*i4,#states7*]);end for;end if;
 end for;
end for;
Append(~counts,subcount);
Write("C3x15/data/diagonal_3_15_D4_arithmetic.txt",Sprint([*
 signatures,degrees,Invariants(G4),[Eltseq(v):v in Rows(RelationMatrix(G4))],
 Invariants(G7),[Eltseq(v):v in Rows(RelationMatrix(G7))],
 [Eltseq(h4(Free.i)):i in [1..20]],[Eltseq(h7(Free.i)):i in [1..20]],
 Invariants(H4),[Eltseq(v):v in Rows(RelationMatrix(H4))],
 [Eltseq(G4!H4.i):i in [1..Ngens(H4)]],[Eltseq(a):a in images],
 states4,states7,pairs,counts *]):Overwrite:=true);
printf "M3N15_DIAGONAL_DONE states4=%o states7=%o patterncounts=%o matches=%o complete=true time=%o\n",#states4,#states7,counts,#pairs,Cputime(t0);quit;
