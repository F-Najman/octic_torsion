/* Fresh relative degree4 diagonal search for6x12. All468 quartic
   source places at7 and all noncuspidal degree4 patterns at13 are included.
   Cubic+cusp and entirelycuspidal p7 branches need the separate arguments. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");
Attach("common/fast_hecke_11.m");
Attach("C6x12/model_6_12_root.m");
Attach("common/class_rational_reduction.m");
Gs:=[* *];Hs:=[* *];hs:=[* *];Xs:=[* *];ms:=[* *];css:=[* *];c0s:=[* *];sigs0:=[];ds0:=[];
for p in [7,13] do
 X:=MDX612(GF(p));C:=Curve(X);FF:=FunctionField(C);
 E:=EllipticCurve([FF!a:a in X`_E]);Q:=E![0,0];assert 12*Q eq E!0;
 px:=FF!X`_P[1];ai:=aInvariants(E);R<Y>:=PolynomialRing(FF);
 ys:=Roots(Y^2+(ai[1]*px+ai[3])*Y-px^3-ai[2]*px^2-ai[4]*px-ai[5]);
 assert #ys eq 2;P:=E![px,ys[1][1]];
 z:=PrimitiveElement(GF(p))^((p-1) div 6);
 if WeilPairing(P,2*Q,6) ne z then P:=-P;end if;
 assert WeilPairing(P,2*Q,6) eq z;
 assert ai[2] eq ai[3];fs:=[-ai[2],1-ai[1],Discriminant(E)];
 anchor:=P[1];
 for i in [0..5] do
  for j in [0..11] do
   A:=i*P+j*Q;if A eq E!0 then continue;end if;
   if A ne Q and A ne -Q then assert A[1] ne 0;Append(~fs,A[1]);end if;
   if A ne P and A ne -P then assert A[1] ne anchor;Append(~fs,A[1]-anchor);end if;
  end for;
 end for;
 assert &and[f ne 0:f in fs];cs:=Cusps(X);groups:=AssociativeArray();
 for c in cs do
  s:=[Valuation(f,c):f in fs];
  if IsDefined(groups,s) then Append(~groups[s],c);else groups[s]:=[c];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));blocks:=[&+[Divisor(c):c in groups[s]]:s in sigs];
 ds:=[Degree(D):D in blocks];
 print "M6N12_D4_UNIT_BLOCKS",p,#sigs,ds,Cputime(t0);
 if p eq 7 then sigs0:=sigs;ds0:=ds;Free:=FreeAbelianGroup(#sigs);
 else assert sigs eq sigs0 and ds eq ds0;end if;
 c0:=Support(blocks[Index(ds,1)])[1];
 G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
 gs:=[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..#blocks]];
 h:=hom<Free->G|gs>;H:=Image(h);QH:=quo<G|H>;
 print "M6N12_D4_CUSP_TRANSPORT_GROUP",p,Invariants(G),Invariants(H),Invariants(QH),Cputime(t0);
 Append(~Gs,G);Append(~Hs,H);Append(~hs,h);Append(~Xs,X);Append(~ms,m2);Append(~css,cs);Append(~c0s,c0);
end for;
for a in Generators(Kernel(hs[1])) do assert hs[2](Free!a) eq Gs[2]!0;end for;
for a in Generators(Kernel(hs[2])) do assert hs[1](Free!a) eq Gs[1]!0;end for;
G7:=Gs[1];G13:=Gs[2];H7:=Hs[1];h7:=hs[1];h13:=hs[2];
QH,qh:=quo<G7|H7>;assert Invariants(QH) eq [4,4,0];
predictions:=AssociativeArray();source_records:=[* *];
X:=Xs[1];m2:=ms[1];c0:=c0s[1];ps7:=NoncuspidalPlaces(X,4);assert #ps7 eq 468;
for i in [1..#ps7] do
 cl:=m2(Divisor(ps7[i])-4*Divisor(c0));r:=Order(qh(cl));assert r in [2,4];
 a:=h13(Free!((r*cl)@@h7));key:=<r,Eltseq(a)>;
 if IsDefined(predictions,key) then Append(~predictions[key],i);else predictions[key]:=[i];end if;
 Append(~source_records,[*i,r,Eltseq(cl),Eltseq(a)*]);
end for;
print "M6N12_D4_DIAGONAL_PREDICTIONS",#source_records,#Keys(predictions),Cputime(t0);
X:=Xs[2];C:=Curve(X);m2:=ms[2];c0:=c0s[2];cs:=css[2];
ps:=[NoncuspidalPlaces(X,d):d in [1..4]];assert #ps[1] eq 0;
print "M6N12_D4_DIAGONAL_COUNTS13",[#s:s in ps],Cputime(t0);
pcs:=[* [G13|] *];
for d in [2..4] do
 temp:=[];
 for i in [1..#ps[d]] do
  Append(~temp,m2(Divisor(ps[d][i])-d*Divisor(c0)));
  if i mod 500 eq 0 then print "M6N12_D4_DIAGONAL_CLASS_PROGRESS",d,i,#ps[d],Cputime(t0);end if;
 end for;
 Append(~pcs,temp);
end for;
rat:=[c:c in cs|Degree(c) eq 1];quad:=[c:c in cs|Degree(c) eq 2];
rc:=[m2(Divisor(c)-Divisor(c0)):c in rat];
qc:=[m2(Divisor(c)-2*Divisor(c0)):c in quad];
complements:=[* [G13!0],rc,qc cat [rc[i]+rc[j]:j in [i..#rc],i in [1..#rc]] *];
matches:=[* *];records:=[* *];counts:=[0,0,0,0];
procedure Check(cl,descr,~matches,~records)
 hits:=[* *];
 for r in [2,4] do
  key:=<r,Eltseq(r*cl)>;
  if IsDefined(predictions,key) then
   for i in predictions[key] do Append(~hits,i);end for;
  end if;
 end for;
 Append(~records,[*descr,Eltseq(cl),hits*]);
 if #hits gt 0 then Append(~matches,[*descr,hits,Eltseq(cl)*]);end if;
end procedure;
for d in [2..4] do
 for i in [1..#pcs[d]] do
  for j in [1..#complements[5-d]] do
   cl:=pcs[d][i]+complements[5-d][j];
   Check(cl,[*d,i,j*],~matches,~records);counts[d-1]+:=1;
  end for;
 end for;
 print "M6N12_D4_DIAGONAL_TYPE",d,counts[d-1],#matches,Cputime(t0);
end for;
for i in [1..#pcs[2]] do
 for j in [i..#pcs[2]] do
  cl:=pcs[2][i]+pcs[2][j];Check(cl,[*22,i,j*],~matches,~records);counts[4]+:=1;
 end for;
end for;
Write("C6x12/data/diagonal_6_12_D4_root_arithmetic.txt",Sprint([*sigs0,ds0,
 [Invariants(G):G in Gs],[[Eltseq(v):v in Rows(RelationMatrix(G))]:G in Gs],
 [[Eltseq(h(Free.i)):i in [1..Ngens(Free)]]:h in hs],source_records,records,matches,counts*]):Overwrite:=true);
printf "M6N12_D4_DIAGONAL_DONE source_places=%o patterns13=%o matches=%o complete=true time=%o\n",#source_records,counts,#matches,Cputime(t0);quit;
