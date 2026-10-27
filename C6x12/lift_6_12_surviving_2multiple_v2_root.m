/* Fresh relative degree4 diagonal search for6x12. All468 quartic
   source places at7 and all noncuspidal degree4 patterns at13 are included.
   Cubic+cusp and entirelycuspidal p7 branches need the separate arguments. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();SetLogFile("C6x12/logs/lift_6_12_surviving_2multiple_v2_root.log");
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");
Attach("common/fast_hecke_11.m");
Attach("C6x12/model_6_12_root.m");
Attach("common/class_rational_reduction.m");
Gs:=[* *];Hs:=[* *];hs:=[* *];Xs:=[* *];ms:=[* *];css:=[* *];c0s:=[* *];sigs0:=[];ds0:=[];
for p in [7] do
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

G:=Gs[1];h:=hs[1];X:=Xs[1];m2:=ms[1];c0:=c0s[1];xs:=NoncuspidalPlaces(X,4);
indices:=[127,128,129,130,131,132,343,344,345,346,347,348];
rows:=[* *];
for i in indices do
 cl:=m2(Divisor(xs[i])-4*Divisor(c0));assert 2*cl in Image(h);
 a:=Free!((2*cl)@@h);assert h(a) eq 2*cl;
 Append(~rows,[*i,Eltseq(a)*]);
 print "M6N12_GLOBAL_L8_CLASS",i,Eltseq(a),Cputime(t0);
end for;
Write("C6x12/data/M6N12_global_L8_cusp_classes_root.m",Sprint([*sigs0,ds0,rows*]):Overwrite:=true);
print "M6N12_GLOBAL_L8_CLASSES_COMPLETE",Cputime(t0);quit;
