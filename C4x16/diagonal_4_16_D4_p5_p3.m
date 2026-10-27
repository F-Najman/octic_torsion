/* Fresh octic M4 two-multiple diagonal. Global C is the subgroup of
   F-rational cusp DIVISORS, not assumed to be the full Jacobian torsion. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_full_aut_hash_m4_16.m");
Attach("common/fast_full_aut_hash_m4_16_relative.m");
Attach("common/class_rational_reduction.m");
function CuspUnitsM4(X)
 FF:=FunctionField(Curve(X));E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
 ai:=aInvariants(E);B:=ai[2];assert ai eq [FF!1,B,B,FF!0,FF!0] and P eq E![0,0];
 assert Discriminant(E) eq B^4*(1-16*B);fs:=[B,1-16*B];
 for k in [1..8] do for j in [0..3] do
  A:=k*Q+j*P;assert A ne E!0 and A[1] ne 0 and A[2] ne 0 and A[1]+B ne 0;
  fs cat:=[A[1],A[2],A[1]+B,A[2]+A[1]+B];
 end for;end for;
 for k in [1..8] do A:=k*Q;for j in [1..3] do R:=A+j*P;assert A ne R and A ne -R;Append(~fs,A[1]-R[1]);end for;end for;
 rawcs:=Cusps(X);groups:=AssociativeArray();
 for c in rawcs do sig:=[Valuation(f,c):f in fs];
  if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 24;
 return [&+[Divisor(c):c in groups[sig]]:sig in sigs],sigs;
end function;
Free:=FreeAbelianGroup(24);Xs:=[* *];Gs:=[* *];hs:=[* *];ms:=[* *];css:=[* *];bases:=[* *];cgs:=[* *];fps:=[];
sigs:=[];ds:=[];
for p in [5,9] do
 X:=MDX11(4,16,GF(p));C:=Curve(X);assert Genus(X) eq 13;cs,ss:=CuspUnitsM4(X);dd:=[Degree(c):c in cs];
 assert {*d:d in dd*} eq {*1^^16,2^^4,4^^4*};
 if p eq 5 then sigs:=ss;ds:=dd;else assert ss eq sigs and dd eq ds;end if;
 c0:=Support(cs[Index(ds,1)])[1];print "M4N16_D4_DIAGONAL_CUSP_SETUP",p,ds,Cputime(t0);
 G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
 J:=JH(64,[17]);assert Dimension(J) eq 13;fp:=FrobeniusPolynomial(J,p eq 9 select 3 else 5);expected:=p eq 9 select Evaluate(fp,1)*Evaluate(fp,-1) else Evaluate(fp,1);assert &*[n:n in Invariants(G)|n ne 0] eq Integers()!expected;
 hg:=[m2(cs[i]-ds[i]*Divisor(c0)):i in [1..24]];h:=hom<Free->G|hg>;H:=Image(h);
 assert Invariants(H) eq [2,2,2,4,4,20,40,40,80];
 if p eq 5 then assert Invariants(quo<G|H>) eq [2,2,2,2,2,0];end if;
 Append(~Xs,X);Append(~Gs,G);Append(~hs,h);Append(~ms,m2);Append(~css,cs);Append(~bases,c0);Append(~cgs,hg);Append(~fps,Coefficients(fp));
 print "M4N16_D4_DIAGONAL_GROUP",p,Invariants(G),Invariants(H),Invariants(quo<G|H>),Cputime(t0);
end for;
h5:=hs[1];h9:=hs[2];G5:=Gs[1];G9:=Gs[2];
for v in Generators(Kernel(h5)) do assert h9(Free!v) eq G9!0;end for;
for v in Generators(Kernel(h9)) do assert h5(Free!v) eq G5!0;end for;
X:=Xs[1];C:=Curve(X);c0:=bases[1];m2:=ms[1];source:=[* *];cache:=AssociativeArray();predictions:=AssociativeArray();
for d in [3,4] do
 xs:=NoncuspidalPlaces(X,d);assert #xs eq (d eq 3 select 64 else 160);
 for i in [1..#xs] do
  a:=m2(Divisor(xs[i])-d*Divisor(c0));inds:=d eq 3 select [j:j in [1..24]|ds[j] eq 1] else [0];
  for j in inds do
   cl:=a+(j eq 0 select G5!0 else cgs[1][j]);tw:=2*cl;key:=Eltseq(tw);
   if not IsDefined(cache,key) then cache[key]:=h9(Free!(tw@@h5));end if;
   pred:=cache[key];Append(~source,[*d,i,j,Eltseq(cl),Eltseq(pred)*]);ky:=Eltseq(pred);
   if IsDefined(predictions,ky) then Append(~predictions[ky],#source);else predictions[ky]:=[#source];end if;
  end for;
 end for;
end for;
assert #source eq 1184;print "M4N16_D4_DIAGONAL_SOURCE",#source,#Keys(predictions),Cputime(t0);
X:=Xs[2];C:=Curve(X);c0:=bases[2];m2:=ms[2];rawcs:=Cusps(X);
assert {*Degree(c):c in rawcs*} eq {*1^^24,2^^8*};
rat:=[c:c in rawcs|Degree(c) eq 1];ratcl:=[m2(Divisor(c)-Divisor(c0)):c in rat];
places:=[NoncuspidalPlaces(X,d):d in [1..4]];assert #places[1] eq 0 and #places[2] eq 0 and #places[3] eq 320;
rs:=[* *];coverage:=[* *];pc:=[* *];
for d in [3..4] do
 rr,labels,keys:=PlacesUpToFullAutM4SixteenRelative(X,places[d]);Append(~rs,rr);
 x:=rr[1];E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);P:=Ls`P;Q:=Ls`Q;orb:={x};
 for par in M4SixteenAutParameters() do a,b,c,dd:=Explode(par);La:=rec<recformat<P:PtEll,Q:PtEll>|P:=a*P+4*b*Q,Q:=c*P+dd*Q>;Include(~orb,ModuliPointFast(X,E,La));end for;
 assert orb eq {places[d][i]:i in [1..#places[d]]|labels[i] eq 1};
 Append(~coverage,[*d,keys,labels,[Index(places[d],x):x in rr]*]);
 Append(~pc,[Eltseq(m2(Divisor(x)-d*Divisor(c0))):x in rr]);
 print "M4N16_D4_DIAGONAL_F9_TARGET_PLACES",d,#places[d],#rr,"reference_orbit",#orb,Cputime(t0);
end for;
targets:=[* *];matches:=[* *];counts:=[0,0];
for i in [1..#rs[1]] do for j in [1..#rat] do
 cl:=G9!pc[1][i]+ratcl[j];Append(~targets,[*3,i,j,Eltseq(cl)*]);ky:=Eltseq(2*cl);
 if IsDefined(predictions,ky) then Append(~matches,[*#targets,predictions[ky]*]);end if;counts[1]+:=1;
end for;end for;
for i in [1..#rs[2]] do
 cl:=G9!pc[2][i];Append(~targets,[*4,i,0,Eltseq(cl)*]);ky:=Eltseq(2*cl);
 if IsDefined(predictions,ky) then Append(~matches,[*#targets,predictions[ky]*]);end if;counts[2]+:=1;
end for;
Write("C4x16/data/diagonal_4_16_D4_p5_p3_arithmetic.txt",Sprint([*sigs,ds,Invariants(G5),[Eltseq(r):r in Rows(RelationMatrix(G5))],Invariants(G9),[Eltseq(r):r in Rows(RelationMatrix(G9))],[[Eltseq(v):v in cg]:cg in cgs],source,targets,matches,[#s:s in places],[#s:s in rs],coverage,[Eltseq(a):a in ratcl],pc,counts,fps*]):Overwrite:=true);
print "M4N16_D4_DIAGONAL_F9_DONE",#source,#targets,counts,#matches,Cputime(t0);quit;
