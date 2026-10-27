/* Fresh bounded M4 genus13 global cusp-image and Jacobian-order probe.
   At5 Frobenius generates Gal(Q(zeta16)/Q(i)), so ALL24 closed cusps
   are degree-preserving lifts of global F-rational cusp orbits. */
SetColumns(0);SetQuitOnError(true);SetVerbose("ClassGroup",1);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_full_aut_hash_m4_16.m");
J:=JH(64,[17]);assert Dimension(J) eq 13;fps:=[];orders:=[];
for p in [5,13,17,29] do
 fp:=FrobeniusPolynomial(J,p);assert Degree(fp) eq 26;Append(~fps,Coefficients(fp));Append(~orders,Integers()!Evaluate(fp,1));
end for;print "M4N16_D4_SPLIT_FROB_ORDERS",orders,GCD(orders),Cputime(t0);
X:=MDX11(4,16,GF(5));C:=Curve(X);assert Genus(X) eq 13;
cs:=Cusps(X);assert {*Degree(c):c in cs*} eq {*1^^16,2^^4,4^^4*};
FF:=FunctionField(C);EE:=EllipticCurve([FF!f:f in X`_E]);PP:=EE![FF!f:f in X`_P];QQ:=EE![FF!f:f in X`_Q];
ai:=aInvariants(EE);B:=ai[2];assert ai eq [FF!1,B,B,FF!0,FF!0] and PP eq EE![0,0];
assert Discriminant(EE) eq B^4*(1-16*B);fs:=[B,1-16*B];
for k in [1..8] do
 for j in [0..3] do
  A:=k*QQ+j*PP;assert A ne EE!0 and A[1] ne 0 and A[2] ne 0 and A[1]+B ne 0;
  fs cat:=[A[1],A[2],A[1]+B,A[2]+A[1]+B];
 end for;
end for;
for k in [1..8] do
 A:=k*QQ;for j in [1..3] do R:=A+j*PP;assert A ne R and A ne -R;Append(~fs,A[1]-R[1]);end for;
end for;
sigs:=[[Valuation(f,c):f in fs]:c in cs];print "M4N16_D4_UNIT_SIGNATURE_COUNT",#Seqset(sigs),#cs,Cputime(t0);
ord:=Sort([<sigs[i],i>:i in [1..#cs]]);cs:=[cs[x[2]]:x in ord];sigs:=[x[1]:x in ord];
ds:=[Degree(c):c in cs];c0:=cs[Index(ds,1)];
print "M4N16_D4_GLOBAL_CUSP_SIGNATURES",#sigs,ds,Cputime(t0);
AF:=AlgorithmicFunctionField(FF);G,_,phi:=ClassGroup(AF:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
assert &*[n:n in Invariants(G)|n ne 0] eq orders[1];
m2:=map<DivisorGroup(C)->G|D:->phi(FunctionFieldDivisor(D))>;
hg:=[m2(Divisor(cs[i])-ds[i]*Divisor(c0)):i in [1..24]];Free:=FreeAbelianGroup(24);h:=hom<Free->G|hg>;H:=Image(h);Q,qp:=quo<G|H>;
ho:=#H;assert GCD(orders) mod ho eq 0;
print "M4N16_D4_GLOBAL_CUSP_GROUP",Invariants(G),Invariants(H),Invariants(Q),"global_order_ratio_bound",GCD(orders) div ho,Cputime(t0);
records:=[* *];contained:=0;targets:={G|};
for dd in [3,4] do
 xs:=NoncuspidalPlaces(X,dd);rs,labels,keys:=PlacesUpToFullAutM4Sixteen(X,xs);
 assert #xs eq (dd eq 3 select 64 else 160) and #rs eq (dd eq 3 select 1 else 4);
 for i in [1..#rs] do
  a:=m2(Divisor(rs[i])-dd*Divisor(c0));
  inds:=dd eq 3 select [j:j in [1..24]|ds[j] eq 1] else [0];
  for j in inds do
   cl:=a+(j eq 0 select G!0 else hg[j]);inside:=cl in H;
   if inside then contained+:=1;Include(~targets,cl);end if;
   Append(~records,[*dd,i,j,Eltseq(cl),inside*]);
   print "M4N16_D4_GLOBAL_CUSP_TEST",dd,i,j,inside,Eltseq(cl),Cputime(t0);
  end for;
 end for;
end for;assert #records eq 20;
short:=[* *];KB:=Matrix(Integers(),[Eltseq(Free!v):v in Generators(Kernel(h))]);LL:=LatticeWithBasis(LLL(KB));
for tv in Sort([Eltseq(t):t in targets]) do
 cl:=G!tv;v:=Eltseq(cl@@h);w:=ClosestVector(LL,Vector(Rationals(),v));sv:=[v[j]-Integers()!w[j]:j in [1..24]];
 assert h(Free!sv) eq cl;Append(~short,[*tv,sv*]);
end for;
Write("C4x16/data/class_4_16_D4_p5_v3_arithmetic.txt",Sprint([*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(g):g in hg],sigs,ds,Index(ds,1),orders,fps,records,short*]):Overwrite:=true);
print "M4N16_D4_GLOBAL_CUSP_PROBE_DONE",#records,contained,#targets,Cputime(t0);quit;
