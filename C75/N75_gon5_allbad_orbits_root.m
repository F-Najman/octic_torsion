/* New canonical degree20 quotient model attempt; no exclusion claim. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
R<T,x,y>:=PolynomialRing(GF(2),3);B:=FieldOfFractions(R);
F25:=eval Read("mdmagma/models_X1_n/FFFc25.txt");
xx:=B!x;yy:=B!y;r:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));s:=(xx*yy-yy+1)/(xx*yy);
b:=r*s*(r-1);c:=s*(r-1);
// The exact degree5 modular unit from the archived vanHoeij gonality table.
f3:=b;f4:=c;f5:=b-c;f6:=s-1;f7:=s-r;f8:=r*s-2*r+1;f9:=s^2-s-r+1;
f10:=xx-yy+1;f12:=xx-yy;f13:=xx^3*yy-xx^2*yy^2-xx^2*yy+xx*yy^2-yy+1;
h:=f3*f4*f6^3*f7^2*f10*f12^3/(f5*f8^2*f9^2*f13^3);
N:=R!Numerator(h);D:=R!Denominator(h);
print "N75_GON5_CANONICAL_RESULTANT_START",TotalDegree(N),TotalDegree(D),Cputime(tt);
res:=Resultant(F25,N-T*D,3);assert res ne 0;
K<t>:=RationalFunctionField(GF(2));S<X>:=PolynomialRing(K);proj:=hom<R->S|[S!t,X,S!0]>;
pol:=proj(res);factors:=Factorization(pol);
print "N75_GON5_CANONICAL_RESULTANT_FACTORS",[<Degree(z[1]),z[2]>:z in factors],Cputime(tt);
fs:=[z[1]:z in factors|Degree(z[1]) eq 5];assert #fs eq 1;
L<a>:=FunctionField(fs[1]);SY<Y>:=PolynomialRing(L);ev:=hom<R->SY|[SY!t,SY!a,Y]>;
gg:=GCD(ev(F25),ev(N-T*D));assert Degree(gg) eq 1;
ya:=-Coefficient(gg,0)/Coefficient(gg,1);assert Evaluate(ev(F25),ya) eq 0;
xa:=a;ra:=(xa^2*ya-xa*ya+ya-1)/(xa*(xa*ya-1));sa:=(xa*ya-ya+1)/(xa*ya);
ba:=ra*sa*(ra-1);ca:=sa*(ra-1);E:=EllipticCurve([1-ca,-ba,-ba,L!0,L!0]);P:=E![0,0];
assert 25*P eq E!0 and 5*P ne E!0 and Genus(L) eq 12;
print "N75_GON5_CANONICAL_BASE_MODEL_VERIFIED",Degree(L),Genus(L),Cputime(tt);
S3<W3>:=PolynomialRing(L);j25:=jInvariant(E);f:=W3^4+W3^3+1/j25;assert Degree(f) eq 4;
assert IsIrreducible(f);M<w>:=FunctionField(f);
print "N75_GON5_CANONICAL_DEGREE20_GENUS_START",Cputime(tt);
assert Genus(M) eq 73;
print "N75_GON5_CANONICAL_GENUS73_MODEL_VERIFIED",Cputime(tt);

function Gon5Value(xx,yy)
 r:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));s:=(xx*yy-yy+1)/(xx*yy);
 b:=r*s*(r-1);c:=s*(r-1);
 return b*c*(s-1)^3*(s-r)^2*(xx-yy+1)*(xx-yy)^3/((b-c)*(r*s-2*r+1)^2*(s^2-s-r+1)^2*(xx^3*yy-xx^2*yy^2-xx^2*yy+xx*yy^2-yy+1)^3);
end function;
AttachSpec("mdmagma/v2/mdmagma.spec");
xy2:=TateNormalForm_xy(E,2*P);h2:=Gon5Value(xy2[1],xy2[2]);
sig0:=hom<K->L|h2>;sig:=hom<L->L|sig0,xy2[1]>;
assert sig(ya) eq xy2[2] and sig(j25) eq j25;
assert &and[(sig^10)(x) eq x:x in [L!t,a]];
dg:=hom<M->M|map<L->M|x:->M!sig(x)>,w>;
assert &and[(dg^10)(x) eq x:x in [M!t,M!a,w]];
print "N75_D2_AUTOMORPHISM_VERIFIED",Cputime(tt);
bcs:=Poles(j25);assert Sort([Degree(p):p in bcs]) eq [1:i in [1..10]] cat [4,4,10];
selected:=[[i:i in [1..#bcs]|Degree(bcs[i]) eq d][1]:d in [1,4]];
sourcecs:=[];
for i in selected do
 p:=bcs[i];ff,gg:=TwoGenerators(p);u:=UniformizingElement(p);
 raw:=Decomposition(M,p);above:=[q:q in raw|Valuation(M!ff,q) gt 0 and Valuation(M!gg,q) gt 0];
 assert #above eq 2 and Sort([Valuation(M!u,q):q in above]) eq [1,3];
 assert &and[Degree(q) eq Degree(p):q in above];sourcecs cat:=above;
 print "N75_ORBIT_REPRESENTATIVE_CUSPS",i,Degree(p),[<-Valuation(M!j25,q),Degree(q)>:q in above],Cputime(tt);
end for;
assert #sourcecs eq 4;
Om,om:=SpaceOfHolomorphicDifferentials(M);bs:=[om(v):v in Basis(Om)];
C:=Matrix(GF(2),[Eltseq(Cartier(omega)@@om):omega in bs]);
I:=IdentityMatrix(GF(2),73);U:=RowSpace(C-I);BM:=BasisMatrix(U);
eta:=Differential(M!t);deta:=Differential(dg(M!t));
Md:=Matrix(GF(2),[Eltseq((dg(omega/eta)*deta)@@om):omega in bs]);
assert Md^10 eq I and Md*C eq C*Md and RowSpace(BM*Md) eq U;
print "N75_CARTIER_DIAMOND_READY",Rank(C),Dimension(U),Cputime(tt);
Write("C75/data/N75_allbad_orbit_matrices_root.m",Sprintf("CartierMatrix:=Matrix(GF(2),73,73,%o);Diamond2Matrix:=Matrix(GF(2),73,73,%o);ImageBasis:=Matrix(GF(2),%o,73,%o);",Eltseq(C),Eltseq(Md),Dimension(U),Eltseq(BM)):Overwrite:=true);
ds:=[];widths:=[];steps:=[];orbs:=[];sourceOf:=[];
for rep in [1..4] do
 q:=sourcecs[rep];d:=Degree(q);ww:=-Valuation(M!j25,q);len:=d eq 1 select 10 else 2;
 assert ww in (d eq 1 select [25,75] else [5,15]);
 orb:=[];
 for pow in [0..len-1] do
  Append(~ds,d);Append(~widths,ww);Append(~steps,d eq 1 and ww eq 25 select 2 else 1);
  Append(~orb,#ds);Append(~sourceOf,rep);
 end for;Append(~orbs,orb);
 // Diamond2 has the stated closed-cusp orbit length by the independent
 // primitive Tate-label permutation. Check its stabilizing kernel later.
end for;
assert #ds eq 24 and Sort(ds) eq [1:i in [1..20]] cat [4,4,4,4];
assert Sort(widths) eq Sort([25:i in [1..10]] cat [75:i in [1..10]] cat [5,5,15,15]);
allowed:=[1..24];perm:=[];
for orb in orbs do for j in [1..#orb] do Append(~perm,orb[(j mod #orb)+1]);end for;end for;
bindex:=sourceOf;
print "N75_ABSTRACT_CUSP_ORBITS",orbs,ds,widths,steps,Cputime(tt);
Write("C75/data/N75_allbad_orbit_kernels_root.m",Sprintf("N:=75;D:=8;gx:=73;udim:=%o;ds:=%o;CuspWidths:=%o;CuspBaseIndices:=%o;CuspPermutation2:=%o;AllowedCusps:=%o;MultiplicitySteps:=%o;kernelspaces:=AssociativeArray();basisjets:=AssociativeArray();localjets:=AssociativeArray();sourcekernels:=AssociativeArray();",Dimension(U),ds,widths,bindex,perm,allowed,steps):Overwrite:=true);
VO:=VectorSpace(GF(2),73);allJets:=AssociativeArray();
for oi in [1..#orbs] do
 orb:=orbs[oi];rep:=orb[1];qsource:=sourcecs[oi];step:=steps[rep];cols:=[];annprev:=sub<VO|>;
 for m in [step..(8 div ds[rep]) by step] do
  print "N75_RR_SOURCE_START",rep,ds[rep],m,Cputime(tt);
  VV,vm:=DifferentialSpace(m*Divisor(qsource));
  W:=sub<VO|[VO!Eltseq(vm(v)@@om):v in Basis(VV)]>;
  ann:=Nullspace(Transpose(BasisMatrix(W)));assert annprev subset ann;
  for v in Basis(ann) do
   if v notin annprev then Append(~cols,VO!v);annprev:=sub<VO|[VO!x:x in Basis(annprev)] cat [VO!v]>;end if;
  end for;
  while #cols lt m*ds[rep] do Append(~cols,VO!0);end while;
  A:=Transpose(Matrix(GF(2),[Eltseq(v):v in cols]));assert Nullspace(A) eq W and Ncols(A) eq m*ds[rep];
  Write("C75/data/N75_allbad_orbit_kernels_root.m",Sprintf("sourcekernels[<%o,%o>]:=Matrix(GF(2),%o,73,%o);",rep,m,Dimension(W),Eltseq(BasisMatrix(W))));
  for target in orb do
   powers:=[];k:=rep;
   for pow in [0..9] do if k eq target then Append(~powers,pow);end if;k:=perm[k];end for;
   assert #powers gt 0;D:=Md^powers[1];WT:=RowSpace(BasisMatrix(W)*D^-1);AT:=D*A;
   assert Nullspace(AT) eq WT;
   assert &and[RowSpace(BasisMatrix(W)*(Md^pow)^-1) eq WT:pow in powers];
   JT:=BM*AT;allJets[<target,m>]:=JT;
   Write("C75/data/N75_allbad_orbit_kernels_root.m",Sprintf("kernelspaces[<%o,%o>]:=Matrix(GF(2),%o,73,%o);basisjets[<%o,%o>]:=Matrix(GF(2),73,%o,%o);localjets[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);",target,m,Dimension(WT),Eltseq(BasisMatrix(WT)),target,m,Ncols(AT),Eltseq(AT),target,m,Nrows(JT),Ncols(JT),Eltseq(JT)));
  end for;
  print "N75_RR_SOURCE_DONE",rep,m,Dimension(W),Rank(BM*A),Cputime(tt);
 end for;
end for;
function AdmissibleVectors(d,i,ds,steps)
 if i gt #ds then return d eq 0 select [[]] else [];end if;
 rows:=[];
 for m in [0..(d div ds[i]) by steps[i]] do
  if ds[i] gt 8 and m gt 0 then continue;end if;
  for v in AdmissibleVectors(d-m*ds[i],i+1,ds,steps) do Append(~rows,[m] cat v);end for;
 end for;return rows;
end function;
vecs:=AdmissibleVectors(8,1,ds,steps);assert #vecs eq 131790;
bad:=[];tested:=0;
for v in vecs do
 A:=ZeroMatrix(GF(2),Dimension(U),0);
 for i in allowed do if v[i] gt 0 then A:=HorizontalJoin(A,allJets[<i,v[i]>]);end if;end for;
 assert Ncols(A) eq 8;tested+:=1;
 if Rank(A) ne 8 then Append(~bad,v);print "N75_FI_EXCEPTION",tested,Rank(A),v;end if;
end for;
print "N75_ALLBAD_FI_DONE",tested,#bad,"complete",tested eq 131790,Cputime(tt);
Write("C75/data/N75_allbad_orbit_result_root.txt",Sprint([*tested,bad,ds,widths,steps,orbs*]):Overwrite:=true);
quit;
