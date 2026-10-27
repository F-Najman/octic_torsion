SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/class_rational_reduction.m");
Attach("common/fast_full_aut_hash_m3_15.m");
lift:=eval Read("C3x15/data/M3N15_D4_exact_cusp_lifts.txt");
assert lift[1..3] eq [*3,15,4*];
X:=MDX11(3,15,GF(7));C:=Curve(X);FF:=FunctionField(C);
E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
ai:=aInvariants(E);assert ai[2] eq 0 and ai[4] eq 0 and ai[5] eq 0;
fs:=[ai[3],ai[1]^3-27*ai[3]];
for k in [1..7] do
 for j in [0..2] do
  A:=k*Q+j*P;assert A[3] ne 0 and A[1] ne 0 and A[2] ne 0;
  Append(~fs,A[1]);Append(~fs,A[2]);
 end for;
 A:=k*Q;B:=A+P;assert A ne B and A ne -B;Append(~fs,A[1]-B[1]);
end for;
groups:=AssociativeArray();cs:=Cusps(X);
assert {*Degree(c):c in cs*} eq {*1^^16,4^^4*};
for c in cs do
 sig:=[Valuation(f,c):f in fs];
 if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
end for;
sigs:=Sort(Setseq(Keys(groups)));assert sigs eq lift[4];
blocks:=[&+[Divisor(c):c in groups[sig]]:sig in sigs];
ds:=[Degree(o):o in blocks];assert ds eq lift[5];c0:=Support(blocks[Index(ds,1)])[1];
print "M3N15_EXACT_ORBIT7_CUSPS",ds,Cputime(t0);

G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);Free:=FreeAbelianGroup(20);
hg:=[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..20]];
h:=hom<Free->G|hg>;H:=Image(h);assert Invariants(H) eq [3,3,24,96,96];
targets:=[h(Free!v):v in lift[7]];assert #Seqset(targets) eq 12;
xs:=NoncuspidalPlaces(X,2);assert #xs eq 24;pcs:=[m2(Divisor(x)-2*Divisor(c0)):x in xs];
pairs:=[];for i in [1..24] do for j in [i..24] do if pcs[i]+pcs[j] eq targets[1] then Append(~pairs,[i,j]);end if;end for;end for;assert #pairs gt 0;
i,j:=Explode(pairs[1]);points:=[xs[i],xs[j]];ecs:=[*EllipticCurve(X,x):x in points*];
levs:=[*LevelStructure(X,x):x in points*];orb:={G|};imageclasses:=[];
for mat in M3FifteenAutMatrices() do
 a,b,c,d:=Explode(mat);
 for u in [1,2] do
  DD:=DivisorGroup(C)!0;
  for k in [1..2] do
   EE:=ecs[k];P:=levs[k]`P;Q:=levs[k]`Q;Q3:=5*Q;Q5:=3*Q;
   La:=rec<recformat<P:PtEll,Q:PtEll>|P:=a*P+b*Q3,Q:=2*(c*P+d*Q3)+2*u*Q5>;
   y:=ModuliPointFast(X,EE,La);assert Degree(y) eq 2;DD+:=Divisor(y);
  end for;
  cl:=m2(DD-4*Divisor(c0));Include(~orb,cl);Append(~imageclasses,Eltseq(cl));
 end for;
end for;
assert orb eq Seqset(targets);
Write("C3x15/data/M3N15_D4_exact_orbit_arithmetic.txt",Sprint([*
 Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(v):v in hg],
 lift[7],[Eltseq(t):t in targets],imageclasses,pairs[1]*]):Overwrite:=true);
print "M3N15_EXACT_GLOBAL_CUSP_CLASSES_ONE_ORBIT",#targets,#orb,#imageclasses,pairs,Cputime(t0);quit;
