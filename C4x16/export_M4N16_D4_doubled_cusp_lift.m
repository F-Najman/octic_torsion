/* Pure exact finite-group arithmetic: no new ClassGroup computation. */
SetColumns(0);SetQuitOnError(true);t0:=Cputime();
a:=eval Read("C4x16/data/diagonal_4_16_D4_p5_p3_arithmetic.txt");
sigs:=a[1];ds:=a[2];rels:=a[4];cg:=a[7][1];source:=a[8];pairs:=a[10];
rows:=[ss:ss in pairs|ss[1] eq 29];assert #rows eq 1 and #rows[1][2] eq 1;idx:=rows[1][2][1];
entry:=source[idx];assert entry[1] eq 4 and entry[2] eq 31 and entry[3] eq 0;
F:=FreeAbelianGroup(#cg[1]);G,pi:=quo<F|[F!r:r in rels]>;
Free:=FreeAbelianGroup(24);h:=hom<Free->G|[pi(F!v):v in cg]>;
assert Invariants(Image(h)) eq [2,2,2,4,4,20,40,40,80];target:=2*pi(F!entry[4]);assert target in Image(h);
v:=Eltseq(target@@h);KB:=Matrix(Integers(),[Eltseq(Free!w):w in Generators(Kernel(h))]);LL:=LatticeWithBasis(LLL(KB));
w:=ClosestVector(LL,Vector(Rationals(),v));short:=[v[i]-Integers()!w[i]:i in [1..24]];assert h(Free!short) eq target;
print "M4N16_D4_DOUBLED_GLOBAL_CUSP_LIFT",short,"degrees",ds,Cputime(t0);
Write("C4x16/data/M4N16_D4_doubled_cusp_lift.txt",Sprint([*4,16,8,sigs,ds,Index(ds,1),short,idx,entry*]):Overwrite:=true);
print "M4N16_D4_DOUBLED_CUSP_LIFT_DONE",Cputime(t0);quit;
