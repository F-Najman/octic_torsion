/* Fresh positive enumeration and absolute-field arithmetic export.
   Override C420PositiveStem before loading to preserve all frozen outputs. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
if not assigned C420PositiveStem then C420PositiveStem:="C4x20/data/C420_positive_fresh";end if;
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_full_aut_hash_m4_16.m");
Attach("common/fast_full_aut_hash_m4_16_relative.m");
Attach("C4x20/fast_full_aut_hash_m4_20.m");
XX:=MDX11(4,20,GF(9));assert Genus(XX) eq 25;
assert #NoncuspidalPlaces(XX,1) eq 0;
cc:=[#Cusps(XX,d):d in [1..4]];assert cc eq [24,12,0,0];
prior:=eval Read("C4x20/data/preflight_4_20_D4_Aut96_arithmetic.txt");
out:=[* *];
for d in [2..4] do
 xs:=NoncuspidalPlaces(XX,d);assert #xs eq [96,320,888][d-1];
 rs,ll,kk:=PlacesUpToFullAutM4Twenty(XX,xs);
 assert kk eq prior[d-1][2] and ll eq prior[d-1][3];
 FF:=BaseRing(XX);KK:=ext<FF|d>;pp:=MinimalPolynomial(KK.1,PrimeField(KK));
 assert Degree(pp) eq 2*d;
 /* This explicit check certifies that exported Eltseq coordinates use powers
    of precisely the exported absolute generator, including its Fi embedding. */
 for ky in kk do
  for offset in [0,2*d] do
   z:=Evaluate(PolynomialRing(KK)![KK!ky[offset+j]:j in [1..2*d]],KK.1);
   assert [Integers()!a:a in Eltseq(z,PrimeField(KK))] eq ky[offset+1..offset+2*d];
  end for;
 end for;
 oo:=[M4TwentyOrbitKeys(XX,x,KK):x in rs];
 assert oo eq prior[d-1][4];
 Append(~out,[*d,[Integers()!a:a in Coefficients(pp)],kk,ll,oo,[Index(xs,x):x in rs]*]);
 print "C420_POSITIVE_COMPLETE_DEGREE",d,#xs,#rs,Cputime(tt);
end for;
Write(C420PositiveStem cat "_arithmetic.txt",Sprint(out):Overwrite:=true);
print "C420_POSITIVE_COMPLETE nc=96,320,888 reps=2,4,12 cusps=24,12,0,0 complete=true",Cputime(tt);
quit;
