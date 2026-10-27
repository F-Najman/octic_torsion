/* Fresh octic55 NC7+cusp1 branch only. Independent exact counts exclude other noncuspidal degrees through8. All F2-rational source cusps liftQ and are diamond translates of the rational Neron55-gon cusp, so A3 kills the complement as a divisor. Allbad2 remains separate. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");
Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
X:=MDX1(55,GF(2));assert Genus(X) eq 81;
places:=NoncuspidalPlaces(X,7);
reps:=PlacesUpToDiamondFast(X,places);
print "N55_D8_NC7_PLUS_CUSP_SETUP",#places,#reps,Cputime(t0);
remaining:=0;
for i in [1..#reps] do
    print "N55_D8_NC7_PLUS_CUSP_CASE_START",i,#reps,Cputime(t0);
    x:=reps[i];T:=HeckeOperatorCommonZeros(X,3,Divisor(x));
    A:=T-3*DiamondOperator(X,3,Divisor(x))-Divisor(x);
    assert Degree(A) eq 0;
    principal:=IsPrincipal(A);
    print "N55_D8_NC7_PLUS_CUSP_CASE",i,"principal",principal,Cputime(t0);
    if principal then remaining+:=1;end if;
end for;
printf "N55_D8_NC7_PLUS_CUSP_DONE total=%o remaining=%o complete=true time=%o\n",#reps,remaining,Cputime(t0);
quit;
