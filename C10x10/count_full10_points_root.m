/* Independent moduli point counts, without the genus13 function field. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
for p in [11,31] do
 for d in [1,2] do
  k:=GF(p^d);mass:=Rationals()!0;curves:=0;groups:={@ @};
  for j in k do
   for E in Twists(EllipticCurveFromjInvariant(j)) do
    if #E mod 100 ne 0 then continue;end if;
    A:=AbelianGroup(E);iv:=Invariants(A);
    if #iv eq 2 and iv[1] mod 10 eq 0 and iv[2] mod 10 eq 0 then
     aut:=#Automorphisms(E);assert aut in [2,4,6];
     mass+:=720/aut;curves+:=1;Include(~groups,<iv,aut>);
    end if;
   end for;
  end for;
  expected:=p eq 11 select 300 else 1620;
  assert mass eq (d eq 1 select 0 else expected);
  print "FULL10_INDEPENDENT_MODULI_COUNT",p,d,"noncuspidal_points",mass,"curves",curves,"groups",groups,Cputime(tt);
 end for;
end for;
print "FULL10_INDEPENDENT_MODULI_COUNTS_DONE";quit;
