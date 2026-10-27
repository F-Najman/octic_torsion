/* Independent exhaustive finite-curve enumeration, not a modular point formula. */
SetColumns(0);SetQuitOnError(true);tt:=Cputime();
for d in [1..8] do
 F:=GF(2^d);good:=0;orders:={};groups:={};n:=0;
 for j in F do
  for E in Twists(EllipticCurveFromjInvariant(j)) do
   n+:=1;
   if #E mod 51 eq 0 then
    inv:=Invariants(AbelianGroup(E));Include(~orders,#E);Include(~groups,inv);
    if #inv gt 0 and inv[#inv] mod 51 eq 0 then good+:=1;end if;
   end if;
  end for;
 end for;
 assert good eq 0;
 print "N51_EXHAUSTIVE_CURVES",d,"twists",n,"orders_divisible51",orders,"groups",groups,"curves_with_order51_point",good,Cputime(tt);
end for;
print "N51_INDEPENDENT_CURVE_ENUMERATION_DONE";quit;
