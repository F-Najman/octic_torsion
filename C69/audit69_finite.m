/* Independent finite-field exhaustion for the octic69 branch audit. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt69:=Cputime();
counts69:=[];
for e69 in [1..8] do
 F69:=GF(2^e69);count69:=0;
 for j69 in F69 do
  for E69 in Twists(EllipticCurveFromjInvariant(j69)) do
   if #E69 mod 69 ne 0 then continue;end if;
   count69+:=1;
   assert e69 in [7,8] and j69 ne 0;
   assert Degree(MinimalPolynomial(j69,GF(2))) eq e69;
   A69,mp69:=AbelianGroup(E69);
   assert Invariants(A69) eq [e69 eq 7 select 138 else 276];
  end for;
 end for;
 Append(~counts69,count69);
 printf "INDEPENDENT69_FINITE_FIELD e=%o qualifying_curves=%o\n",e69,count69;
end for;
assert counts69 eq [0,0,0,0,0,0,21,16];
printf "INDEPENDENT69_FINITE_DONE counts=%o cyclic_groups=true j_fields_exact=true time=%o\n",counts69,Cputime(tt69);
quit;
