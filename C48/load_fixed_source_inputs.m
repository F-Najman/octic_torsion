/* Load only the small fixed integer table; no large source residue fields. */
load "C48/data/fixed_source_inputs.m";
assert C48Degrees eq [3..8] and #C48MasterSHA eq 64;
sourceFields48:=AssociativeArray();sourceRepBC48:=AssociativeArray();
sourceAllBC48:=AssociativeArray();sourceOrbitSizes48:=AssociativeArray();
for ii48 in [1..6] do
 dd48:=C48Degrees[ii48];PP48:=PolynomialRing(GF(5));
 assert IsIrreducible(PP48!C48FieldPolynomials[ii48]);
 KK48:=ext<GF(5)|PP48!C48FieldPolynomials[ii48]>;
 sourceFields48[dd48]:=KK48;
 sourceRepBC48[dd48]:=[ [&+[KK48!v[i]*KK48.1^(i-1):i in [1..#v]]:v in pair]:pair in C48RepresentativePairs[ii48]];
 sourceOrbitSizes48[dd48]:=C48OrbitSizes[ii48];
 assert #sourceRepBC48[dd48] eq #sourceOrbitSizes48[dd48];
 assert &+sourceOrbitSizes48[dd48] eq C48NoncuspidalCounts[ii48];
end for;
for dd48 in [3..5] do
 KK48:=sourceFields48[dd48];allBC:=[];allkeys:={@ @};
 for bc in sourceRepBC48[dd48] do
  b,c:=Explode(bc);EE:=EllipticCurve([1-c,-b,-b,KK48!0,KK48!0]);PP:=EE![0,0];
  assert 48*PP eq EE!0 and 24*PP ne EE!0 and 16*PP ne EE!0;
  for aa in [1..24] do
   if GCD(aa,48) ne 1 then continue;end if;
   pair:=TateNormalForm_bc(EE,aa*PP);key:=MDM2FrobeniusKey(pair,KK48,5,dd48);
   assert key notin allkeys;Include(~allkeys,key);Append(~allBC,pair);
  end for;
 end for;
 assert #allBC eq C48NoncuspidalCounts[dd48-2];sourceAllBC48[dd48]:=allBC;
end for;
print "C48_FIXED_INPUTS_LOADED",[#sourceRepBC48[d]:d in [3..8]],"low_all_counts",[#sourceAllBC48[d]:d in [3..5]],"master_sha",C48MasterSHA;
