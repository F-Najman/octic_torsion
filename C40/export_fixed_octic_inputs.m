/* Retained proven C40 state.  Preserve the interrupted prefix and export a
   single exact representative table before the native Fork snapshot. */
fprintf ar40,"]\ncomplete = False\ninterrupted_prefix = True\n";delete ar40;
delete out40c;
Attach("common/fast_hecke.m");
Attach("common/fast_hecke_11.m");
assert #reps40[8] eq 96 and #places[8] eq 708;
masterlog40:=Open("C40/logs/fixed_octic_inputs.log","w");
Kmaster40:=Parent(Evaluate(X`_coordinates[1],reps40[8][1]));
assert Degree(Kmaster40) eq 8 and #Kmaster40 eq 3^8;
fmaster40:=[Integers()!a:a in Coefficients(DefiningPolynomial(Kmaster40,GF(3)))];
assert #fmaster40 eq 9;
Kround40:=ext<GF(3)|PolynomialRing(GF(3))!fmaster40>;
fixedbc40:=[];fixedplaceindices40:=[];fixedkeys40:=AssociativeArray();
for ii40 in [1..96] do
 P40:=reps40[8][ii40];
 bc40:=MDM2ToCommonField([Evaluate(f,P40):f in X`_coordinates],Kmaster40);
 coeff40:=[[Integers()!a:a in Eltseq(z,GF(3))]:z in bc40];
 assert #coeff40 eq 2 and &and[#v eq 8:v in coeff40];
 bcr40:=[Evaluate(PolynomialRing(Kround40)!v,Kround40.1):v in coeff40];
 assert MDPlaceCommonZeros(X`_coordinates,bcr40) eq P40;
 Append(~fixedbc40,coeff40);Append(~fixedplaceindices40,Index(places[8],P40));
 for key40 in MDCyclicBCHashKeys(X,P40,Kmaster40) do
  if IsDefined(fixedkeys40,key40) then assert fixedkeys40[key40] eq ii40;
  else fixedkeys40[key40]:=ii40;end if;
 end for;
 fprintf masterlog40,"C40_FIXED_OCTIC_ROUNDTRIP index=%o/96 original_place_index=%o exact=true\n",ii40,fixedplaceindices40[ii40];Flush(masterlog40);
end for;
assert #Keys(fixedkeys40) eq 708;
fixedorbits40:=[];
for P40 in places[8] do
 key40:=MDM2FrobeniusKey([Evaluate(f,P40):f in X`_coordinates],Kmaster40,3,8);
 assert IsDefined(fixedkeys40,key40);Append(~fixedorbits40,fixedkeys40[key40]);
end for;
assert Seqset(fixedorbits40) eq Seqset([1..96]);
fixedout40:=Open("C40/data/fixed_octic_inputs.py","w");
fprintf fixedout40,"field_polynomial = %o\nrepresentatives_bc = %o\noriginal_place_indices = %o\nall_place_orbit_indices = %o\n",fmaster40,fixedbc40,fixedplaceindices40,fixedorbits40;
delete fixedout40;
fixedout40:=Open("C40/data/fixed_octic_inputs.m","w");
fprintf fixedout40,"C40FixedFieldPolynomial := %o;\nC40FixedTateBC := %o;\nC40FixedOriginalPlaceIndices := %o;\nC40FixedAllOrbitIndices := %o;\n",fmaster40,fixedbc40,fixedplaceindices40,fixedorbits40;
delete fixedout40;
fprintf masterlog40,"C40_FIXED_OCTIC_MASTER_DONE representatives=96 original_places=708 all_roundtrips=true partition_complete=true\n";
delete masterlog40;
print "C40_FIXED_OCTIC_MASTER_DONE representatives=96 original_places=708 all_roundtrips=true partition_complete=true";
