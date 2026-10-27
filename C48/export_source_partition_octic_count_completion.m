/* Honest octic-only restart after the preserved16GB partition failure.
   Degrees3..7 are the verified completed prefix. All completion writes
   are INSIDE one compiled procedure, so an error cannot fall through. */
load "C48/common_field_by_coefficients.m";
SetMemoryLimit(128*10^9);
assert Keys(sourceRepBC48) eq {3,4,5,6,7};
assert [#sourceRepBC48[d]:d in [3..7]] eq [2,8,82,342,1348];
assert #p8f48 eq 48992 and FR`cld`proven;
tf48:=Cputime();
of48v2:=Open("C48/logs/export_source_partition_octic_count_completion.log","w");
af48:=Open("C48/data/source_partition_inputs_complete.py","a");
procedure PartNote48V3(s)
 fprintf of48v2,"%o\n",s;Flush(of48v2);print s;
end procedure;
procedure C48CompleteOcticByCount(~sourceRepBC48,~sourceAllBC48,~sourceOrbitSizes48,~sourceFields48)
for dd48 in [8] do
 K48:=sourceFields48[dd48];
 seen48:=AssociativeArray();inputs48:={@ @};data48:=[];allbc48:=[];orbitSizes48:=[];
 for ip48 in [1..#sourceFunctions48[dd48]] do
  pp48:=sourceFunctions48[dd48][ip48];ok48,bc48:=DirectSourceBC48(pp48);
  if not ok48 then continue;end if;
  bc48:=C48CommonFieldByCoefficients(bc48,K48);key48:=MDM2FrobeniusKey(bc48,K48,5,dd48);
  assert key48 notin inputs48;Include(~inputs48,key48);Append(~allbc48,bc48);
  if not IsDefined(seen48,key48) then
   Append(~data48,bc48);b48,c48:=Explode(bc48);
   EE48:=EllipticCurve([1-c48,-b48,-b48,K48!0,K48!0]);PP48:=EE48![0,0];oks48:={@ @};
   for aa48 in [1..24] do
    if GCD(aa48,48) ne 1 then continue;end if;
    ky48:=MDM2FrobeniusKey(TateNormalForm_bc(EE48,aa48*PP48),K48,5,dd48);
    Include(~oks48,ky48);
    if IsDefined(seen48,ky48) then assert seen48[ky48] eq #data48;
    else seen48[ky48]:=#data48;end if;
   end for;
   Append(~orbitSizes48,#oks48);
  end if;
  assert IsDefined(seen48,key48);
  if #Keys(seen48) eq nc48[dd48] then break;end if;
  if ip48 mod 2000 eq 0 then
   PartNote48V3(Sprintf("C48_PARTITION_PROGRESS degree=%o input=%o total=%o representatives=%o time=%o",dd48,ip48,#sourceFunctions48[dd48],#data48,Cputime(tf48)));
  end if;
 end for;
 assert #inputs48 le nc48[dd48] and #Keys(seen48) eq nc48[dd48];
 assert &+orbitSizes48 eq nc48[dd48];
 sourceRepBC48[dd48]:=data48;sourceAllBC48[dd48]:=allbc48;sourceOrbitSizes48[dd48]:=orbitSizes48;
 fprintf af48,"field_polynomials[%o] = %o\n",dd48,[Integers()!a:a in Coefficients(DefiningPolynomial(K48,GF(5)))];
 fprintf af48,"representative_pairs[%o] = %o\n",dd48,[[[Integers()!a:a in Eltseq(v,GF(5))]:v in bc]:bc in data48];
 fprintf af48,"noncuspidal_counts[%o] = %o\npartition_orbit_sizes[%o] = %o\n",dd48,nc48[dd48],dd48,orbitSizes48;Flush(af48);
 PartNote48V3(Sprintf("C48_PARTITION_DONE degree=%o places=%o representatives=%o orbit_size_sum=%o complete_by_newton_count=true time=%o",dd48,nc48[dd48],#data48,&+orbitSizes48,Cputime(tf48)));
end for;
assert Keys(sourceRepBC48) eq {3,4,5,6,7,8};
assert #sourceAllBC48[8] le 48992 and &+sourceOrbitSizes48[8] eq 48992;
fprintf af48,"octic_original_places_scanned = %o\ncomplete_by_newton_count = True\n",#sourceAllBC48[8];
fprintf af48,"complete = True\n";Flush(af48);
PartNote48V3(Sprintf("C48_ALL_PARTITIONS_COUNT_DONE complete=true disjoint_octic_places=48992 time=%o",Cputime(tf48)));
end procedure;
C48CompleteOcticByCount(~sourceRepBC48,~sourceAllBC48,~sourceOrbitSizes48,~sourceFields48);
