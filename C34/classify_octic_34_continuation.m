/* Load only after classify_octic_34_setup.m reaches SETUP_DONE. */
assert assigned G and assigned HQ and assigned places;
assert Genus(X) eq 21 and #places[4] eq 26 and #places[8] eq 734;
out34c:=Open("C34/logs/classify_octic_34_continuation.log","w");
procedure Note34C(s)
 fprintf out34c,"%o\n",s;Flush(out34c);print s;
end procedure;
load "C34/known_34_reductions_v2.m";
known34classes:={m2(dd34-8*Divisor(c0)):dd34 in known34divisors};
assert #known34classes eq 7 and &and[a34 in HQ:a34 in known34classes];
Note34C(Sprintf("C34_D8_KNOWN_DIVISORS classes=%o weighted_reduction_patterns=%o time=%o",#known34classes,known34sources,Cputime(t34)));
cc34:=AssociativeArray();for j34 in [0..4] do cc34[j34]:={G|};end for;
cc34[0]:={G!0};
for cusp34 in [c:c in cs|Degree(c) eq 1] do
 a34:=m2(Divisor(cusp34)-Divisor(c0));assert a34 in HQ;
 for j34 in [1..4] do cc34[j34] join:={z34+a34:z34 in cc34[j34-1]};end for;
end for;
ar34:=Open("C34/data/class_arithmetic.py","w");
fprintf ar34,"relations = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G))];
fprintf ar34,"rational_cusp_generators = %o\n",[Eltseq(a):a in qc];
fprintf ar34,"rational_cusp_degrees = %o\n",[Degree(o):o in os];
fprintf ar34,"known_octic_classes = %o\n",[Eltseq(a):a in known34classes];
fprintf ar34,"known_reduction_prime_patterns = %o\n",known34sources;
fprintf ar34,"frobenius_coefficients = %o\n",Coefficients(fp);
fprintf ar34,"tests = [\n";
nc34classes:=AssociativeArray();counts34:=[];passes34:=[];candidates34:={G|};
for d34 in [4..8] do
 aa34:=[G|];pass34:=0;
 for i34 in [1..#places[d34]] do
  a34:=m2(Divisor(places[d34][i34])-d34*Divisor(c0));Append(~aa34,a34);
  member34:=a34 in HQ;
  if member34 then
   pass34+:=1;
   if d34 eq 8 then Include(~candidates34,a34);
   else candidates34 join:={a34+z34:z34 in cc34[8-d34]};end if;
  end if;
  fprintf ar34,"[%o,%o,%o,%o],\n",d34,i34,Eltseq(a34),member34 select 1 else 0;
  Note34C(Sprintf("C34_D8_CLASS_TEST degree=%o index=%o/%o in_rational_image=%o time=%o",d34,i34,#places[d34],member34,Cputime(t34)));
 end for;
 nc34classes[d34]:=aa34;Append(~counts34,#aa34);Append(~passes34,pass34);
end for;
fprintf ar34,"]\nquartic_pair_tests = [\n";
pair_count34:=0;pair_pass34:=0;
for i34 in [1..26] do for j34 in [i34..26] do
 a34:=nc34classes[4][i34]+nc34classes[4][j34];member34:=a34 in HQ;
 pair_count34+:=1;
 if member34 then pair_pass34+:=1;Include(~candidates34,a34);end if;
 fprintf ar34,"[%o,%o,%o,%o],\n",i34,j34,Eltseq(a34),member34 select 1 else 0;
 Note34C(Sprintf("C34_D8_QUARTIC_PAIR i=%o j=%o in_rational_image=%o",i34,j34,member34));
end for;end for;
assert pair_count34 eq 351;
fprintf ar34,"]\ncomplete = True\n";delete ar34;
unknown34:=candidates34 diff known34classes;
assert known34classes subset candidates34;
complete34:=#unknown34 eq 0;
Note34C(Sprintf("C34_OCTIC_CLASSIFICATION_DONE single_counts=%o single_members=%o quartic_pairs=%o pair_members=%o candidate_classes=%o known_octic_classes=7 unknown_octic_classes=%o complete=%o time=%o",counts34,passes34,pair_count34,pair_pass34,#candidates34,#unknown34,complete34,Cputime(t34)));
if complete34 then delete out34c;delete out34;quit;end if;
