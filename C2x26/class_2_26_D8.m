/* Complete fresh octic necessary rational-image test at3.
   NCmin5, so nonempty NC supports are a single closed place ofdegree5..8.
   For degree5 the complement may be ONE cubic cusp (six possibilities)
   or three rational cusps. Higher NC degrees have only rational cusp
   complements. The all-cusp branch is excluded separately by odd-prime
   specialization, rank0 and gon_Q>=12>8. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");
Attach("common/fast_diamond_hash_m2.m");
X:=MDX11(2,26,GF(3));C:=Curve(X);assert Genus(X) eq 25;
assert X`_E[1] eq 0 and X`_E[3] eq 0;
cs:=Cusps(X);assert {*Degree(c):c in cs*} eq {*1^^18,3^^6*};
rat:=[c:c in cs|Degree(c) eq 1];cub:=[c:c in cs|Degree(c) eq 3];c0:=rat[1];
places:=[NoncuspidalPlaces(X,d):d in [1..8]];
assert [#ss:ss in places[1..7]] eq [0,0,0,0,36,122,288];
J:=JH(52,[27]);assert Dimension(J) eq 25;fp:=FrobeniusPolynomial(J,3);
assert Degree(fp) eq 50;traces:=[];closed:=[];
for n in [1..8] do
 sn:=-n*Coefficient(fp,50-n);
 for j in [1..n-1] do sn-:=Coefficient(fp,50-j)*traces[n-j];end for;
 Append(~traces,Integers()!sn);count:=3^n+1-Integers()!sn;
 for d in Divisors(n) do if d lt n then count-:=d*closed[d];end if;end for;
 assert count mod n eq 0;Append(~closed,count div n);
 assert closed[n] eq #places[n]+#[c:c in cs|Degree(c) eq n];
end for;
reps:=[PlacesUpToDiamondM2HashDirect(X,places[d]):d in [5..8]];
assert [#ss:ss in reps[1..3]] eq [6,26,48];
print "OCTIC_N26_CLASS_SETUP",[#ss:ss in places],[#rs:rs in reps],closed,Cputime(t0);
print "OCTIC_N26_CLASS_GROUP_START",Cputime(t0);
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
assert &*[n:n in Invariants(G)|n ne 0] eq Integers()!Evaluate(fp,1);
hgens:=[m2(Divisor(c)-Divisor(c0)):c in rat];H:=sub<G|hgens>;
assert Invariants(H) eq [2,14,266,3990,11970];
cc:=[m2(Divisor(c)-3*Divisor(c0)):c in cub];
print "OCTIC_N26_CLASS_GROUP",Invariants(G),Invariants(H),[Eltseq(a):a in cc],Cputime(t0);
records:=[* *];remaining:=[* *];tested:=0;
for pos in [1..4] do
 d:=pos+4;passed:=0;
 for i in [1..#reps[pos]] do
  x:=reps[pos][i];cl:=m2(Divisor(x)-d*Divisor(c0));matches:=[];
  if cl in H then Append(~matches,0);end if;
  if d eq 5 then
   for j in [1..#cc] do if cl+cc[j] in H then Append(~matches,j);end if;end for;
  end if;
  inside:=#matches gt 0;tested+:=1;
  Append(~records,[*d,i,Eltseq(cl),inside,matches*]);
  if inside then passed+:=1;Append(~remaining,[*d,i,x,matches*]);end if;
  print "OCTIC_N26_CLASS_TEST",d,i,inside,matches,Eltseq(cl),Cputime(t0);
 end for;
 print "OCTIC_N26_CLASS_TYPE",d,#reps[pos],passed,Cputime(t0);
end for;
assert tested eq &+[#rs:rs in reps];export:=[* *];
for rec in remaining do
 d,i,x,matches:=Explode(rec);pt:=RepresentativePoint(x);L:=GF(3^d);
 coords:=MDM2ToCommonField(Eltseq(pt),L);
 def:=[Integers()!a:a in Coefficients(DefiningPolynomial(L,GF(3)))];
 values:=[[Integers()!a:a in Eltseq(z,GF(3))]:z in coords];
 Append(~export,[*d,i,def,values,matches*]);
end for;
Write("C2x26/data/N26_D8_class_arithmetic.txt",Sprint([*
 Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],
 [Eltseq(g):g in hgens],[Eltseq(a):a in cc],records,Coefficients(fp),
 [#ss:ss in places],[#rs:rs in reps]*]):Overwrite:=true);
Write("C2x26/data/N26_D8_class_survivors.m",
 "N26D8ClassSurvivors := " cat Sprint(export) cat ";":Overwrite:=true);
printf "OCTIC_N26_CLASS_DONE tested=%o remaining=%o complete=true time=%o\n",tested,#remaining,Cputime(t0);quit;
