Attach("common/fast_hecke.m");
Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
tl54:=Cputime();allps54:=[[]:d in [1..6]];good54:=[[]:d in [1..6]];inert54:=[];root54:=[];
for d in [2..6] do
 allps54[d]:=NoncuspidalPlaces(X5,d);nrlist:=[];
 for P in allps54[d] do
  EP:=EllipticCurve(X5,P);b2e,b4e,b6e,b8e:=Explode(bInvariants(EP));RE<T>:=PolynomialRing(BaseRing(EP));
  nr:=#Roots(4*T^3+b2e*T^2+2*b4e*T+b6e);assert nr in [0,1,3];Append(~nrlist,nr);
  if nr gt 0 then Append(~good54[d],P);elif d eq 2 then Append(~inert54,P);end if;
 end for;
 Append(~root54,nrlist);print "OCTIC54_LOWER_POINTS",d,#allps54[d],#good54[d],Cputime(tl54);
end for;
assert #good54[2] eq 0 and #inert54 eq 9;
checked54:=0;members54:=0;passing54:=0;hits54:=[];
testData54:="C54/data/C54_lower_classes.txt";
procedure TestClass54(a,e,desc,~checked,~members,~passing,~hits)
 for cc in comps[8-e] do
  q:=a+cc;checked+:=1;member:=q in H5;passes:=false;
  if member then
   members+:=1;v:=q@@h5;assert h5(v) eq q;
   if h2(v) in candidate2 then passes:=true;passing+:=1;Append(~hits,<desc,Eltseq(v)>);end if;
  end if;
  Write(testData54,Sprint([Integers()!z:z in Eltseq(q)] cat [member select 1 else 0,passes select 1 else 0]));
 end for;
end procedure;
for d in [3..6] do
 reps:=PlacesUpToDiamondFast(X5,good54[d]);print "OCTIC54_LOWER_SINGLE_START",d,#reps,Cputime(tl54);
 for i in [1..#reps] do TestClass54(b5(Divisor(reps[i])-d*Divisor(c05)),d,<d,i>,~checked54,~members54,~passing54,~hits54);end for;
 print "OCTIC54_LOWER_SINGLE_DONE",d,checked54,members54,passing54,Cputime(tl54);
end for;
for P in PlacesUpToDiamondFast(X5,inert54) do TestClass54(b5(3*Divisor(P)-6*Divisor(c05)),6,<2,3>,~checked54,~members54,~passing54,~hits54);end for;
lowerQ54:=AssociativeArray();for d in [3,4,5] do lowerQ54[d]:=[b5(Divisor(P)-d*Divisor(c05)):P in good54[d]];end for;
for pair in [<3,3>,<3,4>,<3,5>,<4,4>] do
 d,e:=Explode(pair);reps:=PlacesUpToDiamondFast(X5,good54[d]);
 for i in [1..#reps] do
  a:=b5(Divisor(reps[i])-d*Divisor(c05));
  for b in lowerQ54[e] do TestClass54(a+b,d+e,<d,e,i>,~checked54,~members54,~passing54,~hits54);end for;
 end for;
 print "OCTIC54_LOWER_MIXED_DONE",pair,checked54,members54,passing54,Cputime(tl54);
end for;
printf "OCTIC54_LOWER_SUPPORTS_DONE checked=%o members=%o passing=%o hits=%o complete=true time=%o\n",checked54,members54,passing54,hits54,Cputime(tl54);
