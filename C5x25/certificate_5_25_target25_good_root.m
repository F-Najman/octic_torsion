/* Complete finite good-branch certificate; rank and cuspidal lemmas separate. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
k<z>:=GF(16);K<a>:=GF(256);X:=MDX1(25,k);C:=Curve(X);
assert Genus(X) eq 12;
function Key25(E,P)
 xy:=TateNormalForm_xy(E,P);
 return Minimum([&cat[Eltseq((K!u)^(16^i),GF(2)):u in xy]:i in [0,1]]);
end function;
lookup:=AssociativeArray();curves:=0;points:=0;js:=[];
for j in K do
 for E in Twists(EllipticCurveFromjInvariant(j)) do
  if #E ne 250 then continue;end if;
  A,mp:=AbelianGroup(E);if Invariants(A) ne [5,50] then continue;end if;
  curves+:=1;Append(~js,j);S:=mp(A.1);R:=2*mp(A.2);
  assert Order(S) eq 5 and Order(R) eq 25;
  for d in [d:d in [1..24]|GCD(d,25) eq 1] do for b in [0..4] do
   P:=d*R+b*S;assert Order(P) eq 25;key:=Key25(E,P);
   if not IsDefined(lookup,key) then lookup[key]:=<E,P>;end if;
   points+:=1;
  end for;end for;
 end for;
end for;
printf "M5N25_TARGET25_COMPLETE_FINITE_SCAN curves=%o rawpoints=%o quadratic_places=%o time=%o\n",curves,points,#Keys(lookup),Cputime(tt);
left:=Keys(lookup);reps:=[**];
while #left gt 0 do
 key:=Minimum(Setseq(left));E,P:=Explode(lookup[key]);
 orb:={Key25(E,d*P):d in [d:d in [1..12]|GCD(d,25) eq 1]};
 assert orb subset Keys(lookup);left diff:=orb;Append(~reps,<E,P>);
end while;
printf "M5N25_TARGET25_DIAMOND_ORBITS count=%o time=%o\n",#reps,Cputime(tt);
function EncodeElement25(u)
 return &+[Integers()!v[i]*2^(i-1):i in [1..#v]] where v is Eltseq(K!u,GF(2));
end function;
function Record25(E,P)
 return [EncodeElement25(u):u in aInvariants(E) cat Eltseq(P)[1..2] cat TateNormalForm_bc(E,P) cat TateNormalForm_xy(E,P)];
end function;
allrecords:=[Record25(lookup[key][1],lookup[key][2]):key in Sort(Setseq(Keys(lookup)))];
reprecords:=[Record25(v[1],v[2]):v in reps];
assert #allrecords eq 100 and #reprecords eq 12;
fmt:=recformat<P:PtEll>;remaining:=[];tests:=[**];
for i in [1..#reps] do
 E,P:=Explode(reps[i]);x:=ModuliPoint(X,E,rec<fmt|P:=P>);assert Degree(x) eq 2;
 assert ModuliPoint(X,EllipticCurve(X,x),LevelStructure(X,x)) eq x;
 print "M5N25_TARGET25_TEST_START",i,Cputime(tt);
 H:=HeckeOperator(X,11,x);assert Degree(H) eq 24;
 B:=H-Divisor(x)-11*Divisor(DiamondOperator(X,11,x));assert Degree(B) eq 0;
 BB:=DiamondOperator(X,6,B)-B;
 principal:=IsPrincipal(BB);if principal then Append(~remaining,i);end if;
 Append(~tests,[*i,Degree(x),Degree(H),[Degree(y):y in Support(H)],principal select 1 else 0*]);
 printf "M5N25_TARGET25_ALPHA_A11 index=%o principal=%o time=%o\n",i,principal,Cputime(tt);
end for;
assert remaining eq [];
Write("C5x25/data/M5N25_target25_good_arithmetic_root.txt",Sprint([*Coefficients(DefiningPolynomial(K,GF(2))),allrecords,reprecords,tests*]):Overwrite:=true);
printf "M5N25_TARGET25_GOOD_DIAGNOSTIC_DONE tested=%o survivors=%o time=%o\n",#reps,remaining,Cputime(tt);quit;
