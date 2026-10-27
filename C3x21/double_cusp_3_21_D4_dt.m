/* Fresh native-AF direct local constraints; fallback only if the original Hensel execution is incomplete. */
/* Exact double-cuspidal RR test via intersections of differential spaces.
   For every cusp c and admissible multiplicity m, precompute Omega(2mc)
   inside the common space of holomorphic differentials.  Then
       l(2C) = 2D + 1 - g + dim(intersection_c Omega(2m_c c)).
   This avoids repeating function-field RR calculations for every C.
   The first DirectCheck cases are checked independently with ordinary RR;
   every exceptional case is also checked with ordinary RR. */
SetColumns(0);
SetMemoryLimit(8*10^9);SetQuitOnError(true);
if not assigned N then N:=21;
elif Type(N) eq MonStgElt then N:=StringToInteger(N); end if;
if not assigned D then D:=4;
elif Type(D) eq MonStgElt then D:=StringToInteger(D); end if;
if not assigned Limit then Limit:=0;
elif Type(Limit) eq MonStgElt then Limit:=StringToInteger(Limit); end if;
if not assigned DirectCheck then DirectCheck:=5;
elif Type(DirectCheck) eq MonStgElt then DirectCheck:=StringToInteger(DirectCheck); end if;
if not assigned UseJets then UseJets:=1;
elif Type(UseJets) eq MonStgElt then UseJets:=StringToInteger(UseJets); end if;
AttachSpec("mdmagma/v2/mdmagma.spec");
load "C3x21/dt_fast_differential_expansions.m";
t0:=Cputime();
F:=GF(4);
X:=MDX11(3,N,F);
Sch:=Curve(X);FF:=FunctionField(Sch);C,afmap:=AlgorithmicFunctionField(FF);
g:=Genus(C);assert g eq 25;
jaf:=afmap(jInvariant(EllipticCurve(X`_E)));assert Degree(jaf) eq 576;
cusps:=[c:c in Poles(jaf)|Degree(c) le D];
degs:=[Degree(c):c in cusps];
printf "CUSP_DOUBLE_DIRECT_SETUP N=%o genus=%o cusp_degrees=%o time=%o\n",
 N,g,degs,Cputime(t0);

function WeightedVectors(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else []; end if;
 result:=[];
 for k in [0..d div ds[1]] do
  for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do
   Append(~result,[k] cat v);
  end for;
 end for;
 return result;
end function;
vs:=WeightedVectors(degs,D);
printf "CUSP_DOUBLE_DIRECT_COUNT N=%o effective_divisors=%o time=%o\n",N,#vs,Cputime(t0);
Omega,om:=SpaceOfHolomorphicDifferentials(C);
assert Dimension(Omega) eq g;
W:=VectorSpace(F,g);
localspaces:=AssociativeArray();constraints:=AssociativeArray();
function EncodeF4(z)
 a:=Eltseq(z,GF(2));assert #a eq 2;return Integers()!a[1]+2*(Integers()!a[2]);
end function;
holomorphic:=[om(b):b in Basis(Omega)];
for i in [1..#cusps] do
 if UseJets eq 1 then
  maxmult:=D div degs[i];
  es:=DifferentialExpansionMatrices(holomorphic,cusps[i],2*maxmult);
  assert #es eq 2*maxmult;
  for m in [1..maxmult] do
   B:=HorizontalJoin([es[j]:j in [1..2*m]]);
   assert BaseRing(B) cmpeq F and Nrows(B) eq g and Ncols(B) eq 2*m*degs[i];
   U:=Nullspace(B);constraints[<i,m>]:=[[EncodeF4(a):a in Eltseq(row)]:row in Rows(B)];
   localspaces[<i,m>]:=sub<W|[W!Eltseq(b):b in Basis(U)]>;
   // Independent ordinary RR verifies each individual local kernel dimension.
   rr:=RiemannRochSpace(2*m*Divisor(cusps[i]));
   assert Dimension(rr) eq 2*m*degs[i]+1-g+Dimension(U);
  end for;
 else
  for m in [1..D div degs[i]] do
   U,um:=DifferentialSpace(2*m*Divisor(cusps[i]));
   rows:=[W!Eltseq(Omega!um(b)):b in Basis(U)];
   localspaces[<i,m>]:=sub<W|rows>;
   assert Dimension(localspaces[<i,m>]) eq Dimension(U);
  end for;
 end if;
 printf "CUSP_DOUBLE_DIRECT_LOCAL N=%o cusp=%o/%o time=%o\n",N,i,#cusps,Cputime(t0);
end for;
bad:=[];
tested:=0;
for v in vs do
 if Limit gt 0 and tested ge Limit then break; end if;
 I:=W;
 for i in [1..#cusps] do
  if v[i] gt 0 then I meet:=localspaces[<i,v[i]>]; end if;
 end for;
 ell:=2*D+1-g+Dimension(I);
 assert ell ge 1;
 if tested lt DirectCheck or ell ne 1 then
  E:=&+[v[i]*Divisor(cusps[i]):i in [1..#cusps]];
  assert Degree(E) eq D;
  U:=RiemannRochSpace(2*E);
  assert ell eq Dimension(U);
 end if;
 if ell ne 1 then
  if N eq 15 then
   assert &and[degs[i] eq 1 : i in [1..#cusps] | v[i] gt 0];
   assert Dimension(RiemannRochSpace(E)) eq 1;
  end if;
  Append(~bad,v);
  printf "CUSP_DOUBLE_DIRECT_EXCEPTION N=%o index=%o dim=%o vector=%o\n",N,tested+1,ell,v;
 end if;
 tested+:=1;
 if tested mod 1000 eq 0 then
  printf "CUSP_DOUBLE_DIRECT_PROGRESS N=%o tested=%o exceptions=%o time=%o\n",N,tested,#bad,Cputime(t0);
 end if;
end for;
if N eq 15 then
 print "ALL_EXCEPTIONS_SUPPORTED_ON_RATIONAL_CUSPS_WITH_H0_EQ_1",#bad;
end if;
keys:=[];for i in [1..#cusps] do for m in [1..D div degs[i]] do Append(~keys,<i,m>);end for;end for;
Write("C3x21/data/double_cusp_3_21_D4_dt_arithmetic.txt",Sprint([*
 N,D,g,degs,[[k[1],k[2]]:k in keys],[constraints[k]:k in keys],bad,tested,#vs*]):Overwrite:=true);
printf "CUSP_DOUBLE_DIRECT_DONE N=%o D=%o tested=%o total=%o exceptions=%o complete=%o time=%o\n",
 N,D,tested,#vs,#bad,tested eq #vs,Cputime(t0);
quit;
