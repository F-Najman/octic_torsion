// examples8.m -- independent verification of the 12 realizations used for G subset Phi(8).
// (a) the ten degree-8 entries of van Hoeij's LowDegreePlaces, parsed directly from the raw source file;
// (b) C5 x C10 over Q(zeta20); (c) C3 x C15 over Q(zeta3)(w).
// For each: rebuild K and E, exact order of the given points, full torsion via reductions at good primes
// p with e(P/p) < p-1 (my own choice of primes; torsion injects there), field data (optimized polynomial,
// discriminant, signature, Galois group 8Tn), j-invariant and CM.  Also re-checks the (p,f,h) table of
// Prop. 3.8 and decides which level-34 / level-42 pairs (K,E) are isomorphic or isogenous.
// Markers: EX, TORS, FIELD, JINV, PAIR, TABLE, EXAMPLES8_DONE.
SetColumns(0); SetMemoryLimit(20*10^9);
if not assigned DoTS then DoTS:=1; elif Type(DoTS) eq MonStgElt then DoTS:=StringToInteger(DoTS); end if;
Q:=Rationals(); Qt<t>:=PolynomialRing(Q);
src:="realizations/data/LowDegreePlaces";
F:=Open(src,"r"); L:=[]; while true do s:=Gets(F); if IsEof(s) then break; end if; Append(~L,s); end while;
oct:=[i: i in [1..#L] | Position(L[i],"degv = 8,") gt 0];
printf "octic source lines: %o\n", oct;
assert oct eq [88,102,103,104,112,170,192,268,269,270];

// --- torsion bound from good reductions -------------------------------------------------------
function redbound(E, K, target : pmax:=200, need:=3)
  OK:=MaximalOrder(K); ai:=aInvariants(E); D:=Discriminant(E); g:=0; used:=[];
  for p in [3..pmax] do
    if not IsPrime(p) then continue; end if;
    for pr in Decomposition(OK,p) do
      P:=pr[1]; e:=pr[2]; f:=Degree(P);
      if e ge p-1 then continue; end if;                      // torsion injects only if e < p-1
      if p^f gt 10^12 then continue; end if;
      if exists{a: a in ai | a ne 0 and Valuation(a,P) lt 0} then continue; end if;
      if Valuation(D,P) ne 0 then continue; end if;
      k,red:=ResidueClassField(P);
      Ek:=EllipticCurve([red(a): a in ai]); n:=#Ek;
      g:=Gcd(g,n); Append(~used,<p,f,e,n,Invariants(AbelianGroup(Ek))>);
      if g eq target and #used ge need then return g, used; end if;
    end for;
  end for;
  return g, used;
end function;
function fielddata(K)
  f:=DefiningPolynomial(AbsoluteField(K)); Ka:=AbsoluteField(K);
  Kopt:=OptimizedRepresentation(Ka); fo:=DefiningPolynomial(Kopt);
  G:=GaloisGroup(Ka); d,n:=TransitiveGroupIdentification(G);
  return fo, Discriminant(MaximalOrder(Ka)), Signature(Ka), Sprintf("%oT%o",n,d), #G;
end function;
function jdata(E, K)
  j:=jInvariant(E); mp:=MinimalPolynomial(j, Q);
  cm,D:=HasComplexMultiplication(E);
  if not cm then D:=0; end if;
  return j, mp, cm, D;
end function;

// --- (a) van Hoeij points ------------------------------------------------------------------------
R<x,y>:=PolynomialRing(Q,2);
tab:=[ [<7,8,5760579>,<13,8,815675619>], [<5,2,34>], [<7,4,2448>,<11,2,102>,<13,4,28900>], [<7,4,2482>,<13,4,28424>],
       [<7,8,5763520>,<11,2,105>], [<7,2,39>], [<3,4,80>,<31,3,29560>], [<41,1,42>], [<7,8,5765886>,<13,2,168>], [<7,8,5765886>,<13,2,168>] ];
EXS:=[* *];
for idx in [1..#oct] do
  s:=L[oct[idx]];
  N:=StringToInteger(Split(Split(s,",")[1],"=")[2]);
  a0:=Position(s,"["); b0:=Position(s,"]"); inner:=s[a0+1..b0-1];
  eqs:=[ Split(u,"=")[1] : u in Split(inner,",") ]; assert #eqs eq 2;
  ex:=eval eqs[1]; exy:=eval eqs[2];
  fx:=UnivariatePolynomial(ex);
  if Degree(exy,y) eq 1 then
    assert Degree(fx) eq 8 and IsIrreducible(fx);
    K<a>:=NumberField(fx); xx:=a;
    c1:=Coefficient(exy,y,1); c0:=Coefficient(exy,y,0);
    yy:=-Evaluate(c0,[xx,0])/Evaluate(c1,[xx,0]);
  else
    assert Degree(fx) eq 4 and IsIrreducible(fx) and Degree(exy,y) eq 2;
    K1<a1>:=NumberField(fx); K1u<u>:=PolynomialRing(K1);
    gy:=&+[ Evaluate(Coefficient(exy,y,k),[a1,0])*u^k : k in [0..2]];
    assert IsIrreducible(gy);
    Krel:=NumberField(gy); K:=AbsoluteField(Krel);
    xx:=K!(Krel!a1); yy:=K!(Krel.1);
    mpa:=MinimalPolynomial(xx+yy);
    printf "EX %o: minimal polynomial of x+y = %o (claimed t^8-15t^7+95t^6-331t^5+741t^4-1100t^3+1009t^2-688t+359: %o)\n", idx, mpa,
       mpa eq t^8-15*t^7+95*t^6-331*t^5+741*t^4-1100*t^3+1009*t^2-688*t+359;
  end if;
  assert Evaluate(ex,[xx,yy]) eq 0 and Evaluate(exy,[xx,yy]) eq 0;
  assert Degree(K) eq 8;
  r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx); s2:=(xx*yy-yy+1)/(xx*yy);
  b:=r*s2*(r-1); c:=s2*(r-1);
  E:=EllipticCurve([1-c,-b,-b,0,0]); P:=E![0,0];
  ordP:=Order(P);
  g,used:=redbound(E,K,N);
  printf "EX %o (source line %o): N=%o, field degree %o, Order(P)=%o, gcd of #E(k_P) over good primes with e<p-1: %o\n", idx, oct[idx], N, Degree(K), ordP, g;
  printf "TORS %o: used (p,f,e,#E,group) %o -> torsion %o\n", idx, used, (g eq N and ordP eq N) select Sprintf("exactly C%o",N) else "NOT DETERMINED";
  // re-check the (p,f,h) table of Prop 3.8: counts at all good primes above p with residue degree f
  OK:=MaximalOrder(K);
  for tr in tab[idx] do
    p,f,h:=Explode(tr); cnts:=[];
    for pr in Decomposition(OK,p) do
      if Degree(pr[1]) ne f then continue; end if;
      if exists{aa: aa in aInvariants(E) | aa ne 0 and Valuation(aa,pr[1]) lt 0} or Valuation(Discriminant(E),pr[1]) ne 0 then Append(~cnts,<"bad",pr[2]>); continue; end if;
      k,red:=ResidueClassField(pr[1]); Append(~cnts, <#EllipticCurve([red(aa): aa in aInvariants(E)]), pr[2]>);
    end for;
    printf "TABLE %o: (p,f,h)=(%o,%o,%o): counts (#E, ram.index) at primes above %o of degree %o: %o -> %o\n", idx, p,f,h,p,f, cnts, exists{u: u in cnts | u[1] cmpeq h and u[2] eq 1} select "CONFIRMED" else "NOT FOUND";
  end for;
  fo,dK,sig,gl,og:=fielddata(K);
  j,mp,cm,D:=jdata(E,K);
  printf "FIELD %o: optimized polynomial %o ; disc %o = %o ; signature %o ; Galois group %o (order %o)\n", idx, fo, dK, Factorization(dK), sig, gl, og;
  printf "JINV %o: j has degree %o over Q; minimal polynomial %o ; CM: %o %o ; source says: %o\n", idx, Degree(mp), mp, cm, cm select D else "", Split(s,"[")[1];
  Append(~EXS, <N, K, E, P, xx, yy>);
end for;

// --- isomorphism / isogeny of the level-34 and level-42 pairs (K,E) --------------------------------
function kerpoly(E, T, d)  // kernel polynomial of the cyclic subgroup <T>, T of exact order d
  K:=BaseRing(E); KX<X>:=PolynomialRing(K); xs:={ (k*T)[1]/(k*T)[3] : k in [1..d-1] };
  return &*[X - xv : xv in xs];
end function;
procedure comparepairs(EXS, I)
  for u in [1..#I] do for v in [u+1..#I] do
    i:=I[u]; j:=I[v]; N,K1,E1,P1:=Explode(EXS[i]); _,K2,E2,P2:=Explode(EXS[j]);
    iso,phi:=IsIsomorphic(K1,K2);
    if not iso then
      printf "PAIR %o,%o: fields NOT isomorphic -> distinct pairs (K,E)\n", i, j; continue;
    end if;
    auts:=Automorphisms(K1); res:=[];
    for sig in auts do
      E1s:=EllipticCurve([phi(sig(a)): a in aInvariants(E1)]);
      P1s:=E1s![phi(sig(P1[1])), phi(sig(P1[2]))];
      isoE:=IsIsomorphic(E1s,E2); isog:=[];
      for d in Divisors(N) do
        if d eq 1 then continue; end if;
        T:=(N div d)*P1s; E3:=IsogenyFromKernel(E1s, kerpoly(E1s,T,d));
        if IsIsomorphic(E3,E2) then Append(~isog,d); end if;
      end for;
      Append(~res, <isoE, isog>);
    end for;
    printf "PAIR %o,%o: fields ISOMORPHIC (%o automorphisms); per identification (E1 iso E2?, degrees d of cyclic isogenies E1 -> E1/<(N/d)P> = E2): %o\n", i, j, #auts, res;
  end for; end for;
end procedure;
comparepairs(EXS, [2,3,4]);
comparepairs(EXS, [8,9,10]);
// all ten: which fields coincide at all
for i in [1..10] do for j in [i+1..10] do
  if IsIsomorphic(EXS[i][2], EXS[j][2]) then printf "FIELDISO entries %o and %o have isomorphic fields\n", i, j; end if;
end for; end for;

// --- (b) C5 x C10 over Q(zeta20) -----------------------------------------------------------------
K<z>:=CyclotomicField(20); i:=z^5; assert i^2 eq -1;
E:=EllipticCurve([1-i,-i,-i,0,0]);
P:=E![-1+z^2+z^3-z^4, -1+z^4-z^7]; Qp:=E![1+z^5,-2];
wp:=WeilPairing(P,2*Qp,5);
printf "EX C5xC10: disc(E)=%o (claimed 11-2i: %o), j=%o, Order(P)=%o, Order(Q)=%o, e5(P,2Q)=%o (=-z^6: %o; primitive: %o)\n",
  Discriminant(E), Discriminant(E) eq 11-2*i, jInvariant(E), Order(P), Order(Qp), wp, wp eq -z^6, wp ne 1;
g,used:=redbound(E,K,50);
printf "TORS C5xC10: (p,f,e,#E,group) %o -> gcd %o -> torsion %o\n", used, g, g eq 50 select "exactly C5 x C10" else "NOT DETERMINED";
for pr in [<41,2>,<61,8>] do p:=pr[1]; Fp:=GF(p); r:=Fp!pr[2];
  ok:=Evaluate(CyclotomicPolynomial(20),r) eq 0;
  redz:=func<aa | Evaluate(Polynomial(Fp,[Fp!cc: cc in Eltseq(aa)]), r)>;
  Ek:=EllipticCurve([redz(aa): aa in aInvariants(E)]);
  printf "TABLE C5xC10: z->%o mod %o is a root of Phi20: %o; #E = %o (claimed 50)\n", pr[2], p, ok, #Ek;
end for;
fo,dK,sig,gl,og:=fielddata(K); j,mp,cm,D:=jdata(E,K);
printf "FIELD C5xC10: optimized polynomial %o ; disc %o = %o ; signature %o ; Galois group %o (order %o)\n", fo, dK, Factorization(dK), sig, gl, og;
printf "JINV C5xC10: minimal polynomial of j %o ; CM: %o %o\n", mp, cm, cm select D else "";
if DoTS eq 1 then T,mT:=TorsionSubgroup(E); printf "TS C5xC10: Magma TorsionSubgroup %o\n", Invariants(T); end if;
EX510:=<K,E>;

// --- (c) C3 x C15 over Q(zeta3)(w) ----------------------------------------------------------------
Fz<zz>:=CyclotomicField(3); Fu<u>:=PolynomialRing(Fz);
f:=u^4-5*u^3+8*u^2-6*u+3; assert IsIrreducible(f);
printf "EX C3xC15: disc(f) = %o (claimed -507)\n", Discriminant(t^4-5*t^3+8*t^2-6*t+3);
Kr<w>:=NumberField(f); K:=AbsoluteField(Kr); w:=K!w; z:=K!(Kr!zz);
a1:=((-7*z-8)*w^3+(23*z+34)*w^2+(-8*z-43)*w)/6 + (-3*z+10)/2;
a3:=((17*z+13)*w^3+(-82*z-71)*w^2+(106*z+119)*w)/12 - 2*z - 5;
xQ:=((-z+1)*w^3+(8*z-2)*w^2+(-17*z-1)*w)/6 + 3*z/2;
yQ:=((z+2)*w^3+(-4*z-8)*w^2+(3*z+6)*w+2*z+3)/4;
E:=EllipticCurve([a1,0,a3,0,0]); P:=E![0,0]; Qp:=E![xQ,yQ];
wp:=WeilPairing(P,5*Qp,3);
printf "EX C3xC15: degree %o, disc(E) ne 0: %o, j=%o (claimed 1331/8: %o), Order(P)=%o, Order(Q)=%o, e3(P,5Q)=%o (=z^2: %o)\n",
  Degree(K), Discriminant(E) ne 0, jInvariant(E), jInvariant(E) eq 1331/8, Order(P), Order(Qp), wp, wp eq z^2;
tau:=-yQ^2/xQ^3;
printf "EX C3xC15: tau=-yQ^2/xQ^3 has degree %o over Q (claim: degree 4 over F, i.e. 8 over Q)\n", Degree(MinimalPolynomial(tau,Q));
g,used:=redbound(E,K,45);
printf "TORS C3xC15: (p,f,e,#E,group) %o -> gcd %o -> torsion %o\n", used, g, g eq 45 select "exactly C3 x C15" else "NOT DETERMINED";
fo,dK,sig,gl,og:=fielddata(K); j,mp,cm,D:=jdata(E,K);
printf "FIELD C3xC15: optimized polynomial %o ; disc %o = %o ; signature %o ; Galois group %o (order %o)\n", fo, dK, Factorization(dK), sig, gl, og;
printf "JINV C3xC15: j = %o ; CM: %o %o\n", jInvariant(E), cm, cm select D else "";
// rational curve with j = 1331/8 and the quadratic twist class of E over K
E0:=MinimalModel(EllipticCurveWithjInvariant(Q!1331/8));
printf "JINV C3xC15: rational curve with j=1331/8 (Magma EllipticCurveWithjInvariant, minimal model) %o conductor %o\n", aInvariants(E0), Conductor(E0);
c4,c6:=Explode(cInvariants(E)); d4,d6:=Explode(cInvariants(E0));
dd:=(c6*d4)/(d6*c4);    // E is the quadratic twist of E0_K by dd (up to squares)
ok:=IsIsomorphic(QuadraticTwist(ChangeRing(E0,K),dd), E); printf "JINV C3xC15: E = twist of E0 by d = %o over K: %o; d is a square in K: %o\n", dd, ok, IsSquare(dd);
cands:=[t: t in [-1000..1000] | t ne 0 and IsSquarefree(t)];
rat:=[t: t in cands | IsSquare(dd*t)];
printf "JINV C3xC15: squarefree integers t in [-1000,1000] with d*t a square in K: %o\n", rat;
if DoTS eq 1 then T,mT:=TorsionSubgroup(E); printf "TS C3xC15: Magma TorsionSubgroup %o\n", Invariants(T); end if;
// --- Magma TorsionSubgroup for the ten van Hoeij curves (independent of the reduction bound) --------
if DoTS eq 1 then
  for idx in [1..10] do t0:=Cputime(); T,mT:=TorsionSubgroup(EXS[idx][3]);
    printf "TS %o: N=%o Magma TorsionSubgroup %o (%o s)\n", idx, EXS[idx][1], Invariants(T), Cputime(t0); end for;
end if;
print "EXAMPLES8_DONE"; quit;
