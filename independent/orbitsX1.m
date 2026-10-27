// orbitsX1.m -- closed points of X1(N) over F_q of exact degree e and their diamond orbits, by explicit
// enumeration: every curve over F_{q^e} (one j per Frobenius orbit, all twists), every point P of exact
// order N, Tate normal form (b,c) of (E,P) as canonical coordinate, closed-point key = minimum over the
// Frobenius conjugates, diamond <a> acting by (E,P) -> (E,aP).  Parameters q, D, Ns="45,55", optional M
// (restrict to curves with E[M] rational, e.g. M:=5 for X1(25) over F16).  Markers ORB, STAB, ORBITS_DONE.
SetColumns(0); SetMemoryLimit(12*10^9);
if Type(q) eq MonStgElt then q:=StringToInteger(q); end if;
if Type(D) eq MonStgElt then D:=StringToInteger(D); end if;
if not assigned Dmin then Dmin:=1; elif Type(Dmin) eq MonStgElt then Dmin:=StringToInteger(Dmin); end if;
if not assigned M then M:=1; elif Type(M) eq MonStgElt then M:=StringToInteger(M); end if;
Nlist:=[StringToInteger(s): s in Split(Ns,",")];
p:=Factorization(q)[1][1];
function tate(E,Q)
  k:=BaseRing(E); a1,a2,a3,a4,a6:=Explode(aInvariants(E)); x0:=Q[1]; y0:=Q[2];
  R<X,Y>:=PolynomialRing(k,2); f:=Y^2+a1*X*Y+a3*Y-(X^3+a2*X^2+a4*X+a6);
  g:=Evaluate(f,[X+x0,Y+y0]); assert MonomialCoefficient(g,R!1) eq 0;
  A3:=MonomialCoefficient(g,Y); A4:=-MonomialCoefficient(g,X); assert A3 ne 0;
  g:=Evaluate(g,[X,Y+(A4/A3)*X]); assert MonomialCoefficient(g,X) eq 0;
  A1:=MonomialCoefficient(g,X*Y); A2:=-MonomialCoefficient(g,X^2); A3:=MonomialCoefficient(g,Y);
  return -A2^3/A3^2, 1-A1*A2/A3;
end function;
function ptdeg(b,c,e) for f in Divisors(e) do if b^(q^f) eq b and c^(q^f) eq c then return f; end if; end for; end function;
function key(b,c,e)
  return Min([ [Integers()!z : z in Eltseq(b^(q^i),GF(p))] cat [Integers()!z : z in Eltseq(c^(q^i),GF(p))] : i in [0..e-1]]); end function;
for N in Nlist do
  dia:=[a: a in [1..N div 2] | Gcd(a,N) eq 1];
  for e in [Dmin..D] do
    K:=GF(q^e); visited:={K|}; seen:={}; orbs:=[]; stabs:=[]; jorb:={* *};
    for j in K do
      if j in visited then continue; end if;
      conj:={j^(q^i): i in [0..e-1]}; visited join:=conj;
      E0:=EllipticCurveWithjInvariant(j);
      Ts:=(j eq 0 or j eq K!1728) select Twists(E0) else [E0, QuadraticTwist(E0)];
      for T in Ts do
        n:=#T; if n mod (N*M) ne 0 then continue; end if;
        G,mG:=AbelianGroup(T); inv:=Invariants(G);
        if M gt 1 and (#inv lt 2 or inv[1] mod M ne 0) then continue; end if;
        HN:=sub<G | [(inv[i] div Gcd(inv[i],N))*G.i : i in [1..#inv]]>;
        for g in HN do
          if Order(g) ne N then continue; end if;
          P:=mG(g); b,c:=tate(T,P);
          if ptdeg(b,c,e) ne e then continue; end if;
          k0:=key(b,c,e); if k0 in seen then continue; end if;
          orb:={}; st:=[];
          for a in dia do b2,c2:=tate(T,a*P); ka:=key(b2,c2,e); Include(~orb,ka); if ka eq k0 then Append(~st,a); end if; end for;
          seen join:=orb; Append(~orbs,#orb); Append(~stabs,st); Include(~jorb, <#conj, n, inv>);
        end for;
      end for;
    end for;
    printf "ORB X1(%o)%o over F_%o, degree %o: %o closed points, %o diamond orbits, orbit sizes %o\n", N, M gt 1 select Sprintf("[E[%o] rational]",M) else "", q, e, #seen, #orbs, {* x: x in orbs *};
    if #orbs gt 0 and #orbs le 30 then printf "STAB X1(%o) F_%o deg %o: diamond stabilizers (reps a in 1..N/2) per orbit %o ; curves (j-orbit length, #E, group) %o\n", N, q, e, stabs, jorb; end if;
  end for;
end for;
print "ORBITS_DONE"; quit;
