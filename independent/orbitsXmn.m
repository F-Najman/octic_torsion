// orbitsXmn.m -- closed points of X1(m,n) (m|n, Weil pairing e_m(P,(n/m)Q) fixed when m>2) over F_q of exact
// degree e, and their orbits under the explicit groups of level-structure automorphisms used in the draft.
// Canonical coordinates of a point (E,P,Q): Tate normal form (b,c) of (E,Q) and the coordinates of P in that
// (rigid) model; closed-point key = minimum over Frobenius conjugates.  Each group element is given as
// integers [al,be,ga,de] meaning (P,Q) -> (al*P + be*Q, ga*P + de*Q); every image is checked to be a level
// structure of the same type (orders, independence, pairing).  Parameters q, D, Dmin, m, n, grp.
// Markers ORBMN, ORBITSMN_DONE.
SetColumns(0); SetMemoryLimit(16*10^9);
if Type(q) eq MonStgElt then q:=StringToInteger(q); end if; if Type(D) eq MonStgElt then D:=StringToInteger(D); end if;
if Type(Dmin) eq MonStgElt then Dmin:=StringToInteger(Dmin); end if;
if Type(m) eq MonStgElt then m:=StringToInteger(m); end if; if Type(n) eq MonStgElt then n:=StringToInteger(n); end if;
p:=Factorization(q)[1][1];
U:=func<N | [d: d in [1..N] | Gcd(d,N) eq 1]>;
// ---- the groups, as in the draft --------------------------------------------------------------------
case grp:
  when "diamond2":   // X1(2,n): (P, dQ), d in (Z/n)^*/+-1
    GR:=[[1,0,0,d]: d in U(n) | d le n div 2];
  when "aut2":       // X1(2,n), 4|n: (P + (n/2)bQ, cP + dQ), b,c in {0,1}, d in (Z/n)^*/+-1   (Lemmas 3.52, 3.56; Prop 3.59)
    GR:=[[1,(n div 2)*b,c,d]: b,c in [0,1], d in U(n) | d le n div 2];
  when "aut2_42":    // X1(2,42): (P,Q) -> (S_i, S_j + 2uQ), S1=P, S2=21Q, S3=P+21Q, i ne j, u in {1,2,4,5,8,10}
    S:=[[1,0],[0,21],[1,21]];
    GR:=[[S[i][1],S[i][2],S[j][1],S[j][2]+2*u]: i,j in [1..3], u in [1,2,4,5,8,10] | i ne j];
  when "aut3_18":    // (aP + 6bQ, cP + dQ), b,c in Z/3, d in {1,5,7}, a = d^-1 mod 3
    GR:=[[Integers()!((GF(3)!d)^-1),6*b,c,d]: b,c in [0..2], d in [1,5,7]];
  when "aut3_21":    // Q3=7Q, Q7=3Q: (aP + bQ3, cP + dQ3 + 5uQ7), (a b; c d) in SL2(F3), u in (Z/7)^*, mod (-1,-1)
    GR:=[[a,7*b,c,7*d+15*u]: a,b,c,d in [0..2], u in [1..6] | (a*d-b*c) mod 3 eq 1];
  when "aut3_24":    // R3=8Q, R8=3Q, Q=2R3+3R8: (aP + bR3, 2(cP + dR3) + 3uR8), SL2(F3), u in (Z/8)^*
    GR:=[[a,8*b,2*c,16*d+9*u]: a,b,c,d in [0..2], u in [1,3,5,7] | (a*d-b*c) mod 3 eq 1];
  when "aut4_16":    // (aP + 4bQ, cP + dQ), b,c in Z/4, d in (Z/16)^*, a = d^-1 mod 4
    GR:=[[Integers()!((Integers(4)!d)^-1),4*b,c,d]: b,c in [0..3], d in U(16)];
  when "aut4_24":    // (aP + 6bQ, cP + dQ), b,c in Z/4, d in {1,5,7,11}, a = d^-1 (1+6bc) mod 4
    GR:=[[Integers()!((Integers(4)!d)^-1*(1+6*b*c)),6*b,c,d]: b,c in [0..3], d in [1,5,7,11]];
  when "aut4_20":    // Q4=5Q, Q5=4Q: P -> aP + bQ4, Q -> cP + dQ4 - uQ5, (a b; c d) in SL2(Z/4), u in (Z/5)^*
    GR:=[[a,5*b,c,5*d-4*u]: a,b,c,d in [0..3], u in [1..4] | (a*d-b*c) mod 4 eq 1];
  when "inv6_12":    // diamond involution (P,Q) -> (5P,5Q)
    GR:=[[1,0,0,1],[5,0,0,5]];
end case;
printf "group %o for X1(%o,%o): %o transformations (before identification of -1)\n", grp, m, n, #GR;
function tatemap(E,Q) // (b,c) of (E,Q) and the map E -> Tate model
  k:=BaseRing(E); a1,a2,a3,a4,a6:=Explode(aInvariants(E)); x0:=Q[1]; y0:=Q[2];
  R<X,Y>:=PolynomialRing(k,2); f:=Y^2+a1*X*Y+a3*Y-(X^3+a2*X^2+a4*X+a6);
  g:=Evaluate(f,[X+x0,Y+y0]); A3:=MonomialCoefficient(g,Y); A4:=-MonomialCoefficient(g,X); s:=A4/A3;
  g:=Evaluate(g,[X,Y+s*X]); A1:=MonomialCoefficient(g,X*Y); A2:=-MonomialCoefficient(g,X^2); A3:=MonomialCoefficient(g,Y);
  lam:=A3/A2;
  return -A2^3/A3^2, 1-A1*A2/A3, x0, y0, s, lam;
end function;
function coords(E,P,Q)
  b,c,x0,y0,s,lam:=tatemap(E,Q);
  X:=P[1]-x0; Y:=P[2]-y0-s*X;
  return [b,c,X/lam^2,Y/lam^3];
end function;
function cdeg(v,e) for f in Divisors(e) do if forall{z: z in v | z^(q^f) eq z} then return f; end if; end for; end function;
function key(v,e) return Min([ &cat[ [Integers()!z: z in Eltseq(w^(q^i),GF(p))] : w in v] : i in [0..e-1]]); end function;
for e in [Dmin..D] do
  K:=GF(q^e); zeta:=K!1; if m gt 2 then zeta:=K!(PrimitiveElement(GF(q))^((q-1) div m)); end if;
  visited:={K|}; seen:={}; orbs:=[]; npts:=0;
  for j in K do
    if j in visited then continue; end if;
    visited join:={j^(q^i): i in [0..e-1]};
    E0:=EllipticCurveWithjInvariant(j);
    Ts:=(j eq 0 or j eq K!1728) select Twists(E0) else [E0, QuadraticTwist(E0)];
    for T in Ts do
      N:=#T; if N mod (m*n) ne 0 then continue; end if;
      G,mG:=AbelianGroup(T); inv:=Invariants(G);
      if #inv lt 2 or inv[1] mod m ne 0 then if m gt 1 then continue; end if; end if;
      Hn:=sub<G | [(inv[i] div Gcd(inv[i],n))*G.i : i in [1..#inv]]>;
      Hm:=sub<G | [(inv[i] div Gcd(inv[i],m))*G.i : i in [1..#inv]]>;
      Qs:=[mG(g): g in Hn | Order(g) eq n]; Ps:=[mG(g): g in Hm | Order(g) eq m];
      for Q in Qs do for P in Ps do
        // independence and pairing
        if m gt 2 and WeilPairing(P,(n div m)*Q,m) ne zeta then continue; end if;
        if m eq 2 and P eq (n div 2)*Q then continue; end if;
        if m gt 2 and Order(WeilPairing(P,(n div m)*Q,m)) ne m then continue; end if;
        v:=coords(T,P,Q); if cdeg(v,e) ne e then continue; end if;
        k0:=key(v,e); if k0 in seen then continue; end if;
        orb:={};
        for gg in GR do
          P2:=gg[1]*P+gg[2]*Q; Q2:=gg[3]*P+gg[4]*Q;
          assert Order(P2) eq m and Order(Q2) eq n;
          if m eq 2 then assert P2 ne (n div 2)*Q2; else assert WeilPairing(P2,(n div m)*Q2,m) eq zeta; end if;
          Include(~orb, key(coords(T,P2,Q2),e));
        end for;
        seen join:=orb; Append(~orbs,#orb);
      end for; end for;
    end for;
  end for;
  printf "ORBMN X1(%o,%o) over F_%o, degree %o, group %o: %o closed points, %o orbits, sizes %o\n", m, n, q, e, grp, #seen, #orbs, {* x: x in orbs *};
end for;
print "ORBITSMN_DONE"; quit;
