// misc8.m -- small independent checks of explicit statements in the draft.  Marker MISC, MISC8_DONE.
SetColumns(0); SetMemoryLimit(8*10^9);
Q:=Rationals();
// (1) Lemma 3.95: E0: y^2 = x^3 - 7x - 6, the map to Z, ranks over Q and of the (-1)-twist, torsion over Q(i)
E0:=EllipticCurve([0,0,0,-7,-6]);
r1:=Rank(E0); r2:=Rank(QuadraticTwist(E0,-1));
printf "MISC L3.95: rank E0(Q) = %o (proved bound %o), rank E0^(-1)(Q) = %o (bound %o)\n", r1, RankBound(E0), r2, RankBound(QuadraticTwist(E0,-1));
K<i>:=QuadraticField(-1); EK:=ChangeRing(E0,K); T,mT:=TorsionSubgroup(EK);
printf "MISC L3.95: E0(Q(i))_tors = %o ; points %o\n", Invariants(T), [mT(g): g in T];
E9:=ChangeRing(E0,GF(9)); printf "MISC L3.95: E0(F9) = %o\n", Invariants(AbelianGroup(E9));
FF<X,Y>:=FunctionField(EK);
t:=(X-3)^2*(X+2)/(X+1)^2; w:=Y*(X+1)^2/(4*(X+2));
printf "MISC L3.95: map identity 16(w^2+16)^3 t = (t^2+10t+5)^3 w^4 on E0: %o\n", 16*(w^2+16)^3*t eq (t^2+10*t+5)^3*w^4;
for P in [EK![-1+2*i,4-2*i], EK![-1+2*i,-(4-2*i)], EK![-1-2*i,4+2*i], EK![-1-2*i,-(4+2*i)]] do
  tv:=(P[1]-3)^2*(P[1]+2)/(P[1]+1)^2; wv:=P[2]*(P[1]+1)^2/(4*(P[1]+2));
  printf "MISC L3.95: point %o -> t=%o, w=%o, j=%o\n", P, tv, wv, 16*(wv^2+16)^3/wv^4;
end for;
// (2) Prop 3.41: J0(49) model
E49:=EllipticCurve([1,-1,0,-2,-1]); T49:=E49![2,-1];
printf "MISC P3.41: conductor %o, disc %o, Order(T)=%o, reduction of T mod 2: %o, E49(Q)_tors=%o\n", Conductor(E49), Discriminant(E49), Order(T49),
   ChangeRing(E49,GF(2))![2,-1], Invariants(TorsionSubgroup(E49));
// (3) Lemma 3.48: the X(10) / X(5) data
Fz<z>:=CyclotomicField(5); Ft<tt>:=FunctionField(Fz);
a:=tt^20-228*tt^15+494*tt^10+228*tt^5+1; b:=tt^30+522*tt^25-10005*tt^20-10005*tt^10-522*tt^5+1; h:=tt^5*(tt^10+11*tt^5-1)^5;
printf "MISC L3.48: b^2-a^3 = 1728h: %o\n", b^2-a^3 eq 1728*h;
Et:=EllipticCurve([Ft| -27*a, 54*b]);
printf "MISC L3.48: j(E) = -a^3/h: %o ; deg j = %o\n", jInvariant(Et) eq -a^3/h, Degree(Numerator(-a^3/h));
R<u>:=PolynomialRing(Ft); f5:=R!DivisionPolynomial(Et,5); fac:=Factorization(f5);
printf "MISC L3.48: 5-division polynomial degree %o factors as degrees %o\n", Degree(f5), [Degree(x[1]): x in fac];
// (4) Prop 3.84: twist class of the C3 x C15 curve relative to LMFDB 338.d2 = [1,1,1,3,-5]
Fz3<zz>:=CyclotomicField(3); Fu<uu>:=PolynomialRing(Fz3); Kr<w>:=NumberField(uu^4-5*uu^3+8*uu^2-6*uu+3); KK:=AbsoluteField(Kr); w:=KK!w; z:=KK!(Kr!zz);
a1:=((-7*z-8)*w^3+(23*z+34)*w^2+(-8*z-43)*w)/6 + (-3*z+10)/2;
a3:=((17*z+13)*w^3+(-82*z-71)*w^2+(106*z+119)*w)/12 - 2*z - 5;
E:=EllipticCurve([a1,0,a3,0,0]); E338:=EllipticCurve([1,1,1,3,-5]);
printf "MISC P3.84: j(338.d2) = %o, conductor %o\n", jInvariant(E338), Conductor(E338);
c4,c6:=Explode(cInvariants(E)); d4,d6:=Explode(cInvariants(E338)); dd:=(c6*d4)/(d6*c4);
printf "MISC P3.84: E = twist of 338.d2 by d: %o ; Norm_{K/Q}(d) = %o, a square in Q: %o\n", IsIsomorphic(QuadraticTwist(ChangeRing(E338,KK),dd),E), Norm(dd), IsSquare(Norm(dd));
OK:=MaximalOrder(KK); fd:=Factorization(dd*OK); printf "MISC P3.84: support of (d): %o\n", [<Norm(x[1]), x[2]>: x in fd];
// is E a base change from F = Q(zeta3) of a twist of 338.d2?  d in F^* K^*2  <=>  N_{K/F}(d) ... test directly on small F-elements
Fs:=[x+y*z : x,y in [-12..12] | x ne 0 or y ne 0];
hit:=[s: s in Fs | IsSquare(dd*s)];
printf "MISC P3.84: elements s = x+y*zeta3, |x|,|y|<=12, with d*s a square in K: %o\n", hit;
// (5) Prop 3.32: automorphism multipliers on the j=0 curves over F256 with 273 points
F256:=GF(2^8); Ts:=Twists(EllipticCurveWithjInvariant(F256!0));
for T in [T: T in Ts | #T eq 273] do
  G,mG:=AbelianGroup(T); P:=mG(G.1); P:=(#G div 91)*P; A:=AutomorphismGroup(T);
  mult:={}; for u in Automorphisms(T) do Pu:=u(P); for k in [1..90] do if k*P eq Pu then Include(~mult,k); end if; end for; end for;
  printf "MISC P3.32: curve %o, group %o, #Aut %o, multiplier set on a point of order 91: %o\n", aInvariants(T), Invariants(G), #A, mult;
end for;
print "MISC8_DONE"; quit;
