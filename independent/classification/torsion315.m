// full torsion of the C3 x C15 curve of Prop. 3.84 over the octic field K = F(w), F = Q(zeta_3)
SetColumns(0);
QQz<Zz>:=PolynomialRing(Rationals()); Fz<zz>:=NumberField(Zz^2+Zz+1); Fw<W>:=PolynomialRing(Fz);
fw:=W^4-5*W^3+8*W^2-6*W+3; printf "f irreducible over F: %o\n", IsIrreducible(fw);
Krel<w>:=NumberField(fw);
a1:=((-7*zz-8)*w^3+(23*zz+34)*w^2+(-8*zz-43)*w)/6 + (-3*zz+10)/2;
a3:=((17*zz+13)*w^3+(-82*zz-71)*w^2+(106*zz+119)*w)/12 - 2*zz - 5;
xQ:=((-zz+1)*w^3+(8*zz-2)*w^2+(-17*zz-1)*w)/6 + 3*zz/2;
yQ:=((zz+2)*w^3+(-4*zz-8)*w^2+(3*zz+6)*w+2*zz+3)/4;
K:=AbsoluteField(Krel); EK:=EllipticCurve([K!a1,0,K!a3,0,0]); P:=EK![0,0]; Q:=EK![K!xQ,K!yQ];
printf "orders %o %o, independent (5Q not in <P>): %o, j = %o\n", Order(P), Order(Q), 5*Q notin {EK!0,P,2*P}, jInvariant(EK);
OK:=MaximalOrder(K); g:=0; used:=[];
for pq in PrimesInInterval(5,500) do if #used ge 30 then break; end if;
  for pp in Factorization(pq*OK) do if pp[2] eq 1 and Valuation(K!a1,pp[1]) ge 0 and Valuation(K!a3,pp[1]) ge 0 and Valuation(Discriminant(EK),pp[1]) eq 0 then
     kk,red:=ResidueClassField(pp[1]); Eb:=EllipticCurve([red(K!a1),0,red(K!a3),0,0]); h:=#Eb; g:=Gcd(g,h);
     Append(~used,<pq,Degree(kk),h,Invariants(AbelianGroup(Eb))>); end if; end for; end for;
printf "#E(k_v) at %o good unramified odd places: %o\n", #used, used;
printf "gcd = %o (so E(K)_tors has order dividing it; it contains C3 x C15 of order 45)\n", g;
printf "TORSION315_DONE\n"; quit;
