// full torsion of the known octic point on X_1(40) (van Hoeij): 2-division field and reductions at many primes
SetColumns(0); load "independent/classification/helpers8.m";
Qz<Z>:=PolynomialRing(Rationals());
K<a>:=NumberField(Z^8-3*Z^7+3*Z^6-2*Z^5+13/5*Z^4-8/5*Z^3-7/5*Z^2+2*Z-2/5);
yy:=(5455*a^7-25010*a^6+27370*a^5-10105*a^4+24648*a^3-20135*a^2-10773*a+14311)/6653;
bK,cK:=xy2bc(a,yy); E:=TateCurve(bK,cK); P:=E![0,0]; assert Order(P) eq 40;
f2:=DivisionPolynomial(E,2); printf "2-division polynomial over K factors with degrees %o\n", [Degree(t[1]): t in Factorization(f2)];
printf "points of order 2 over K: %o; is 20P halvable over K (point Q with 2Q = P): %o\n", #[t: t in Factorization(f2) | Degree(t[1]) eq 1], IsDivisibleBy(P,2);
O:=MaximalOrder(K); g:=0; used:=[];
for q in PrimesInInterval(3,3000) do if q eq 5 or #used ge 40 then continue; end if;
  for pp in Factorization(q*O) do if pp[2] eq 1 and Min(Valuation(bK,pp[1]),Valuation(cK,pp[1])) ge 0 and Valuation(Discriminant(E),pp[1]) eq 0 then
     kk,red:=ResidueClassField(pp[1]); h:=#TateCurve(red(bK),red(cK)); g:=Gcd(g,h); Append(~used,<q,Degree(kk),h>); end if; end for;
end for;
printf "gcd of #E(k_v) over %o good unramified odd places: %o\n", #used, g;
printf "places with #E(k_v) not divisible by 80: %o\n", [u: u in used | u[3] mod 80 ne 0];
printf "TORSION40_DONE\n"; quit;
