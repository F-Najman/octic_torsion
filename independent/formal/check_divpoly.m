// Checks the division-polynomial convention used by XMultNumerator (formal_common.m):
// for odd n, x([n]R) = XMultNumerator(E,n)(x(R)) / f_n(x(R))^2 on random points of
// random curves over F_{2^13}.
SetColumns(0);
load "independent/formal/formal_common.m";
K := GF(2^13);
ok := true; tests := 0;
for trial in [1..20] do
  E := EllipticCurve([K| Random(K), Random(K), Random(K), Random(K), Random(K)]) where _ := 0;
  for n in [5,7] do
    phi := XMultNumerator(E, n); fn := DivisionPolynomial(E, n);
    assert Degree(phi) eq n^2;
    for k in [1..5] do
      R := Random(E); if R eq E!0 or n*R eq E!0 then continue; end if;
      nR := n*R; X := R[1]/R[3];
      tests +:= 1;
      if nR[1]/nR[3] ne Evaluate(phi, X)/Evaluate(fn, X)^2 then ok := false; end if;
    end for;
  end for;
end for;
printf "DIVPOLY_CHECK tests %o all agree %o\n", tests, ok;
quit;
