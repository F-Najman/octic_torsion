// jacorders8.m -- independent orders #J_H(N)(F_p) = det(1 - T_p + p<p>) on the sign +1 cuspidal modular
// symbols of Gamma_H(N) (Eichler-Shimura), for the finite-field Jacobian orders quoted in the draft.
// X1(m,n) over Q(zeta_m) at a split prime p = X_H(mn) at p with H = {a = +-1 mod n}.  Markers JAC, JACORDERS_DONE.
SetColumns(0); SetMemoryLimit(16*10^9);
function jacorder(N, Hgens, p)
  M:=ModularSymbolsH(N, Hgens, 2, +1); S:=CuspidalSubspace(M);
  T:=HeckeOperator(S,p); T2:=HeckeOperator(S,p^2);   // p<p> = T_p^2 - T_{p^2} for p not dividing N
  return Determinant(IdentityMatrix(Rationals(),Dimension(S)) - T + (T^2 - T2)), Dimension(S);
end function;
function Hmn(m,n) L:=m*n; return [a: a in [2..L-1] | Gcd(a,L) eq 1 and (a mod n eq 1 or a mod n eq n-1)]; end function;
function H1(N) return [N-1]; end function;
// validation: X1(13) is the genus-2 curve y^2 = x^6-2x^5+x^4-2x^3+6x^2-4x+1
P<x>:=PolynomialRing(GF(3)); C:=HyperellipticCurve(x^6-2*x^5+x^4-2*x^3+6*x^2-4*x+1);
o,g:=jacorder(13,H1(13),3); printf "JAC validation X1(13) at 3: modular symbols %o (genus %o), hyperelliptic model %o\n", o, g, #Jacobian(C);
o,g:=jacorder(11,[2],2); printf "JAC validation X0(11) at 2: %o (expected 5)\n", o;
claims:=[
 <"J1(23)(F3)", 23, H1(23), 3, 9406793>,
 <"J1(25)(F3)", 25, H1(25), 3, 2503105>,
 <"J1(25)(F11)", 25, H1(25), 11, 5014032203125>,
 <"J1(27)(F2)", 27, H1(27), 2, 3*3*52497>,
 <"J1(27)(F5)", 27, H1(27), 5, 3*3*342175446>,
 <"J1(28)(F3)", 28, H1(28), 3, 4*4*24*936>,
 <"J1(33)(F2)", 33, H1(33), 2, 15*42373650>,
 <"J1(34)(F3)", 34, H1(34), 3, 4^3*8760*1191360>,
 <"J1(35)(F2)", 35, H1(35), 2, 18446099880>,
 <"J1(40)(F3)", 40, H1(40), 3, 31579088486400>,
 <"J1(21)(F5)", 21, H1(21), 5, 2184>,
 <"J1(2,14)(F5)", 28, Hmn(2,14), 5, 2*2*12*36>,
 <"J1(2,24)(F5)", 48, Hmn(2,24), 5, 2^2*4^2*12*24*480^2>,
 <"J1(2,26)(F3)", 52, Hmn(2,26), 3, 124501550040000>,
 <"J1(3,15)(F7)", 45, Hmn(3,15), 7, 149299200>,
 <"J1(3,15)(F13)", 45, Hmn(3,15), 13, 12230590464>,
 <"J1(3,15)(F19)", 45, Hmn(3,15), 19, 247669456896>,
 <"J1(3,18)(F7)", 54, Hmn(3,18), 7, 3^3*9^2*63^2*189>,
 <"J1(4,16)(F5)", 64, Hmn(4,16), 5, 2^6*4^2*40^2*80^2>,
 <"J1(6,12)(F7)", 72, Hmn(6,12), 7, 6*12*24^2*48^2>,
 <"J(X1(57)/C9)(F2)", 57, [25,56], 2, 105*210>,
 <"J(X1(39)/<16>)(F2)", 39, [16,38], 2, 6*336>,
 <"J0(57)(F2)", 57, [a: a in [2..56] | Gcd(a,57) eq 1], 2, 15*30>,
 <"J0(55)(F2)", 55, [a: a in [2..54] | Gcd(a,55) eq 1], 2, 5*20>,
 <"J0(77)(F2)", 77, [a: a in [2..76] | Gcd(a,77) eq 1], 2, 30*60>,
 <"J0(69)(F2)", 69, [a: a in [2..68] | Gcd(a,69) eq 1], 2, 22*44>,
 <"J0(39)(F2)", 39, [a: a in [2..38] | Gcd(a,39) eq 1], 2, 28>,
 <"J(X(10))(F11)", 100, [11,99], 11, 2*6*30^4*90^4>,
 <"J(X(10))(F31)", 100, [11,99], 31, 2^9*6*30^3*90*270^4>

];
for c in claims do
  t0:=Cputime(); o,g:=jacorder(c[2],c[3],c[4]);
  printf "JAC %o: genus %o, det(1-T_p+p<p>) = %o = %o ; claimed %o -> %o (%o s)\n", c[1], g, o, Factorization(Integers()!o), c[5],
     c[5] eq 0 select "(no order claimed; see text)" else (o eq c[5] select "CONFIRMED" else "CONTRADICTED"), Cputime(t0);
end for;
print "JACORDERS_DONE"; quit;
