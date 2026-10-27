// placecounts_ms.m -- number of places of degree d <= 8 of X_1(N)_{F_p}, from the Eichler--Shimura Frobenius
// (modular symbols), independent of the function-field enumeration in classify8.m.
// Frobenius on H^1 restricted to S_2(N,eps): char poly det(X^2 - T_p X + p eps(p)); #X(F_{p^n}) = p^n + 1 - sum alpha^n.
SetColumns(0); load "independent/classification/helpers8.m";
function placecounts(N, p, maxd)
  G:=DirichletGroup(N, CyclotomicField(EulerPhi(N))); s:=[0: n in [1..maxd]];
  for eps0 in GaloisConjugacyRepresentatives(G) do
    if Evaluate(eps0,-1) ne 1 then continue; end if; eps:=MinimalBaseRingCharacter(eps0);
    M:=CuspidalSubspace(ModularSymbols(eps,2,+1)); if Dimension(M) eq 0 then continue; end if;
    Tp:=HeckeOperator(M,p); K:=BaseRing(Tp); n:=Nrows(Tp); e:=K!Evaluate(eps,p);
    Z:=ZeroMatrix(K,n,n); I:=IdentityMatrix(K,n);
    Fr:=BlockMatrix(2,2,[Tp, -p*e*I, I, Z]); A:=Fr;
    for k in [1..maxd] do tr:=Trace(A); if Type(tr) ne FldRatElt then tr:=Trace(tr); end if; s[k]+:=Integers()!tr; A:=A*Fr; end for;   // Trace_{K/Q}
  end for;
  Npts:=[p^k+1-s[k]: k in [1..maxd]];
  P:=[ (&+[MoebiusMu(k div d)*Npts[d]: d in Divisors(k)]) div k : k in [1..maxd]];
  return P, Npts;
end function;
for pr in [<33,2>,<34,3>,<35,2>,<39,2>,<40,3>] do
  P,Np:=placecounts(pr[1],pr[2],8);
  printf "MSPLACES N=%o p=%o: #X(F_{p^n}) n=1..8: %o; places of degree 1..8: %o\n", pr[1], pr[2], Np, P;
end for;
printf "PLACECOUNTS_DONE\n"; quit;
