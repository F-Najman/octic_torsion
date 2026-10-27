SetColumns(0);SetQuitOnError(true);tt:=Cputime();
J:=JH(25,[1]);assert Dimension(J) eq 12;
for p in [2,11,31,41,61] do
 f:=FrobeniusPolynomial(J,p);
 if p eq 2 then A:=CompanionMatrix(f);n:=Determinant(IdentityMatrix(Rationals(),24)-A^4);
 else n:=Evaluate(f,1);end if;
 printf "J25_SPLIT_REDUCTION p=%o order=%o odd=%o frobenius=%o time=%o\n",p,Integers()!n,IsOdd(Integers()!n),Coefficients(f),Cputime(tt);
end for;quit;
