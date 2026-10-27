SetColumns(0);SetQuitOnError(true);
SetMemoryLimit(2*10^9);
F:=GF(256);
for n in [91] do
    qualifying:=0;
    for j in F do
        E0:=EllipticCurveFromjInvariant(j);
        for E in Twists(E0) do
            o:=#E;
            if o mod n ne 0 then continue; end if;
            qualifying+:=1;
            assert o eq 3*n and j eq 0;
            A,mp:=AbelianGroup(E);
            assert Invariants(A) eq [3*n];
            P:=3*mp(A.1);assert Order(P) eq n;
            auts:=Automorphisms(E);
            stabilizers:=[a:a in [1..n-1]|GCD(a,n) eq 1 and &or[u(P) eq a*P:u in auts]];
            print "GOOD91_STABILIZER_CHECK",stabilizers;
            assert 16 in stabilizers;
            printf "GOOD2 N=%o j=%o curve=%o order=%o aut_count=%o stabilizer_multipliers=%o\n",n,j,aInvariants(E),o,#auts,stabilizers;
        end for;
    end for;
    printf "GOOD2_ALL_CURVES N=%o qualifying_twists=%o\n",n,qualifying;
end for;
assert qualifying eq 2;
print "OCTIC_GOOD91_AUTOMORPHISM_DONE complete=true";
quit;
