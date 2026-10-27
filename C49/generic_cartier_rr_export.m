/* Direct Y=X1(63)/<8> model and dyadic Cartier FI certificate.
   Y parametrizes the genus-zero X1(9) universal family and x(P7), with
   P7 of exact order7. Its function field is psi7(E9)(u)=0, degree24.
   All6542 degree6 cusp divisors are tested (stronger than image-only).
   Model checks: irreducibility, genus49, degree(j)=864, exact cusp signature.
   Rank JY(Q)=0 and good reduction at2 are theoretical prerequisites.
   TestN=21 runs the identical AF-only differential/jet code on known X21.
*/
SetColumns(0);
if not assigned D then D:=6; elif Type(D) eq MonStgElt then D:=StringToInteger(D); end if;
if not assigned TestN then TestN:=0; elif Type(TestN) eq MonStgElt then TestN:=StringToInteger(TestN); end if;
if not assigned Limit then Limit:=0; elif Type(Limit) eq MonStgElt then Limit:=StringToInteger(Limit); end if;
if not assigned DirectCheck then DirectCheck:=1; elif Type(DirectCheck) eq MonStgElt then DirectCheck:=StringToInteger(DirectCheck); end if;
if not assigned StopAfterExceptions then StopAfterExceptions:=1; elif Type(StopAfterExceptions) eq MonStgElt then StopAfterExceptions:=StringToInteger(StopAfterExceptions); end if;
if not assigned MemGB then MemGB:=24; elif Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
if not assigned JetFile then JetFile:=""; end if;
SetMemoryLimit(MemGB*10^9); t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
if TestN gt 0 then
    N:=TestN;
    X:=MDX1(TestN,GF(2)); C:=Curve(X); FF:=FunctionField(C);
    AF,afmap:=AlgorithmicFunctionField(FF);
    jaf:=afmap(jInvariant(EllipticCurve(X`_E)));
else
    assert TestN eq 0;
    N:=63;
    F0<t>:=FunctionField(GF(2));
    b:=t^2*(t-1)*(t^2-t+1); c:=t^2*(t-1);
    E9:=EllipticCurve([1-c,-b,-b,F0!0,F0!0]); P9:=E9![0,0];
    assert 9*P9 eq E9!0 and 3*P9 ne E9!0;
    f:=DivisionPolynomial(E9,7);
    assert Degree(f) eq 24 and GCD(f,Derivative(f)) eq 1;
    printf "Y63_EQUATION degree=24 separable=true time=%o\n",Cputime(t0);
    assert IsIrreducible(f);
    printf "Y63_IRREDUCIBLE time=%o\n",Cputime(t0);
    AF<u>:=FunctionField(f);
    assert Degree(AF) eq 24;
    jaf:=AF!jInvariant(E9);
end if;
gx:=Genus(AF);
if N eq 63 then assert gx eq 49; end if;
printf "RANK0_CARTIER_SETUP N=%o D=%o genus=%o direct_AF_model=true time=%o\n",N,D,gx,Cputime(t0);
Omega,om:=SpaceOfHolomorphicDifferentials(AF);
B:=[om(v):v in Basis(Omega)];
CM:=Matrix(GF(2),[Eltseq(Cartier(w) @@ om):w in B]);
U:=RowSpace(CM-IdentityMatrix(GF(2),gx)); udim:=Dimension(U);
printf "RANK0_CARTIER_FORMS N=%o image_dimension=%o time=%o\n",N,udim,Cputime(t0);
if udim lt D then print "RANK0_CARTIER_INSUFFICIENT_COTANGENT_DIMENSION"; quit; end if;
cs_all:=Poles(jaf);
printf "Y63_ALL_CUSPS pole_degree=%o signature=%o time=%o\n",Degree(jaf),{*Degree(c):c in cs_all*},Cputime(t0);
if N eq 63 then
    assert Degree(jaf) eq 864;
    assert {*Degree(c):c in cs_all*} eq {*1^^9,2^^3,3^^9,6^^1*};
end if;
cs:=[c:c in cs_all|Degree(c) le D]; ds:=[Degree(c):c in cs];
printf "RANK0_CARTIER_CUSPS N=%o signature=%o time=%o\n",N,{*d:d in ds*},Cputime(t0);
if N eq 63 and D eq 6 then assert gx eq 49 and {*d:d in ds*} eq {*1^^9,2^^3,3^^9,6^^1*}; end if;
function WeightedVectors(ds,d)
    if #ds eq 0 then return d eq 0 select [[]] else []; end if;
    result:=[];
    for k in [0..d div ds[1]] do
        for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do Append(~result,[k] cat v); end for;
    end for;
    return result;
end function;
// Export every local kernel through degree D, then check one actual full divisor.
// Exhaustive rank enumeration is performed independently by replay_fi_jets.py.
ci:=[i:i in [1..#ds]|ds[i] eq 1][1];
vv:=[0:i in [1..#ds]];vv[ci]:=D;vs:=[vv];
printf "CARTIER_EXPORT_CUSP_SIGNATURE N=%o D=%o degrees=%o\n",N,D,ds;
/* Exact local restriction kernels from DifferentialSpace, no series. */
if #JetFile eq 0 then JetFile:="work/degree8/cyclic_sieves/N45_cartier_kernels.m";end if;
Write(JetFile,Sprintf("// Direct differential-space restriction certificate.\nN:=%o; As:=[]; D:=%o; gx:=%o; o:=1; udim:=%o; JetFactor:=1; ds:=%o; perms:=[];\nlocaljets:=AssociativeArray(); basisjets:=AssociativeArray(); kernelspaces:=AssociativeArray();",N,D,gx,udim,ds):Overwrite:=true);
Write(JetFile,Sprintf("CartierMatrix:=Matrix(GF(2),%o,%o,%o); ImageBasis:=Matrix(GF(2),%o,%o,%o);",gx,gx,Eltseq(CM),udim,gx,Eltseq(BasisMatrix(U))));
VO:=VectorSpace(GF(2),gx);localjets:=AssociativeArray();
for i in [1..#cs] do
    maxm:=D div ds[i];cols:=[];annprev:=sub<VO|>;
    for m in [1..maxm] do
        printf "GENERIC_RR_KERNEL_START cusp=%o/%o multiplicity=%o time=%o\n",i,#cs,m,Cputime(t0);
        VV,vvmap:=DifferentialSpace(m*Divisor(cs[i]));
        W:=sub<VO|[VO!Eltseq(vvmap(b) @@ om):b in Basis(VV)]>;
        Write(JetFile,Sprintf("kernelspaces[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);",i,m,Dimension(W),gx,Eltseq(BasisMatrix(W))));
        ann:=Nullspace(Transpose(BasisMatrix(W)));assert annprev subset ann;
        extra:=[];
        for v in Basis(ann) do
            if v notin annprev then Append(~extra,v);annprev:=sub<VO|Basis(annprev) cat [v]>;end if;
        end for;
        assert #extra le ds[i];cols cat:=extra;
        while #cols lt m*ds[i] do Append(~cols,VO!0);end while;
        assert annprev eq ann;
        A:=Transpose(Matrix(GF(2),[Eltseq(v):v in cols]));
        assert Nullspace(A) eq W;
        J:=BasisMatrix(U)*A;localjets[<i,m>]:=J;
        Write(JetFile,Sprintf("localjets[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);",i,m,Nrows(J),Ncols(J),Eltseq(J)));
        printf "GENERIC_RR_KERNEL_DONE cusp=%o/%o multiplicity=%o kernel_dim=%o rank_on_U=%o time=%o\n",i,#cs,m,Dimension(W),Rank(J),Cputime(t0);
    end for;
    Write(JetFile,Sprintf("basisjets[%o]:=Matrix(GF(2),%o,%o,%o);",i,Nrows(A),Ncols(A),Eltseq(A)));
end for;
Write(JetFile,"print \"FI_JET_DATA_READY\";");

VO:=VectorSpace(GF(2),gx); bad:=[]; tested:=0;
for v in vs do
    if Limit gt 0 and tested ge Limit then break; end if;
    A:=ZeroMatrix(GF(2),udim,0);
    for i in [1..#cs] do if v[i] gt 0 then A:=HorizontalJoin(A,localjets[<i,v[i]>]); end if; end for;
    assert Ncols(A) eq D; r:=Rank(A);
    if tested lt DirectCheck or r ne D then
        Ediv:=&+[v[i]*Divisor(cs[i]):i in [1..#cs]]; assert Degree(Ediv) eq D;
        V,vm:=DifferentialSpace(Ediv);
        V_original:=sub<VO|[VO!Eltseq(vm(b) @@ om):b in Basis(V)]>;
        direct_rank:=udim-Dimension(U meet V_original); assert direct_rank eq r;
        printf "RANK0_CARTIER_DIRECT N=%o index=%o jet_rank=%o direct_rank=%o time=%o\n",N,tested+1,r,direct_rank,Cputime(t0);
    end if;
    tested+:=1;
    if r ne D then Append(~bad,v); printf "RANK0_CARTIER_EXCEPTION N=%o index=%o rank=%o vector=%o time=%o\n",N,tested,r,v,Cputime(t0); end if;
    if tested mod 1000 eq 0 then printf "RANK0_CARTIER_PROGRESS N=%o tested=%o total=%o exceptions=%o time=%o\n",N,tested,#vs,#bad,Cputime(t0); end if;
    if StopAfterExceptions gt 0 and #bad ge StopAfterExceptions then print "RANK0_CARTIER_EXCEPTION_LIMIT_REACHED"; break; end if;
end for;
printf "RANK0_CARTIER_EXPORT_DONE N=%o D=%o direct_tested=%o direct_total=%o exceptions=%o export_complete=%o time=%o\n",N,D,tested,#vs,#bad,tested eq #vs,Cputime(t0);
quit;
