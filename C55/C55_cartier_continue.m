SetColumns(0);SetQuitOnError(true);N:=55;D:=8;gx:=17;gy:=17;t0:=Cputime();Limit:=0;DirectCheck:=1;StopAfterExceptions:=0;JetFile:="C55/data/C55_cartier_kernels.m";
AFY:=CY;
Omega,om:=SpaceOfHolomorphicDifferentials(AFY);
B:=[om(v):v in Basis(Omega)];
CM:=Matrix(GF(2),[Eltseq(Cartier(w) @@ om):w in B]);
U:=RowSpace(CM-IdentityMatrix(GF(2),gy)); udim:=Dimension(U);
printf "RANK0_CARTIER_FORMS N=%o image_dimension=%o time=%o\n",N,udim,Cputime(t0);
if udim lt D then print "RANK0_CARTIER_INSUFFICIENT_COTANGENT_DIMENSION"; quit; end if;
cs:=[cys[i]:i in [1..#cys]|cywidths[i] in [11,55]];ds:=[Degree(p):p in cs];
assert {*d:d in ds*} eq {*1^^4,4^^1*};
function WeightedVectors(ds,d)
    if #ds eq 0 then return d eq 0 select [[]] else []; end if;
    result:=[];
    for k in [0..d div ds[1]] do
        for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do Append(~result,[k] cat v); end for;
    end for;
    return result;
end function;
vs:=WeightedVectors(ds,D);
printf "RANK0_CARTIER_COUNT N=%o total=%o time=%o\n",N,#vs,Cputime(t0);
assert #vs eq 201;
/* Exact local restriction kernels from DifferentialSpace, no series. */

Write(JetFile,Sprintf("// Direct differential-space restriction certificate.\nN:=%o; As:=[]; D:=%o; gx:=%o; o:=1; udim:=%o; JetFactor:=1; ds:=%o; perms:=[];\nlocaljets:=AssociativeArray(); basisjets:=AssociativeArray(); kernelspaces:=AssociativeArray();",N,D,gx,udim,ds));
Write(JetFile,Sprintf("CartierMatrix:=Matrix(GF(2),%o,%o,%o); ImageBasis:=Matrix(GF(2),%o,%o,%o);",gx,gx,Eltseq(CM),udim,gx,Eltseq(BasisMatrix(U))));
VO:=VectorSpace(GF(2),gx);localjets:=AssociativeArray();
for i in [1..#cs] do
    maxm:=D div ds[i];cols:=[];annprev:=sub<VO|>;
    for m in [1..maxm] do
        printf "C55_RR_KERNEL_START cusp=%o/%o multiplicity=%o time=%o\n",i,#cs,m,Cputime(t0);
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
        printf "C55_RR_KERNEL_DONE cusp=%o/%o multiplicity=%o kernel_dim=%o rank_on_U=%o time=%o\n",i,#cs,m,Dimension(W),Rank(J),Cputime(t0);
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
printf "RANK0_CARTIER_DONE N=%o D=%o tested=%o total=%o exceptions=%o complete=%o time=%o\n",N,D,tested,#vs,#bad,tested eq #vs,Cputime(t0);

