SetColumns(0);SetQuitOnError(true);SetMemoryLimit(24*10^9);tt:=Cputime();
CQ:=SmallModularCurve(11);jraw:=jInvariant(CQ,11);
CC:=ChangeRing(CQ,GF(2));assert IsNonsingular(CC) and Genus(CC) eq 1;
RC:=CoordinateRing(Ambient(CC));jn:=RC!Numerator(jraw);jd:=RC!Denominator(jraw);
K,kmap:=AlgorithmicFunctionField(FunctionField(CC));jK:=kmap(FunctionField(CC)!(jn/jd));
// The classical X0(7) parameter has j=(v^2+v+1)^4/v^7 in characteristic2.
Fv<v>:=FunctionField(GF(2));j7:=(v^2+v+1)^4/v^7;j7dual:=(v^2+v+1)^4/v;
Phi7:=ClassicalModularPolynomial(7);
assert Evaluate(Phi7,[j7,j7dual]) eq 0;
R<x>:=PolynomialRing(K);f:=x^8+jK*x^7+x^4+1;
assert IsIrreducible(f);AF<u>:=FunctionField(f);assert Genus(AF) eq 7;
jaf:=AF!jK;
cs:=Poles(jaf);widths:=[-Valuation(jaf,c):c in cs];
assert #cs eq 4 and &and[Degree(c) eq 1:c in cs] and Seqset(widths) eq {1,7,11,77};
printf "OCTIC77_X0_MODEL genus7=true modular_polynomial_verified=true cusp_widths=%o time=%o\n",widths,Cputime(tt);
load "C77/C77_AS_component.m";
RF<Y>:=PolynomialRing(AF);f2:=Y^2+Y+a2;assert IsIrreducible(f2);
// Direct cyclic cubic: the fine X1(7) family over its genus-zero parameter.
FB<t>:=FunctionField(GF(2));b7:=t^3-t^2;c7:=t^2-t;
E7u:=EllipticCurve([1-c7,-b7,-b7,FB!0,FB!0]);P7u:=E7u![0,0];assert Order(P7u) eq 7;
ju:=jInvariant(E7u);vu:=(t^2+t)/(t^3+t+1);
assert vu^8+ju*vu^7+vu^4+1 eq 0;
assert Evaluate(Numerator(ju),1/(t+1))/Evaluate(Denominator(ju),1/(t+1)) eq ju;
// j on X011 is j0(x)+j1(x)y, with y quadratic over F2(x).
jparts:=Eltseq(jK);assert #jparts eq 2 and jparts[2] ne 0;
kpoly:=DefiningPolynomial(K);assert Degree(kpoly) eq 2;
FX<X>:=FunctionField(FB);rx:=hom<BaseField(K)->FX|X>;
j0:=rx(jparts[1]);j1:=rx(jparts[2]);aa:=rx(Coefficient(kpoly,1));bb:=rx(Coefficient(kpoly,0));
normrel:=((FX!ju)+j0)^2+aa*j1*((FX!ju)+j0)+bb*j1^2;
f3:=Numerator(normrel);assert Degree(f3) eq 12 and IsIrreducible(f3);
E3<s>:=FunctionField(f3);
kx:=hom<BaseField(K)->E3|s>;
yval:=((E3!ju)+kx(jparts[1]))/kx(jparts[2]);
kemb:=hom<K->E3|kx,yval>;
afemb:=hom<AF->E3|kemb,E3!vu>;
assert afemb(jaf) eq E3!ju;
tau:=hom<E3->E3|hom<FB->E3|1/((E3!t)+1)>,s>;
assert tau(tau(tau(E3!t))) eq E3!t and tau(tau(tau(s))) eq s;
assert &and[tau(afemb(v)) eq afemb(v):v in [AF.1,AF!K.1,AF!K!BaseField(K).1]];
assert Genus(E3) eq 19;
printf "OCTIC77_DIRECT_CUBIC_VERIFIED degree12=true genus19=true exact_base_embedding=true deck_order3=true time=%o\n",Cputime(tt);
QY<z>:=PolynomialRing(E3);fy:=QY![afemb(c):c in Coefficients(f2)];assert IsIrreducible(fy);
CY<zz>:=FunctionField(fy);jY:=CY!E3!ju;
assert Genus(CY) eq 37;
printf "OCTIC77_DIRECT_GENUS37_VERIFIED time=%o\n",Cputime(tt);
cys:=Poles(jY);cywidths:=[-Valuation(jY,P):P in cys];
assert {*<cywidths[i],Degree(cys[i])>:i in [1..#cys]*} eq {*<77,1>^^6,<7,2>^^3,<11,3>^^2,<1,6>^^1*};
assert &+[cywidths[i]*Degree(cys[i]):i in [1..#cys]] eq 576;
printf "OCTIC77_GLOBAL_CUSPS_VERIFIED signature=%o width_degree_pairs=%o jdegree576=true time=%o\n",[Degree(P):P in cys],[<cywidths[i],Degree(cys[i])>:i in [1..#cys]],Cputime(tt);
print "OCTIC77_DIRECT_GLOBAL_GENUS37_MODEL_DONE";
SetColumns(0);SetQuitOnError(true);N:=77;D:=8;gx:=37;gy:=37;t0:=Cputime();Limit:=0;DirectCheck:=1;StopAfterExceptions:=0;if not assigned JetFile then JetFile:="C77/data/C77_full_RR_kernels.m";end if;
AFY:=CY;
Omega,om:=SpaceOfHolomorphicDifferentials(AFY);
B:=[om(v):v in Basis(Omega)];
CM:=Matrix(GF(2),[Eltseq(Cartier(w) @@ om):w in B]);
U:=RowSpace(CM*(CM-IdentityMatrix(GF(2),gy))); udim:=Dimension(U);
printf "RANK0_CARTIER_FORMS N=%o image_dimension=%o time=%o\n",N,udim,Cputime(t0);
if udim lt D then print "RANK0_CARTIER_INSUFFICIENT_COTANGENT_DIMENSION"; quit; end if;
cs:=[cys[i]:i in [1..#cys]|cywidths[i] eq 77];ds:=[Degree(p):p in cs];
assert {*d:d in ds*} eq {*1^^6*};
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
assert #vs eq 1287;
/* Exact local restriction kernels from DifferentialSpace, no series. */

Write(JetFile,Sprintf("// Direct differential-space restriction certificate.\nN:=%o; As:=[]; D:=%o; gx:=%o; o:=1; udim:=%o; JetFactor:=1; ds:=%o; perms:=[];\nlocaljets:=AssociativeArray(); basisjets:=AssociativeArray(); kernelspaces:=AssociativeArray();",N,D,gx,udim,ds) : Overwrite:=true);
Write(JetFile,"CartierImagePolynomial:=\"C*(C-1)\";");
Write(JetFile,Sprintf("CartierMatrix:=Matrix(GF(2),%o,%o,%o); ImageBasis:=Matrix(GF(2),%o,%o,%o);",gx,gx,Eltseq(CM),udim,gx,Eltseq(BasisMatrix(U))));
VO:=VectorSpace(GF(2),gx);localjets:=AssociativeArray();
for i in [1..#cs] do
    maxm:=D div ds[i];cols:=[];annprev:=sub<VO|>;
    for m in [1..maxm] do
        printf "C77_RR_KERNEL_START cusp=%o/%o multiplicity=%o time=%o\n",i,#cs,m,Cputime(t0);
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
        printf "C77_RR_KERNEL_DONE cusp=%o/%o multiplicity=%o kernel_dim=%o rank_on_U=%o time=%o\n",i,#cs,m,Dimension(W),Rank(J),Cputime(t0);
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

