/* Exact cusp-space certificate; global rank and odd-torsion lemmas separate. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
k<z>:=GF(16);X:=MDX1(25,k);C:=Curve(X);FF:=FunctionField(C);
E:=EllipticCurve(X`_E);P:=E!X`_P;
xy6:=TateNormalForm_xy(E,6*P);xy21:=TateNormalForm_xy(E,21*P);
phi:=iso<C->C|[xy6[1],xy6[2],FF!1],[xy21[1],xy21[2],FF!1]:Check:=false>;
F,mp:=AlgorithmicFunctionField(FF);assert Genus(F) eq 12;
sig:=map<F->F|a:->mp(Pullback(phi,a@@mp))>;
fieldgens:=[F!BaseField(F).1,F.1];assert &and[(sig^5)(a) eq a:a in fieldgens];
Om,om:=SpaceOfHolomorphicDifferentials(F);basis:=[om(v):v in Basis(Om)];
eta:=Differential(fieldgens[1]);sige:=Differential(sig(fieldgens[1]));
M:=Matrix(k,[Eltseq((sig(w/eta)*sige)@@om):w in basis]);I:=IdentityMatrix(k,12);
assert M^5 eq I;
cs:=[];
for w in basis do v:=w;for i in [1..4] do v:=Cartier(v);end for;Append(~cs,Eltseq(v@@om));end for;
C4:=Matrix(k,cs);assert C4*M eq M*C4;
U:=RowSpace((C4-I)*(M-I));
printf "M5N25_TARGET25_COTANGENT genus=12 diamond_image=%o invariant_dimension=%o cartier4_rank=%o joint_image=%o time=%o\n",Rank(M-I),12-Rank(M-I),Rank(C4),Dimension(U),Cputime(tt);
cusps:=[p:p in Poles(F!mp(jInvariant(E)))];
// Poles of the actual universal j-map are exactly the cusps.
assert &+[-Valuation(F!mp(jInvariant(E)),p)*Degree(p):p in cusps] eq 300;
printf "M5N25_TARGET25_CUSP_DEGREES %o time=%o\n",Sort([Degree(p):p in cusps]),Cputime(tt);
rats:=[p:p in cusps|Degree(p) eq 1];assert #rats eq 18;
V:=VectorSpace(k,12);Ws:=AssociativeArray();
for i in [1..#rats] do for m in [1,2] do
 VV,vm:=DifferentialSpace(m*Divisor(rats[i]));
 W:=sub<V|[V!Eltseq(vm(v)@@om):v in Basis(VV)]>;Ws[<i,m>]:=W;
 assert Dimension(W) eq 12-m;
end for;end for;
bad:=[];tested:=0;
for i in [1..#rats] do for j in [i..#rats] do
 W:=i eq j select Ws[<i,2>] else Ws[<i,1>] meet Ws[<j,1>];
 r:=Dimension(U)-Dimension(U meet W);tested+:=1;
 if r ne 2 then Append(~bad,[i,j,r]);end if;
end for;end for;
assert tested eq 171 and bad eq [];
function EncodeElement16(u)
 return &+[Integers()!v[i]*2^(i-1):i in [1..#v]] where v is Eltseq(k!u,GF(2));
end function;
spaces:=[*[[[EncodeElement16(a):a in Eltseq(v)]:v in Basis(Ws[<i,m>])]:m in [1,2]]:i in [1..18]*];
Write("C5x25/data/M5N25_target25_cusp_spaces_root.txt",Sprint([*Coefficients(DefiningPolynomial(k,GF(2))),spaces*]):Overwrite:=true);
printf "M5N25_TARGET25_ALLBAD_DIAGNOSTIC_DONE tested=%o exceptions=%o data=%o time=%o\n",tested,#bad,bad,Cputime(tt);quit;
