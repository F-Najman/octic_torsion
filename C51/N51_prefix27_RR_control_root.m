SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
X:=MDX1(51,GF(2));C:=Curve(X);FF:=FunctionField(C);AF,afmap:=AlgorithmicFunctionField(FF);
jaf:=afmap(jInvariant(EllipticCurve(X`_E)));gx:=Genus(AF);assert gx eq 65;
Omega,om:=SpaceOfHolomorphicDifferentials(AF);B:=[om(v):v in Basis(Omega)];
CM:=Matrix(GF(2),[Eltseq(Cartier(w) @@ om):w in B]);U:=RowSpace(CM-IdentityMatrix(GF(2),gx));assert Dimension(U) eq 62;
cs_all:=Poles(jaf);cs:=[c:c in cs_all|Degree(c) le 8];ds:=[Degree(c):c in cs];
path:="C51/data/N51_prefix27_RR_control_root_data.m";
Write(path,Sprintf("ds:=%o; CartierMatrix:=Matrix(GF(2),65,65,%o); ImageBasis:=Matrix(GF(2),62,65,%o); kernelspaces:=AssociativeArray();",ds,Eltseq(CM),Eltseq(BasisMatrix(U))):Overwrite:=true);
print "N51_PREFIX27_RR_CONTROL_SETUP",Cputime(t0);
VO:=VectorSpace(GF(2),gx);
for par in [[1,1],[1,8],[15,1]] do
 i,m:=Explode(par);VV,vvmap:=DifferentialSpace(m*Divisor(cs[i]));W:=sub<VO|[VO!Eltseq(vvmap(b) @@ om):b in Basis(VV)]>;
 Write(path,Sprintf("kernelspaces[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);",i,m,Dimension(W),gx,Eltseq(BasisMatrix(W))));
 print "N51_PREFIX27_RR_CONTROL_KERNEL",i,m,Dimension(W),Cputime(t0);
end for;
print "N51_PREFIX27_RR_CONTROLS_COMPLETE",Cputime(t0);quit;
