assert Dimension(Omega65) eq 61 and Dimension(U65) eq 48;
data65:="C65/data/C65_root_allbad_matrices.m";
Write(data65,Sprintf("CartierMatrix:=Matrix(GF(2),61,61,%o);",Eltseq(CM65)));
Write(data65,Sprintf("Diamond5Matrix:=Matrix(GF(2),61,61,%o);",Eltseq(M5)));
Write(data65,Sprintf("Diamond13Matrix:=Matrix(GF(2),61,61,%o);",Eltseq(M13)));
Write(data65,Sprintf("ImageBasis:=Matrix(GF(2),48,61,%o);",Eltseq(BasisMatrix(U65))));
print "C65_GENUS61_RETAINED_MATRICES_EXPORTED";
