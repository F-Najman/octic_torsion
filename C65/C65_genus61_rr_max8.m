try
 rrmax65:=Cputime();VO65:=VectorSpace(GF(2),61);
 ratrep65:=[o[1]:o in orbits65|qdegree65[o[1]] eq 1][1];
 printf "C65_GENUS61_RR_MAX8_START cusp=%o time=%o\n",ratrep65,Cputime(rrmax65);
 VV65,vmap65:=DifferentialSpace(8*Divisor(cusps65[ratrep65]));
 Wmax65:=sub<VO65|[VO65!Eltseq(vmap65(v)@@om65):v in Basis(VV65)]>;
 Amax65:=Transpose(BasisMatrix(Nullspace(Transpose(BasisMatrix(Wmax65)))));
 assert Nullspace(Amax65) eq Wmax65;
 sourceW65:=AssociativeArray();sourceW65[<ratrep65,8>]:=Wmax65;
 Write("C65/data/C65_root_allbad_max8.m",Sprintf("SourceCusp:=%o;SourceKernel:=Matrix(GF(2),%o,61,%o);RestrictionMatrix:=Matrix(GF(2),61,%o,%o);",ratrep65,Dimension(Wmax65),Eltseq(BasisMatrix(Wmax65)),Ncols(Amax65),Eltseq(Amax65)));
 printf "C65_GENUS61_RR_MAX8_DONE kernel_dim=%o full_rank=%o rank_on_U=%o time=%o\n",Dimension(Wmax65),Rank(Amax65),Rank(BM65*Amax65),Cputime(rrmax65);
catch err
 print "C65_GENUS61_RR_MAX8_FAILED";print err`Object;
end try;
