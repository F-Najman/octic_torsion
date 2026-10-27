try
 KB57:=BaseField(K);AssignNames(~KB57,["x0"]);AssignNames(~K,["y0"]);AssignNames(~AF,["u0"]);AssignNames(~E3,["v0"]);AssignNames(~CY,["z0"]);
 Write("C57/data/C57_fixed_AS.m","// Exact AS coefficient over AF.\nx0:=AF!K!BaseField(K).1;y0:=AF!K.1;u0:=AF.1;":Overwrite:=true);
 Write("C57/data/C57_fixed_AS.m",Sprintf("aAS:=%o;assert Parent(aAS) cmpeq AF;",aAS));
 Write("C57/data/C57_fixed_basis_ratios.m","// Exact holomorphic basis relative to dx0 on the fixed model.\nx0:=CY!E3!K!BaseField(K).1;y0:=CY!E3!K.1;v0:=CY!E3.1;z0:=CY.1;BasisRatios:=[];":Overwrite:=true);
 for w in B do Write("C57/data/C57_fixed_basis_ratios.m",Sprintf("Append(~BasisRatios,CY!(%o));",w/etaX));end for;
 Write("C57/data/C57_fixed_basis_ratios.m",Sprintf("SnapshotCartierMatrix:=Matrix(GF(2),25,25,%o);SnapshotU3Matrix:=Matrix(GF(2),25,25,%o);",Eltseq(CM),Eltseq(H3)));
 print "C57_MODEL_AND_HOLOMORPHIC_BASIS_SNAPSHOT_DONE";
catch err
 print "C57_SNAPSHOT_FAILED";print err`Object;
end try;
