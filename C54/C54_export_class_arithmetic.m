assert #candidate2 eq 531 and Kernel(h2) eq Kernel(h5);
data54out:="C54/data/C54_class_arithmetic_base.py";
Write(data54out,Sprintf("relations2 = %o\nrelations5 = %o\ncusp2 = %o\ncusp5 = %o\ncandidate2 = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G2))],[Eltseq(r):r in Rows(RelationMatrix(G5))],[Eltseq(v):v in qc2],[Eltseq(v):v in qc5],[Eltseq(v):v in candidate2]));
rawtests54:=Split(Read("C54/data/C54_lower_classes.txt"),"\n");fulltests54:=[* *];
for line in rawtests54 do
 if #line eq 0 then continue;end if;
 row:=eval line;assert #row eq 6;q:=G5!row[1..4];member:=q in H5;assert (member select 1 else 0) eq row[5];im2:=[];preimage:=[];passes:=false;
 if member then
  v:=q@@h5;assert h5(v) eq q;im2:=Eltseq(h2(v));preimage:=Eltseq(v);passes:=h2(v) in candidate2;
 end if;
 assert (passes select 1 else 0) eq row[6];Append(~fulltests54,[*row[1..4],row[5],row[6],im2,preimage*]);
end for;
classout54:=Open("C54/data/C54_replayed_testclasses.py","w");fprintf classout54,"tests = [\n";
for row in fulltests54 do fprintf classout54,"[%o,%o,%o,%o,%o],\n",row[1],row[2],row[3],row[4],row[5];end for;
fprintf classout54,"]\ncomplete = True\n";delete classout54;
printf "OCTIC54_CLASS_ARITHMETIC_EXPORTED tests=%o members=%o\n",#fulltests54,#[r:r in fulltests54|r[2] eq 1];
