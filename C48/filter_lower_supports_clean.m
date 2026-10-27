/* Complete low supports for a rational effective degree8 trace:
   singles3..5; (3,3)+cusps2; (3,4)+cusp1; (3,5); (4,4).
   First support only is reduced by diamonds, all second supports retained. */
tl48:=Cputime();ol48:=Open("C48/logs/clean_run/filter_lower_supports.log","w");
al48:=Open("C48/data/clean_run/lower_support_arithmetic.py","w");
fprintf al48,"relations = %o\nrational_cusp_generators = %o\nfinite_cusp_classes = %o\nfinite_cusp_degrees = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(a):a in gens],[Eltseq(a):a in cuspclasses],[Degree(c):c in cs];
fprintf al48,"all_source_pairs = {}\nall_source_classes = {}\nrepresentative_classes = {}\n";
allClasses48:=AssociativeArray();repClasses48:=AssociativeArray();sourceClassCache48:=AssociativeArray();
for dd48 in [3..5] do
 allClasses48[dd48]:=[];
 for i48 in [1..#sourceAllBC48[dd48]] do
  bc48:=sourceAllBC48[dd48][i48];aa48:=G!0;
  C48TargetClass(bc48,dd48,~targetClassCache48,~aa48);
  Append(~allClasses48[dd48],aa48);
  kk48:=[dd48] cat MDM2FrobeniusKey(bc48,sourceFields48[dd48],5,dd48);
  assert not IsDefined(sourceClassCache48,kk48);sourceClassCache48[kk48]:=aa48;
 end for;
 repClasses48[dd48]:=[sourceClassCache48[[dd48] cat MDM2FrobeniusKey(bc,sourceFields48[dd48],5,dd48)]:bc in sourceRepBC48[dd48]];
 fprintf al48,"all_source_pairs[%o] = %o\nall_source_classes[%o] = %o\nrepresentative_classes[%o] = %o\n",dd48,[[[Integers()!a:a in Eltseq(v,GF(5))]:v in bc]:bc in sourceAllBC48[dd48]],dd48,[Eltseq(a):a in allClasses48[dd48]],dd48,[Eltseq(a):a in repClasses48[dd48]];Flush(al48);
 fprintf ol48,"C48_LOW_CLASSES_DONE degree=%o all=%o representatives=%o time=%o\n",dd48,#allClasses48[dd48],#repClasses48[dd48],Cputime(tl48);Flush(ol48);
end for;
fprintf al48,"tests = [\n";Flush(al48);
lowTested48:=0;lowAdmitted48:=0;
procedure C48LowTest(a,lambda,indices,~tested,~admitted)
 ee:=&+lambda;assert ee le 8;
 pass:=-pi(a) in comps48[8-ee];tested+:=1;if pass then admitted+:=1;end if;
 fprintf al48,"[%o,%o,%o,%o],\n",lambda,indices,Eltseq(a),pass select 1 else 0;
 if pass then fprintf ol48,"C48_LOW_ADMITTED support=%o indices=%o class=%o time=%o\n",lambda,indices,Eltseq(a),Cputime(tl48);Flush(ol48);end if;
end procedure;
for dd48 in [3..5] do
 before48:=lowAdmitted48;
 for i48 in [1..#repClasses48[dd48]] do
  C48LowTest(repClasses48[dd48][i48],[dd48],[i48],~lowTested48,~lowAdmitted48);
 end for;
 fprintf ol48,"C48_LOW_TYPE support=%o tested=%o admitted=%o time=%o\n",[dd48],#repClasses48[dd48],lowAdmitted48-before48,Cputime(tl48);Flush(ol48);
end for;
for tp48 in [<3,3>,<3,4>,<3,5>,<4,4>] do
 d1,d2:=Explode(tp48);before48:=lowAdmitted48;nt48:=0;
 for i48 in [1..#repClasses48[d1]] do
  for j48 in [1..#allClasses48[d2]] do
   C48LowTest(repClasses48[d1][i48]+allClasses48[d2][j48],[d1,d2],[i48,j48],~lowTested48,~lowAdmitted48);nt48+:=1;
  end for;
 end for;
 assert nt48 eq #repClasses48[d1]*#allClasses48[d2];
 fprintf ol48,"C48_LOW_TYPE support=%o tested=%o admitted=%o time=%o\n",[d1,d2],nt48,lowAdmitted48-before48,Cputime(tl48);Flush(ol48);
end for;
assert lowTested48 eq &+[#repClasses48[d]:d in [3..5]]+#repClasses48[3]*(#allClasses48[3]+#allClasses48[4]+#allClasses48[5])+#repClasses48[4]*#allClasses48[4];
fprintf al48,"]\ncomplete = True\n";delete al48;
fprintf ol48,"C48_LOW_SUPPORT_FILTER_DONE tested=%o admitted=%o complete=true time=%o\n",lowTested48,lowAdmitted48,Cputime(tl48);Flush(ol48);
printf "C48_LOW_SUPPORT_FILTER_DONE tested=%o admitted=%o complete=true time=%o\n",lowTested48,lowAdmitted48,Cputime(tl48);
