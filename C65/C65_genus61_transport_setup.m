try
 tr65:=Cputime();assert #cusps65 eq 22 and Dimension(U65) eq 48;
 s13inverse:=map<K13->K13|a:->m13(Pullback(phi13^-1,a@@m13))>;
 dg13inverse:=hom<Y65->Y65|map<K13->Y65|a:->Y65!s13inverse(a)>,b65>;
 assert &and[dg13(dg13inverse(a)) eq a:a in fieldgens65];
 baseperm65:=[];
 for P in basecusps65 do
  f,g:=TwoGenerators(P);
  hit:=[i:i in [1..#basecusps65]|Degree(basecusps65[i]) eq Degree(P) and Valuation(s13inverse(f),basecusps65[i]) gt 0 and Valuation(s13inverse(g),basecusps65[i]) gt 0];
  assert #hit eq 1;Append(~baseperm65,hit[1]);
 end for;
 assert Sort(baseperm65) eq [1..7];
 qbase65:=[];off65:=0;
 for i in [1..#basecusps65] do
  n:=Degree(basecusps65[i]) eq 1 select 3 else 4;
  for k in [1..n] do
   Append(~qbase65,i);assert Valuation(Y65!UniformizingElement(basecusps65[i]),cusps65[off65+k]) gt 0;
  end for;
  off65+:=n;
 end for;
 assert off65 eq 22;
 qwidth65:=[-Valuation(j65,P):P in cusps65];qdegree65:=[Degree(P):P in cusps65];
 qbval65:=[Valuation(b65,P):P in cusps65];
 pp5:=[];pp13:=[];
 for i in [1..22] do
  same:=[k:k in [1..22]|qbase65[k] eq qbase65[i] and qwidth65[k] eq qwidth65[i]];
  if qbval65[i] ne 0 then
   hs:=[k:k in same|qbval65[k] eq -qbval65[i]];
  elif #same eq 1 then hs:=same;
  else hs:=[k:k in same|k ne i];end if;
  assert #hs eq 1;Append(~pp5,hs[1]);
  otherbase:=baseperm65[qbase65[i]];
  hs:=[k:k in [1..22]|qbase65[k] eq otherbase and qwidth65[k] eq qwidth65[i] and qbval65[k] eq qbval65[i]];
  if #hs ne 1 then
   f,g:=TwoGenerators(cusps65[i]);
   hs:=[k:k in hs|Valuation(dg13inverse(f),cusps65[k]) gt 0 and Valuation(dg13inverse(g),cusps65[k]) gt 0];
  end if;
  assert #hs eq 1;Append(~pp13,hs[1]);
 end for;
 assert Sort(pp5) eq [1..22] and Sort(pp13) eq [1..22];
 assert &and[pp5[pp5[i]] eq i and pp5[pp13[i]] eq pp13[pp5[i]]:i in [1..22]];
 for i in [1..22] do k:=i;for r in [1..6] do k:=pp13[k];end for;assert k eq i;end for;
 allowed65:=[i:i in [1..22]|qwidth65[i] ne 5];
 steps65:=[qdegree65[i] eq 2 select 2 else 1:i in [1..22]];
 remaining65:=Seqset(allowed65);orbits65:=[];
 while #remaining65 gt 0 do
  r:=Minimum(remaining65);orb:={r};
  repeat old:=orb;orb join:={pp5[i]:i in orb} join {pp13[i]:i in orb};until orb eq old;
  assert orb subset remaining65;Append(~orbits65,Sort(Setseq(orb)));remaining65 diff:=orb;
 end while;
 assert Sort([#o:o in orbits65]) eq [2,6,12];
 BM65:=BasisMatrix(U65);assert RowSpace(BM65*M5) eq U65 and RowSpace(BM65*M13) eq U65;
 printf "C65_GENUS61_TRANSPORT_SETUP_DONE orbits=%o source_representatives=%o source_degrees=%o steps=%o time=%o\n",[#o:o in orbits65],[o[1]:o in orbits65],[qdegree65[o[1]]:o in orbits65],[steps65[o[1]]:o in orbits65],Cputime(tr65);
 print "C65_GENUS61_CUSP_PERMUTATIONS",pp5,pp13;
catch err
 print "C65_GENUS61_TRANSPORT_SETUP_FAILED";print err`Object;
end try;
