SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt6:=Cputime();
records6:=Open("C65/data/C65_good6_D8_records.txt","w");delete records6; // start a fresh records file (Write appends)
F0<z>:=FunctionField(GF(2));j0:=(z^16+z^4+1)/(z^28+z^25+z^24+z^21+z^20+z^17+z^16+z^13);
R<T>:=PolynomialRing(F0);AF<t>:=FunctionField((T+1)^12+j0*(T^7+T^6+T^5));
assert Genus(AF) eq 21;
pinf:=Poles(z)[1];pfour:=Zeros(z^4+z+1)[1];
ss21:=[];
for p in [pinf,pfour] do
 above:=Decomposition(AF,p);u:=UniformizingElement(p);
 assert &+[Valuation(AF!u,q)*Degree(q):q in above] eq 12*Degree(p);
 assert &and[Valuation(AF!j0,q) gt 0:q in above];ss21 cat:=above;
 printf "C65_GOOD6_SS21_FIBRE base_degree=%o places=%o time=%o\n",Degree(p),[<Degree(q),Valuation(AF!u,q)>:q in above],Cputime(tt6);
end for;
assert #Seqset(ss21) eq #ss21;
print "C65_GOOD6_SS21_PREFLIGHT_DONE",[Degree(p):p in ss21],Cputime(tt6);

// Retained genus21 supersingular model from C65_good6_genus21_preflight_v3.m.
try
 st6:=Cputime();
 sigT:=hom<AF->AF|map<F0->AF|a:->AF!a>,1/t>;
 sigU:=hom<AF->AF|hom<F0->AF|AF!(z+1)>,t>;
 assert sigT(sigU(t)) eq sigU(sigT(t)) and sigT(sigU(AF!z)) eq sigU(sigT(AF!z));
 assert sigT(AF!j0) eq AF!j0 and sigU(AF!j0) eq AF!j0;
 perms6:=[];
 for sig in [sigT,sigU] do
  perm:=[];
  for P in ss21 do
   f,g:=TwoGenerators(P);
   hits:=[i:i in [1..#ss21]|Valuation(sig(f),ss21[i]) gt 0 and Valuation(sig(g),ss21[i]) gt 0];
   assert #hits eq 1 and Degree(ss21[hits[1]]) eq Degree(P);Append(~perm,hits[1]);
  end for;
  assert Sort(perm) eq [1..#ss21];Append(~perms6,perm);
 end for;
 pT6:=perms6[1];pU6:=perms6[2];
 assert &and[pT6[pT6[i]] eq i and pU6[pU6[i]] eq i and pT6[pU6[i]] eq pU6[pT6[i]]:i in [1..#ss21]];
 inputs6:=[i:i in [1..#ss21]|6 mod Degree(ss21[i]) eq 0];
 targets6:=[i:i in [1..#ss21]|12 mod Degree(ss21[i]) eq 0 and Valuation(AF!z,ss21[i]) ge 0];
 assert #inputs6 eq 3 and #targets6 eq 4 and &and[pT6[pU6[i]] eq i:i in [1..#ss21]];
 printf "OCTIC65_GOOD6_SS_COVERAGE input_indices=%o target_indices=%o degrees=%o Tperm=%o Uperm=%o\n",inputs6,targets6,[Degree(p):p in ss21],pT6,pU6;
 total6:=0;sTests6:=0;uTests6:=0;hTests6:=0;bad6:=[];
 for inp in inputs6 do
  P6:=(6 div Degree(ss21[inp]))*Divisor(ss21[inp]);
  sP6:=(6 div Degree(ss21[inp]))*Divisor(ss21[pT6[inp]]);
  uP6:=(6 div Degree(ss21[inp]))*Divisor(ss21[pU6[inp]]);
  usP6:=(6 div Degree(ss21[inp]))*Divisor(ss21[pU6[pT6[inp]]]);
  for ii in [1..#targets6] do for jj in [ii..#targets6] do
   a:=targets6[ii];b:=targets6[jj];
   T6:=(12 div Degree(ss21[a]))*Divisor(ss21[a])+(12 div Degree(ss21[b]))*Divisor(ss21[b]);
   sT6:=(12 div Degree(ss21[a]))*Divisor(ss21[pT6[a]])+(12 div Degree(ss21[b]))*Divisor(ss21[pT6[b]]);
   uT6:=(12 div Degree(ss21[a]))*Divisor(ss21[pU6[a]])+(12 div Degree(ss21[b]))*Divisor(ss21[pU6[b]]);
   As6:=sT6-T6+2*(sP6-P6);assert Degree(As6) eq 0;
   ps:=IsPrincipal(As6);sTests6+:=1;pu:=false;
   if ps then
    Au6:=uT6-T6-uP6+P6-3*usP6+3*sP6;assert Degree(Au6) eq 0;
    pu:=IsPrincipal(Au6);uTests6+:=1;
   end if;
   ph:=false;
   if ps and pu then
    A2six:=4*(T6-P6-3*sP6);assert Degree(A2six) eq 0;ph:=IsPrincipal(A2six);hTests6+:=1;
   end if;
   Write("C65/data/C65_good6_D8_records.txt",Sprint([inp,a,b,ps select 1 else 0,pu select 1 else 0,ph select 1 else 0]));
   if ps and pu and ph then Append(~bad6,[inp,a,b]);end if;
   total6+:=1;
  end for;end for;
  printf "OCTIC65_GOOD6_SS_INPUT_DONE input=%o cases=%o survivors=%o time=%o\n",inp,total6,#bad6,Cputime(st6);
 end for;
 printf "OCTIC65_GOOD6_SS_DONE cases=%o s_tests=%o u_tests=%o hecke2_tests=%o survivors=%o data=%o complete=true time=%o\n",total6,sTests6,uTests6,hTests6,#bad6,bad6,Cputime(st6);
catch err
 print "OCTIC65_GOOD6_SS_FAILED";print err`Object;
end try;

quit;
