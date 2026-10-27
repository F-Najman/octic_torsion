/* Fresh full10 relative degree2 diagonal test over Q(zeta5).
   All36 global cusps matched by integral-unit valuations; only necessary
   multiples are used, with no unproved cuspidal generation assertion. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
F<zeta>:=CyclotomicField(5);RF<T>:=PolynomialRing(F);
rr:=Roots(T^10+11*T^5-1);assert #rr eq 10;
global_roots:=[r[1]:r in rr];
Gs:=[* *];Hs:=[* *];hs:=[* *];phis:=[* *];Ks:=[* *];c0s:=[* *];css:=[* *];
sigs0:=[];Free:=FreeAbelianGroup(36);
for p in [11,31] do
 k:=GF(p);zp:=PrimitiveElement(k)^((p-1) div 5);
 function Red(a)
  aa:=Eltseq(a);return &+[k!aa[i]*zp^(i-1):i in [1..#aa]];
 end function;
 A<t>:=FunctionField(k);P<x>:=PolynomialRing(A);
 aa:=t^20-228*t^15+494*t^10+228*t^5+1;
 bb:=t^30+522*t^25-10005*t^20-10005*t^10-522*t^5+1;
 den:=t^5*(t^10+11*t^5-1)^5;assert bb^2-aa^3 eq 1728*den;
 j:=-aa^3/den;f:=x^3-27*aa*x+54*bb;assert IsIrreducible(f);
 B<r>:=ext<A|f>;PB<y>:=PolynomialRing(B);ff:=y^2-(-3*r^2+108*(B!aa));
 assert IsIrreducible(ff);K<s>:=ext<B|ff>;
 assert Genus(K) eq 13 and Degree(K,A) eq 6 and Degree(K!j) eq 360;
 cs:=Poles(K!j);assert #cs eq 36 and &and[Degree(c) eq 1:c in cs];
 r1:=K!r;r2:=(-r1+s)/2;r3:=(-r1-s)/2;
 assert r1+r2+r3 eq 0 and r1*r2+r2*r3+r3*r1 eq -27*(K!aa);
 fs:=[K!t] cat [K!(t-Red(a)):a in global_roots] cat [r1-r2,r1-r3,r2-r3];
 groups:=AssociativeArray();
 for c in cs do
  sig:=[Valuation(f,c):f in fs];
  if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 36 and &and[#groups[s] eq 1:s in sigs];
 cs:=[groups[s][1]:s in sigs];
 if p eq 11 then sigs0:=sigs;else assert sigs eq sigs0;end if;
 c0:=cs[1];G,gm,phi:=ClassGroup(K:ReductionDivisor:=Divisor(c0),Proof:=true);
 gs:=[phi(Divisor(c)-Divisor(c0)):c in cs];h:=hom<Free->G|gs>;H:=Image(h);
 print "FULL10_D2_DIAGONAL_GROUP",p,Invariants(G),Invariants(H),Invariants(quo<G|H>),Cputime(tt);
 Append(~Gs,G);Append(~Hs,H);Append(~hs,h);Append(~phis,phi);Append(~Ks,K);Append(~c0s,c0);Append(~css,cs);
end for;
for a in Generators(Kernel(hs[1])) do assert hs[2](Free!a) eq Gs[2]!0;end for;
for a in Generators(Kernel(hs[2])) do assert hs[1](Free!a) eq Gs[1]!0;end for;
G11:=Gs[1];G31:=Gs[2];H11:=Hs[1];h11:=hs[1];h31:=hs[2];
Q,q:=quo<G11|H11>;assert Invariants(Q) eq [3,3,3,3,0];
predictions:=AssociativeArray();source:=[* *];target:=[* *];matches:=[* *];
for idx in [1,2] do
 K:=Ks[idx];cs:=css[idx];c0:=c0s[idx];phi:=phis[idx];
 assert #[c:c in Places(K,1)|not c in cs] eq 0;
 ps:=[c:c in Places(K,2)|not c in cs];
 print "FULL10_D2_DIAGONAL_PLACES",idx,#ps,Cputime(tt);
 for i in [1..#ps] do
  cl:=phi(Divisor(ps[i])-2*Divisor(c0));
  if idx eq 1 then
   assert Order(q(cl)) eq 3;a:=h31(Free!((3*cl)@@h11));key:=Eltseq(a);
   if IsDefined(predictions,key) then Append(~predictions[key],i);else predictions[key]:=[i];end if;
   Append(~source,[*i,Eltseq(cl),key*]);
  else
   key:=Eltseq(3*cl);hit:=IsDefined(predictions,key);Append(~target,[*i,Eltseq(cl),hit*]);
   if hit then Append(~matches,[*i,predictions[key]*]);end if;
  end if;
 end for;
end for;
Write("C10x10/data/diagonal_full10_D2_root_arithmetic.txt",Sprint([*sigs0,
 [Invariants(G):G in Gs],[[Eltseq(v):v in Rows(RelationMatrix(G))]:G in Gs],
 [[Eltseq(h(Free.i)):i in [1..36]]:h in hs],source,target,matches*]) : Overwrite:=true);
printf "FULL10_D2_DIAGONAL_DONE source=%o target=%o matches=%o complete=true time=%o\n",#source,#target,#matches,Cputime(tt);quit;
