/* Clean optional native replay of ALL210 support cases, not a continuation.
   Set C420FullStem to a NEW path before loading. Approximately 800CPU sec. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
assert assigned C420FullStem;
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_full_aut_hash_m4_16.m");
Attach("common/fast_full_aut_hash_m4_16_relative.m");
Attach("C4x20/fast_full_aut_hash_m4_20.m");
X:=MDX11(4,20,GF(9));assert Genus(X) eq 25;
assert #NoncuspidalPlaces(X,1) eq 0;
assert [#Cusps(X,d):d in [1..4]] eq [24,12,0,0];
prior:=eval Read("C4x20/data/preflight_4_20_D4_Aut96_arithmetic.txt");
allxs:=[* *];allrs:=[* *];tests:=[* *];remaining:=[* *];
for d in [2..4] do
 xs:=NoncuspidalPlaces(X,d);rs,labels,keys:=PlacesUpToFullAutM4Twenty(X,xs);
 assert #xs eq [96,320,888][d-1] and #rs eq [2,4,12][d-1];
 assert keys eq prior[d-1][2] and labels eq prior[d-1][3];
 Append(~allxs,xs);Append(~allrs,rs);
 for i in [1..#rs] do
  killed:=false;
  for q in (d le 3 select [29] else [13,17]) do
   A:=HeckeOperatorFast(X,q,rs[i])-Divisor(rs[i])-q*Divisor(DiamondOperator(X,q,rs[i]));
   pr:=IsPrincipal(A);Append(~tests,[*1,d,i,0,q,pr select 1 else 0*]);
   print "C420_FRESH_SINGLE",d,i,q,pr,Cputime(tt);
   if not pr then killed:=true;break;end if;
  end for;
  if not killed then Append(~remaining,[*1,d,i,0*]);end if;
 end for;
end for;
cache:=AssociativeArray();xs:=allxs[1];rs:=allrs[1];
for i in [1..2] do
 a:=Index(xs,rs[i]);assert a eq prior[1][5][i];
 for j in [1..96] do
  killed:=false;
  for q in [13,17] do
   for k in Seqset([a,j]) do
    ky:=[q,k];
    if not IsDefined(cache,ky) then
     cache[ky]:=HeckeOperatorFast(X,q,xs[k])-Divisor(xs[k])-q*Divisor(DiamondOperator(X,q,xs[k]));
    end if;
   end for;
   pr:=IsPrincipal(cache[[q,a]]+cache[[q,j]]);
   Append(~tests,[*2,2,i,j,q,pr select 1 else 0*]);
   print "C420_FRESH_PAIR",i,j,q,pr,Cputime(tt);
   if not pr then killed:=true;break;end if;
  end for;
  if not killed then Append(~remaining,[*2,2,i,j*]);end if;
  Write(C420FullStem cat "_arithmetic.txt",Sprint([*tests,remaining*]):Overwrite:=true);
 end for;
end for;
assert remaining eq [*[*1,4,4,0*],[*1,4,12,0*],[*2,2,1,82*],[*2,2,2,84*]*];
records:=[* *];
for d in [2,4] do
 xs:=allxs[d-1];rs:=allrs[d-1];K:=ext<GF(9)|d>;keys:=prior[d-1][2];
 if d eq 4 then
  for i in [4,12] do
   a:=Index(xs,rs[i]);ims:=M4TwentyOrbitKeys(X,rs[i],K);
   st:=[j:j in [1..96]|ims[j] eq keys[a]];assert #st eq 8;
   Append(~records,[*1,d,i,a,keys[a],ims,st*]);
  end for;
 else
  for pair in [[1,82],[2,84]] do
   i,j:=Explode(pair);a:=Index(xs,rs[i]);key:=Sort([keys[a],keys[j]]);
   ka:=M4TwentyOrbitKeys(X,xs[a],K);kj:=M4TwentyOrbitKeys(X,xs[j],K);
   ims:=[Sort([ka[g],kj[g]]):g in [1..96]];
   st:=[g:g in [1..96]|ims[g] eq key];assert #st eq 4;
   Append(~records,[*2,d,i,j,key,ims,st*]);
  end for;
 end if;
end for;
oldst:=eval Read("C4x20/data/M4N20_stabilizers_root_arithmetic.txt");
assert [*M4TwentyAutParameters(),records*] eq oldst;
Write(C420FullStem cat "_stabilizers.txt",Sprint(oldst):Overwrite:=true);
print "C420_FRESH_COMPLETE singles=18 pairs=192 survivors=4 stabilizers=8,8,4,4 complete=true",Cputime(tt);
quit;
