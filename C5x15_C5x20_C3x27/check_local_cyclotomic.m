/* Exhaustive finite elliptic-curve obstruction, no earlier-degree input. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);t0:=Cputime();
for row in [<5,15,11,1,2>,<5,20,11,1,2>,<3,27,2,2,4>] do
 m,n,p,f,dmax:=Explode(row);
 for d in [1..dmax] do
  F:=GF(p^(f*d));hit:=0;cardhits:=0;invs:={};curves:=0;
  for j in F do
   E0:=EllipticCurveFromjInvariant(j);
   for E in Twists(E0) do
    curves+:=1;
    if #E mod (m*n) ne 0 then continue;end if;
    cardhits+:=1;G,mp:=AbelianGroup(E);v:=Invariants(G);Include(~invs,v);
    if #v eq 2 and v[1] mod m eq 0 and v[2] mod n eq 0 then hit+:=1;end if;
   end for;
  end for;
  assert hit eq 0;
  print "OCTIC_LOCAL_CYCLOTOMIC",m,n,p,f,d,#F,curves,cardhits,invs,"hits",hit,Cputime(t0);
 end for;
end for;
print "OCTIC_LOCAL_CYCLOTOMIC_DONE cases=8 hits=0 complete=true",Cputime(t0);quit;
