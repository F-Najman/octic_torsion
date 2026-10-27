/* Independent zeta/closed-place exhaustion check at all three primes. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
J:=JH(45,[16]);assert Dimension(J) eq 9;
for q in [4,7,13] do
 p:=q eq 4 select 2 else q;f:=q eq 4 select 2 else 1;
 fp:=FrobeniusPolynomial(J,p);assert Degree(fp) eq 18;
 traces:=[];
 for n in [1..4*f] do
  sn:=-n*Coefficient(fp,18-n);
  for j in [1..n-1] do sn-:=Coefficient(fp,18-j)*traces[n-j];end for;
  Append(~traces,Integers()!sn);
 end for;
 X:=MDX11(3,15,GF(q));cs:=Cusps(X);closed:=[];nc:=[];
 for d in [1..4] do
  cnt:=q^d+1-traces[f*d];
  for e in Divisors(d) do if e lt d then cnt-:=e*closed[e];end if;end for;
  assert cnt mod d eq 0;Append(~closed,cnt div d);
  x:=#NoncuspidalPlaces(X,d);Append(~nc,x);
  assert cnt div d eq x+#[c:c in cs|Degree(c) eq d];
 end for;
 print "M3N15_FROB_COUNTS_VERIFIED",q,nc,closed,Coefficients(fp),Cputime(t0);
end for;
print "M3N15_FROB_COUNTS_DONE complete=true",Cputime(t0);quit;
