/* Fresh degree-eight relative-divisor local coverage. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
if Type(M) eq MonStgElt then M:=StringToInteger(M);end if;
if Type(N) eq MonStgElt then N:=StringToInteger(N);end if;
if Type(p) eq MonStgElt then p:=StringToInteger(p);end if;
if not assigned D then D:=8 div EulerPhi(M);elif Type(D) eq MonStgElt then D:=StringToInteger(D);end if;
f:=Order(Integers(M)!p);assert GCD(p,M*N) eq 1;
AttachSpec("mdmagma/v2/mdmagma.spec");
F:=GF(p^f);X:=MDX11(M,N,F);C:=Curve(X);cs:=Cusps(X);
print "OCTIC_MIXED_PREFLIGHT_MODEL",M,N,p,f,D,Genus(X),Sort([Degree(c):c in cs]),Cputime(t0);
counts:=[];
for d in [1..D] do
 ps:=NoncuspidalPlaces(X,d);Append(~counts,#ps);
 print "OCTIC_MIXED_PREFLIGHT_COUNT",M,N,p,f,d,#ps,Cputime(t0);
end for;
print "OCTIC_MIXED_PREFLIGHT_DONE",M,N,p,f,D,counts,"complete=true",Cputime(t0);quit;
