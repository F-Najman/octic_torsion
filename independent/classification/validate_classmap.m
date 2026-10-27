// validate_classmap.m -- spot checks of Magma's class-group map used in classify8.m, independent of the
// class-group machinery: for every cusp and for random noncuspidal places Q of degree <= 8,
//   Q - deg(Q) c0 - m(cl(Q)) must be principal (Riemann--Roch, IsPrincipal), and
//   Q - deg(Q) c0 - m(cl(Q) + h) must NOT be principal for a random h != 0.
// Also: #J(F_p) from the class group vs. the Eichler--Shimura value from modular symbols.
// Usage: magma -n -b N:=33 p:=2 validate_classmap.m
SetColumns(0);
if not assigned N then N:=33; elif Type(N) eq MonStgElt then N:=StringToInteger(N); end if;
if not assigned p then p:=2; elif Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
if not assigned Samples then Samples:=10; elif Type(Samples) eq MonStgElt then Samples:=StringToInteger(Samples); end if;
SetMemoryLimit(16*10^9); SetSeed(1);
load "independent/classification/helpers8.m"; t0:=Cputime();
procedure note(s) printf "PROGRESS [validate N=%o p=%o] %o (cpu %o s)\n", N, p, s, Cputime(t0); end procedure;
F,b,c:=rawmodel(N,GF(p)); j:=jInvariant(TateCurve(b,c));
cusps:=[Q: Q in Support(Divisor(j)) | Valuation(j,Q) lt 0];
Cl,m:=ClassGroup(F); T:=TorsionSubgroup(Cl); hms:=JacOrderMS(N,p);
note(Sprintf("#J(F_p): class group %o, modular symbols %o, equal %o", #T, hms, #T eq hms));
c0:=[Q: Q in cusps | Degree(Q) eq 1][1];
test:=cusps; for d in [1..8] do L:=[Q: Q in Places(F,d) | Valuation(j,Q) ge 0]; if #L gt 0 then test cat:=[Random(L): i in [1..Min(Samples,#L)]]; end if; end for;
okpos:=0; okneg:=0;
for Q in test do D:=1*Q - Degree(Q)*c0; x:=D @@ m;
  if IsPrincipal(D - m(x)) then okpos+:=1; end if;
  repeat h:=Random(T); until h ne T!0;
  if not IsPrincipal(D - m(x + Cl!h)) then okneg+:=1; end if;
end for;
note(Sprintf("VALIDATE: %o divisors (all %o cusps + random noncuspidal places up to degree 8): representative principal %o/%o, wrong class non-principal %o/%o",
   #test, #cusps, okpos, #test, okneg, #test));
printf "VALIDATE_DONE N=%o p=%o\n", N, p; quit;
