// #J_1(N)(F_p) via Eichler--Shimura on modular symbols, and the torsion bound gcd_{p odd, p not | N} #J_1(N)(F_p).
SetColumns(0); load "independent/classification/helpers8.m"; t0:=Cputime();
DEHMZ:=AssociativeArray(); DEHMZ[33]:=[5,42373650]; DEHMZ[34]:=[8760,595680]; DEHMZ[35]:=[13,109148520];
DEHMZ[39]:=[7,31122,3236688]; DEHMZ[40]:=[2,2,2,8,120,895440];
for N in [33,34,35,39,40] do
  tab:=[]; g:=0;
  for p in [q: q in PrimesUpTo(43) | N mod q ne 0] do h,dim:=JacOrderMS(N,p); Append(~tab,<p,h>); if p ne 2 then g:=Gcd(g,h); end if; end for;
  printf "MS N=%o genus %o: #J(F_p) = %o\n", N, dim, tab;
  printf "MS N=%o: gcd over odd good p up to 43 = %o = %o; DEHMZ #J(Q) = %o; equal: %o\n", N, g, Factorization(g), &*DEHMZ[N], g eq &*DEHMZ[N];
end for;
printf "MS_ORDERS_DONE (%o s)\n", Cputime(t0); quit;
