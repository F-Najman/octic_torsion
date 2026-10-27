// classify8.m -- independent check (review of the Phi(8) draft, Props 3.87-3.92) of the classification of
// rational effective degree-8 divisors with noncuspidal support on X_1(N), by reduction at a good prime p.
// Own code; model = Sutherland's raw model FFFcN (via the mdmagma copy); known points from van Hoeij's
// LowDegreePlaces (~/degree5_torsion/refs/vanhoeij/LowDegreePlaces), coordinates (x,y) as there.
// Usage:  magma -n -b N:=33 p:=2 MemGB:=16 classify8.m
SetColumns(0);
if not assigned N then N:=33; elif Type(N) eq MonStgElt then N:=StringToInteger(N); end if;
if not assigned p then p:=2; elif Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
if not assigned MemGB then MemGB:=16; elif Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
load "independent/classification/helpers8.m"; t0:=Cputime(); w0:=Realtime();
procedure note(s) printf "PROGRESS [N=%o p=%o] %o (cpu %o s, wall %o s)\n", N, p, s, Cputime(t0), Realtime()-w0; end procedure;

DEHMZ:=AssociativeArray(); DEHMZ[33]:=[5,42373650]; DEHMZ[34]:=[8760,595680]; DEHMZ[35]:=[13,109148520];
DEHMZ[39]:=[7,31122,3236688]; DEHMZ[40]:=[2,2,2,8,120,895440];
GENUS:=AssociativeArray(); GENUS[33]:=21; GENUS[34]:=21; GENUS[35]:=25; GENUS[39]:=33; GENUS[40]:=25;

// ---- known points (van Hoeij), as <minimal data, label> ; each is [K, x, y] ----
Qz<Z>:=PolynomialRing(Rationals());
KNOWN:=[* *];   // octic points
LOWER:=[* *];   // lower-degree points (only N=33 has one of degree < 8)
if N eq 33 then
  K<a>:=NumberField(Z^7-6*Z^6+9*Z^5-5*Z^3-Z+1); Append(~LOWER, <K, a, -(3*a^6-16*a^5+17*a^4+9*a^3-11*a^2-3*a-4)/3, "septic">);
  K4<u>:=NumberField(Z^4-3*Z^3+3*Z^2-3*Z+1); K4y<Y>:=PolynomialRing(K4); K8r<v>:=NumberField(Y^2+(u^3-3*u^2+u-3)*Y+2*u^2-u+2);
  K:=AbsoluteField(K8r); Append(~KNOWN, <K, K!u, K!v, "octic, degj=4">);
elif N eq 34 then
  K<a>:=NumberField(Z^8-6*Z^7+16*Z^6-28*Z^5+42*Z^4-44*Z^3+27*Z^2-8*Z+1);
  Append(~KNOWN, <K, a, -(3*a^7+53*a^6-320*a^5+716*a^4-1036*a^3+1516*a^2-1197*a+64)/157, "f1, j=1728">);
  K<a>:=NumberField(Z^8+Z^7-Z^5-Z^4-Z^3+Z+1); Append(~KNOWN, <K, a, -(a^7-a^5-a^4-a^3+1), "f2, degj=2">);
  K<a>:=NumberField(Z^8-2*Z^6+8*Z^4-7*Z^2+1); Append(~KNOWN, <K, a, -(3*a^7-5*a^5+23*a^3-16*a)/3, "f3, degj=4">);
elif N eq 35 then
  K<a>:=NumberField(Z^8-10*Z^7+36*Z^6-55*Z^5+19*Z^4+35*Z^3-25*Z^2-5*Z+5);
  Append(~KNOWN, <K, a, -(6*a^7-65*a^6+268*a^5-531*a^4+473*a^3-49*a^2-177*a+74), "degj=4">);
elif N eq 39 then
  K<a>:=NumberField(Z^8-4*Z^7+10*Z^6-25*Z^5+37*Z^4-25*Z^3+10*Z^2-4*Z+1);
  Append(~KNOWN, <K, a, (2*a^7-6*a^6+13*a^5-34*a^4+33*a^3+a^2+3*a)/3, "j=0">);
elif N eq 40 then
  K<a>:=NumberField(Z^8-3*Z^7+3*Z^6-2*Z^5+13/5*Z^4-8/5*Z^3-7/5*Z^2+2*Z-2/5);
  Append(~KNOWN, <K, a, (5455*a^7-25010*a^6+27370*a^5-10105*a^4+24648*a^3-20135*a^2-10773*a+14311)/6653, "degj=4">);
end if;
Dreps:=[d: d in [1..N div 2] | Gcd(d,N) eq 1];

// ---- characteristic zero: exact order, distinct diamond translates, full torsion bound ----
function ptkeyQ(bb,cc) K:=Parent(bb); S:=sub<K|bb,cc>; d:=Degree(S);
  for l in [0..20] do g:=bb+l*cc; if Degree(MinimalPolynomial(g)) eq d then
     V:=Matrix(Rationals(),[Eltseq(K!g^e): e in [0..d-1]]); ok1,s1:=IsConsistent(V,Vector(Rationals(),Eltseq(bb))); ok2,s2:=IsConsistent(V,Vector(Rationals(),Eltseq(cc)));
     assert ok1 and ok2; return <d,l,Eltseq(MinimalPolynomial(g)),Eltseq(s1),Eltseq(s2)>; end if; end for;
  error "ptkeyQ"; end function;
for kn in KNOWN cat LOWER do K,xx,yy,lab:=Explode(kn); bK,cK:=xy2bc(xx,yy); EK:=TateCurve(bK,cK); PK:=EK![0,0];
  assert Order(PK) eq N;
  keys:={ ptkeyQ(b1,c1) where b1,c1:=tate(EK,dd*PK) : dd in Dreps };
  assert #{k[1]: k in keys} eq 1;
  // torsion bound from good unramified odd primes of degree 1
  O:=MaximalOrder(K); g:=0; used:=[];
  for q in PrimesInInterval(3,400) do if N mod q eq 0 or #used ge 8 then continue; end if;
    for pp in Factorization(q*O) do if pp[2] eq 1 and Valuation(Discriminant(EK),pp[1]) eq 0 and Min(Valuation(bK,pp[1]),Valuation(cK,pp[1])) ge 0 then
      kk,red:=ResidueClassField(pp[1]); h:=#TateCurve(red(bK),red(cK)); g:=Gcd(g,h); Append(~used,<q,Degree(kk),h>); break; end if; end for;
  end for;
  note(Sprintf("KNOWN %o: degree %o, exact order %o; %o distinct closed points among the %o diamond translates; gcd of #E(k_v) at good unramified odd v %o = %o", lab, Degree(K), N, #keys, #Dreps, used, g));
end for;

// ---- model over F_p ----
k:=GF(p); F,b,c:=rawmodel(N,k); E:=TateCurve(b,c); P0:=E![0,0]; j:=jInvariant(E);
note(Sprintf("genus %o (expected %o), deg j %o (expected %o)", Genus(F), GENUS[N], Degree(j), jdegree(N)));
assert Genus(F) eq GENUS[N] and Degree(j) eq jdegree(N);
for q in PrimeDivisors(N) do assert (N div q)*P0 ne E!0; end for; assert N*P0 eq E!0;
note("P=(0,0) has exact order N over the function field");

// ---- cusps and their grouping into reductions of rational cusp orbits ----
cusps:=[Q: Q in Support(Divisor(j)) | Valuation(j,Q) lt 0]; cdeg:=[Degree(Q): Q in cusps];
D2:=16*b^2+(1-20*c-8*c^2)*b+c*(c-1)^3; X4:=(4*P0)[1]/(4*P0)[3]; X5:=(5*P0)[1]/(5*P0)[3];
sig:=func<Q | <Valuation(b,Q),Valuation(c,Q),Valuation(b-c,Q),Valuation(D2,Q),Valuation(X4,Q),Valuation(X5,Q)> >;
sigs:=SetToSequence({sig(Q): Q in cusps}); groups:=[[i: i in [1..#cusps] | sig(cusps[i]) eq s]: s in sigs];
gdeg:=Sort([&+[cdeg[i]: i in G]: G in groups]);
note(Sprintf("F_p cusp degrees %o; rational orbit degrees %o; %o signature groups with degree sums %o, members %o",
  Sort(cdeg), cusporbitdegrees(N), #groups, gdeg, Sort([Sort([cdeg[i]: i in G]): G in groups])));
assert gdeg eq cusporbitdegrees(N);   // #groups = #rational orbits, so each group is the reduction of exactly one orbit

// ---- class group and the image H of J_1(N)(Q) (= rational cuspidal group, DEHMZ) ----
Cl,m:=ClassGroup(F); T:=TorsionSubgroup(Cl);
note(Sprintf("class group %o, #J(F_p) = %o", Invariants(Cl), #T));
c0i:=[G[1]: G in groups | #G eq 1 and cdeg[G[1]] eq 1][1]; c0:=cusps[c0i];
cl:=func<D | (D - Degree(D)*c0) @@ m>;
cc:=[cl(1*Q): Q in cusps];
H:=sub<Cl | [&+[cc[i]: i in G]: G in groups]>;
note(Sprintf("H (image of rational cuspidal classes) %o, order %o; DEHMZ J_1(N)(Q) = %o of order %o; equal orders: %o; index [J(F_p):H] = %o",
   Invariants(H), #H, DEHMZ[N], &*DEHMZ[N], #H eq &*DEHMZ[N], #T div #H));
injective := #H eq &*DEHMZ[N];
Qt,q:=quo<Cl|H>;
vcl:=func<v | &+[Cl| v[i]*cc[i] : i in [1..#v] | v[i] ne 0]>;          // class of a cuspidal multiplicity vector
liftable:=func<v | forall{G: G in groups | #{v[i]: i in G} eq 1}>;

// ---- noncuspidal places of degree <= 8 ----
nonc:=[]; for d in [1..8] do Append(~nonc,[Q: Q in Places(F,d) | Valuation(j,Q) ge 0]); note(Sprintf("noncuspidal places of degree %o: %o", d, #nonc[d])); end for;
note(Sprintf("NONCUSPIDAL_COUNTS %o", [#L: L in nonc]));
keyof:=AssociativeArray(); for d in [1..8] do for Q in nonc[d] do kq:=Sprint(placekey(Q,b,c)); assert not IsDefined(keyof,kq); keyof[kq]:=Q; end for; end for;

// ---- specialisations of the known points and their diamond translates ----
function reddiv(K,xx,yy,dd) bK,cK:=xy2bc(xx,yy); sp,cus:=specialise(K,bK,cK,N,p : Diamond:=dd); assert not cus;
  return &+[ t[2]*keyof[Sprint(t[1])] : t in sp ]; end function;
KnownOct:=[]; for kn in KNOWN do K,xx,yy,lab:=Explode(kn); L:={ reddiv(K,xx,yy,dd) : dd in Dreps };
  note(Sprintf("KNOWN %o: reductions of the diamond translates: %o distinct divisors, support degrees and multiplicities %o",
       lab, #L, [ [<Degree(Q),Valuation(D,Q)>: Q in Support(D)] : D in L]));
  KnownOct cat:= SetToSequence(L); end for;
KnownLow:=[]; for kn in LOWER do K,xx,yy,lab:=Explode(kn); L:={ reddiv(K,xx,yy,dd) : dd in Dreps };
  note(Sprintf("LOWER %o: reductions of diamond translates: %o distinct divisors of degrees %o", lab, #L, {Degree(D): D in L})); KnownLow cat:= SetToSequence(L); end for;
ratcusps:=[cusps[G[1]]: G in groups | #G eq 1 and cdeg[G[1]] eq 1];
// known global rational effective degree-8 divisors with noncuspidal support: octic translates, and lower-degree translates plus rational cuspidal complements
KnownNCdivs:=KnownOct;
for D in KnownLow do e:=8-Degree(D); L:=[]; effvecs([1: i in ratcusps],1,e,[0: i in ratcusps],~L);
   for v in L do Append(~KnownNCdivs, D + &+[DivisorGroup(F)| v[i]*ratcusps[i]: i in [1..#ratcusps]]); end for; end for;
KnownNCcl:={cl(D): D in KnownNCdivs};
note(Sprintf("known global noncuspidal-support degree-8 divisors: %o (octic reductions %o), distinct classes %o; all in H: %o",
   #KnownNCdivs, #KnownOct, #KnownNCcl, forall{x: x in KnownNCcl | x in H}));

// ---- cuspidal parts: all effective cuspidal divisors of degree e <= 8, indexed by class mod H ----
CV:=[]; CI:=[]; Lift8:={};
for e in [0..8] do L:=[]; effvecs(cdeg,1,e,[0: i in cusps],~L); A:=AssociativeArray();
  for v in L do x:=vcl(v); key:=Eltseq(q(x)); if IsDefined(A,key) then Append(~A[key],v); else A[key]:=[v]; end if;
     if e eq 8 and liftable(v) then Include(~Lift8,x); end if; end for;
  Append(~CV,L); Append(~CI,A); note(Sprintf("effective cuspidal divisors of degree %o: %o (classes mod H: %o)", e, #L, #Keys(A)));
end for;
// entirely cuspidal reductions
zkey:=Eltseq(Qt!0); inH8:=IsDefined(CI[9],zkey) select CI[9][zkey] else [];
badc:=[v: v in inH8 | not liftable(v) and vcl(v) notin Lift8];
note(Sprintf("ENTIRELY_CUSPIDAL degree 8: %o divisors, %o liftable classes, %o with class in H, %o of these non-liftable, %o non-liftable in H with class outside the liftable classes",
   #CV[9], #Lift8, #inH8, #[v: v in inH8 | not liftable(v)], #badc));

// ---- noncuspidal parts: all multisets of noncuspidal places of total degree 1..8 ----
NC:=&cat nonc; ncl:=[]; i:=0;
for Q in NC do i+:=1; Append(~ncl, cl(1*Q)); if i mod 100 eq 0 then note(Sprintf("classes of noncuspidal places: %o of %o", i, #NC)); end if; end for;
note(Sprintf("classes of all %o noncuspidal places computed", #NC));
parts:=[]; // multisets of noncuspidal places, as lists of <index into NC, multiplicity>
procedure rec(start, rem, cur, ~parts)
  for i in [start..#NC] do dQ:=Degree(NC[i]); if dQ gt rem then break; end if;
    c2:=cur; if #c2 gt 0 and c2[#c2][1] eq i then c2[#c2][2]+:=1; else Append(~c2,<i,1>); end if;
    Append(~parts,c2); rec(i, rem-dQ, c2, ~parts);
  end for;
end procedure;
rec(1, 8, [], ~parts);
note(Sprintf("noncuspidal parts of degree <= 8: %o", #parts));
tested:=AssociativeArray(); surv:=[];
for A in parts do dA:=&+[t[2]*Degree(NC[t[1]]): t in A]; xA:=&+[t[2]*ncl[t[1]]: t in A];
  typ:=Sort(&cat[[Degree(NC[t[1]]): r in [1..t[2]]]: t in A]); e:=8-dA;
  tt:=IsDefined(tested,typ) select tested[typ] else [0,0]; tt[1]+:=#CV[e+1];
  key:=Eltseq(-q(xA)); if IsDefined(CI[e+1],key) then for v in CI[e+1][key] do tt[2]+:=1; Append(~surv,<A,v>); end for; end if;
  tested[typ]:=tt;
end for;
for typ in Sort(SetToSequence(Keys(tested))) do note(Sprintf("TEST noncuspidal type %o: %o divisors tested, %o with class in H", typ, tested[typ][1], tested[typ][2])); end for;

// ---- identification of the survivors ----
DG:=DivisorGroup(F);
survdiv:=[]; for s in surv do A:=s[1]; v:=s[2]; Append(~survdiv, &+[DG| t[2]*NC[t[1]]: t in A] + &+[DG| v[i]*cusps[i]: i in [1..#cusps] | v[i] ne 0]); end for;
survcl:=[cl(D): D in survdiv];
expl:=[x in KnownNCcl or x in Lift8 : x in survcl];
note(Sprintf("SURVIVORS: %o divisors (noncuspidal support types %o), %o distinct classes; every survivor class is the reduction class of a known global effective degree-8 divisor: %o",
   #surv, {* Sort(&cat[[Degree(NC[t[1]]): r in [1..t[2]]]: t in s[1]]) : s in surv *}, #Set(survcl), &and expl));
octcl:={cl(D): D in KnownOct};
note(Sprintf("survivor classes equal to classes of reductions of known octic translates: %o of %o known octic classes; survivors equal (as divisors) to a known octic reduction: %o",
   #(Set(survcl) meet octcl), #octcl, #[D: D in survdiv | D in KnownOct]));
note(Sprintf("CONCLUSION N=%o: injective=%o, bad entirely-cuspidal=%o, unexplained survivors=%o", N, injective, #badc, #[i: i in [1..#surv] | not expl[i]]));
printf "CLASSIFY8_DONE N=%o p=%o (cpu %o s)\n", N, p, Cputime(t0);
quit;
