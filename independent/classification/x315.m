// x315.m -- independent check of Lemmas 3.82, 3.83 and Prop. 3.84 of the Phi(8) draft:
// X = X_1(3,15) over F = Q(zeta_3), F-rational effective degree-4 divisors with noncuspidal support.
// Own code. Model: Sutherland's X1_3_15 ("Xz") from the mdmagma copy, with
//   q = u+1, t = v*u/(u+1), E: y^2 + ((1+q)t+2-q) xy + (q t^2+(1-q)t) y = x^3, P3 = (0,0) (order 3), P15 = (-t,t^2) (order 15).
// (q,t) are canonical moduli invariants (P3 at the origin as flex with tangent y=0; scaling fixed by P15 = (-t,t^2)).
SetColumns(0);
if not assigned MemGB then MemGB:=16; elif Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
load "independent/classification/helpers8.m"; t0:=Cputime(); w0:=Realtime();
procedure note(s) printf "PROGRESS [X1(3,15)] %o (cpu %o s, wall %o s)\n", s, Cputime(t0), Realtime()-w0; end procedure;

str:=Read("mdmagma/models_X1_m_n/X1_3_15.txt");
xzline:=[l: l in Split(str,"\n") | #l ge 2 and l[1..2] eq "Xz"][1]; xzstr:=Split(Split(xzline,"=")[2],";")[1];
R3<u,v,z>:=PolynomialRing(Integers(),3); f3:=eval xzstr;

// canonical (q,t) of a triple (E, P3, P15) over a field
function qt_of(E, P3, P15)
  k:=BaseRing(E); a1,a2,a3,a4,a6:=Explode(aInvariants(E)); x0:=P3[1]; y0:=P3[2];
  R<X,Y>:=PolynomialRing(k,2); f:=Y^2+a1*X*Y+a3*Y-(X^3+a2*X^2+a4*X+a6);
  g:=Evaluate(f,[X+x0,Y+y0]); A3:=MonomialCoefficient(g,Y); A4:=-MonomialCoefficient(g,X); assert A3 ne 0;
  g:=Evaluate(g,[X,Y+(A4/A3)*X]);
  A1:=MonomialCoefficient(g,X*Y); A2:=-MonomialCoefficient(g,X^2); A3:=MonomialCoefficient(g,Y);
  assert A2 eq 0 and MonomialCoefficient(g,X) eq 0 and MonomialCoefficient(g,R!1) eq 0;
  x1:=P15[1]-x0; y1:=(P15[2]-y0)-(A4/A3)*x1;
  t:=-y1^2/x1^3; mu:=x1^2/y1; a1n:=A1/mu; a3n:=A3/mu^3; q:=(a1n-t-2)/(t-1);
  assert a3n eq q*t^2+(1-q)*t; return q,t;
end function;
function curve_qt(q,t) E:=EllipticCurve([(1+q)*t+(2-q),0,q*t^2+(1-q)*t,0,0]); return E, E![0,0], E![-t,t^2]; end function;

// the curve over a finite field k with zeta_3 -> zk
function model(k, zk)
  Fu<U>:=RationalFunctionField(k); PV<V>:=PolynomialRing(Fu);
  g:=Evaluate(f3,[PV!U, V, PV!zk]); assert IsIrreducible(g);
  F:=FunctionField(g); uu:=F!U; vv:=F.1; q:=uu+1; t:=vv*uu/(uu+1); return F,q,t;
end function;

// 48 = (SL_2(F_3) x (Z/5)^*)/{+-1} automorphisms, acting on triples: P3' = a P3 + b 5P15, P15' = 2(c P3 + d 5P15) + 2w 3P15
AUT:=[<a,b,c,d,w>: a,b,c,d in [0..2], w in [1..4] | (a*d-b*c) mod 3 eq 1];
function act(g, E, P3, P15) a,b,c,d,w:=Explode(g); P5:=5*P15; P3b:=3*P15;
  return a*P3+b*P5, 2*(c*P3+d*P5)+2*w*P3b; end function;

// ---------------- one residue field ----------------
function setup(k, zk, label)
  F,q,t:=model(k,zk); E,P3,P15:=curve_qt(q,t); j:=jInvariant(E);
  gF:=Genus(F); dj:=Degree(j);
  assert 3*P3 eq E!0 and P3 ne E!0 and 15*P15 eq E!0 and 5*P15 ne E!0 and 3*P15 ne E!0;
  // signature functions: Tate normal form (b,c) of (E, P15 + i P3), i = 0,1,2 (units on the open curve)
  funs:=[]; for i in [0..2] do bb,cc:=tate(E, P15+i*P3); funs cat:=[bb,cc,bb-cc,16*bb^2+(1-20*cc-8*cc^2)*bb+cc*(cc-1)^3]; end for;
  cusps:=[Q: Q in Support(Divisor(j)) | Valuation(j,Q) lt 0];
  sig:=func<Q | [Valuation(f,Q): f in funs]>;
  return F,q,t,E,P3,P15,j,gF,dj,cusps,[sig(Q): Q in cusps];
end function;

F4:=GF(4); om:=F4.1; assert om^2+om+1 eq 0;
F4m,q4,t4,E4,P34,P154,j4,g4,dj4,cusps4,sig4:=setup(F4,om,"F4");
note(Sprintf("F4 model: genus %o, deg j %o (expected 3*96=288); P3 order 3, P15 order 15", g4, dj4));
F7m,q7,t7,E7,P37,P157,j7,g7,dj7,cusps7,sig7:=setup(GF(7),GF(7)!2,"F7");
note(Sprintf("F7 model (zeta_3 -> 2): genus %o, deg j %o", g7, dj7));
assert g4 eq g7 and dj4 eq 288 and dj7 eq 288;

// F-rational cusp orbits: at 7 Frobenius generates Gal(Q(zeta_15)/F), so F7 cusps = F-orbits; match with F4 by signatures
S7:=SetToSequence(Set(sig7)); assert #S7 eq #cusps7;   // signatures separate the F7 cusps
note(Sprintf("F7 cusp degrees %o (= F-rational cusp orbit degrees); F4 cusp degrees %o", Sort([Degree(Q): Q in cusps7]), Sort([Degree(Q): Q in cusps4])));
assert Set(sig4) eq Set(sig7);
grp4:=[[i: i in [1..#cusps4] | sig4[i] eq s]: s in S7]; grp7:=[[Index(sig7,s)]: s in S7];
assert forall{i: i in [1..#S7] | &+[Degree(cusps4[r]): r in grp4[i]] eq Degree(cusps7[grp7[i][1]])};
note(Sprintf("%o F-rational cusp orbits, matched between F4 and F7 by unit signatures; F4 splitting: %o", #S7, [[Degree(cusps4[r]): r in G]: G in grp4]));

// class groups
Cl4,m4:=ClassGroup(F4m); note(Sprintf("Cl(F4 model) %o, #J(F4) = %o", Invariants(Cl4), #TorsionSubgroup(Cl4)));
Cl7,m7:=ClassGroup(F7m); note(Sprintf("Cl(F7 model) %o, #J(F7) = %o", Invariants(Cl7), #TorsionSubgroup(Cl7)));
i0:=[i: i in [1..#S7] | #grp4[i] eq 1 and Degree(cusps4[grp4[i][1]]) eq 1][1];
c04:=cusps4[grp4[i0][1]]; c07:=cusps7[grp7[i0][1]];
cl4:=func<D | (D - Degree(D)*c04) @@ m4>; cl7:=func<D | (D - Degree(D)*c07) @@ m7>;
orb4:=[&+[1*cusps4[r]: r in G]: G in grp4]; orb7:=[1*cusps7[G[1]]: G in grp7];
L:=FreeAbelianGroup(#S7);
f4:=hom<L -> Cl4 | [cl4(D): D in orb4]>; f7:=hom<L -> Cl7 | [cl7(D): D in orb7]>;
H4:=Image(f4); H7:=Image(f7);
note(Sprintf("cuspidal group C: image at F7 %o (order %o), image at F4 %o (order %o); kernels of Z^20 -> J(F4), J(F7) equal: %o",
   Invariants(H7), #H7, Invariants(H4), #H4, Kernel(f4) eq Kernel(f7)));
assert Kernel(f4) eq Kernel(f7);
psi:=func<y | f7(y @@ f4)>;   // H4 -> H7, the identification through C

// #J at further split primes (13: zeta->3, 19: zeta->7)
ords:=[#TorsionSubgroup(Cl7)];
if not assigned SkipJ then SkipJ:=0; elif Type(SkipJ) eq MonStgElt then SkipJ:=StringToInteger(SkipJ); end if;
for pr in [pr0: pr0 in [<13,3>,<19,7>] | SkipJ ne 1] do Fm,_,_:=model(GF(pr[1]),GF(pr[1])!pr[2]); Clp:=ClassGroup(Fm); Append(~ords,#TorsionSubgroup(Clp));
  note(Sprintf("#J(F_%o) = %o (%o)", pr[1], ords[#ords], Invariants(Clp))); end for;
gg:=Gcd(ords); note(Sprintf("#J(F_7), #J(F_13), #J(F_19) = %o; gcd %o = %o * #C; so [J(F):C] divides %o (given rank J(F)=0)", ords, gg, gg/#H7, gg div #H7));

// ---------------- enumeration of effective degree-4 divisors ----------------
function allplaces(Fm, j, cusps)
  NC:=[]; for d in [1..4] do NC cat:=[Q: Q in Places(Fm,d) | Valuation(j,Q) ge 0]; end for; return NC; end function;
NC4:=allplaces(F4m,j4,cusps4); NC7:=allplaces(F7m,j7,cusps7);
note(Sprintf("noncuspidal place counts deg 1..4: F4 %o, F7 %o", [#[Q: Q in NC4 | Degree(Q) eq d]: d in [1..4]], [#[Q: Q in NC7 | Degree(Q) eq d]: d in [1..4]]));
// all effective degree-4 divisors: multisets of noncuspidal places + cuspidal complement
procedure recnc(NC, start, rem, cur, ~parts)
  for i in [start..#NC] do dQ:=Degree(NC[i]); if dQ gt rem then continue; end if;
    c2:=cur+1*NC[i]; Append(~parts,c2); recnc(NC, i, rem-dQ, c2, ~parts); end for;
end procedure;
function alldivs(NC, cusps, wantnc)
  out:=[]; DG:=Parent(1*cusps[1]); parts:=[DG!0]; recnc(NC,1,4,DG!0,~parts);
  for A in parts do if wantnc and A eq DG!0 then continue; end if; L:=[]; effvecs([Degree(Q): Q in cusps],1,4-Degree(A),[0: Q in cusps],~L);
    for vv in L do Append(~out, A + &+[DG| vv[i]*cusps[i]: i in [1..#cusps] | vv[i] ne 0]); end for; end for;
  return out;
end function;
D4:=alldivs(NC4,cusps4,false); D7:=alldivs(NC7,cusps7,true);
note(Sprintf("effective degree-4 divisors: F4 all %o (entirely cuspidal %o); F7 with noncuspidal support %o",
   #D4, #[D: D in D4 | forall{Q: Q in Support(D) | Q in cusps4}], #D7));
C4:=[cl4(D): D in D4]; note("F4 classes computed");
note(Sprintf("GONALITY_CHECK: the %o effective degree-4 divisors over F4 have %o distinct classes (equal => gon over F4, hence over F, exceeds 4)", #D4, #Set(C4)));
C7:=[cl7(D): D in D7]; note("F7 classes computed");
// condition 3x in C:  3*c4 in H4 and 3*c7 = psi(3*c4)
A4:=AssociativeArray(); for i in [1..#D4] do y:=3*C4[i]; if y in H4 then key:=Eltseq(psi(y)); if IsDefined(A4,key) then Append(~A4[key],i); else A4[key]:=[i]; end if; end if; end for;
pairs:=[]; for i in [1..#D7] do key:=Eltseq(3*C7[i]); if IsDefined(A4,key) then for i4 in A4[key] do Append(~pairs,<i4,i>); end for; end if; end for;
clpairs:={<C4[pp[1]],C7[pp[2]]>: pp in pairs};
note(Sprintf("pairs (R4,R7) with 3[x] in C at both primes: %o divisor pairs, %o class pairs; all their F4 classes in H4 (=> [D-4c0] in C): %o",
   #pairs, #clpairs, forall{x: x in clpairs | x[1] in H4}));
fin:={x: x in clpairs | x[1] in H4 and x[2] eq psi(x[1])};
xs:={x[1]: x in fin};
note(Sprintf("TWELVE: class pairs with [x] in C at both primes: %o, distinct elements of C: %o", #fin, #xs));
R7reps:=[D7[pp[2]]: pp in pairs | <C4[pp[1]],C7[pp[2]]> in fin];
note(Sprintf("noncuspidal support types at F7 of the surviving R7: %o", {* Sort([Degree(Q): Q in Support(D) | Q notin cusps7]) : D in R7reps *}));

// ---------------- automorphism action at F7 (on noncuspidal places, via moduli) ----------------
keyof7:=AssociativeArray(); for Q in NC7 do kk,mk:=ResidueClassField(Q); key,dd:=ptkey(mk(q7),mk(t7)); assert dd eq Degree(Q); keyof7[Sprint(key)]:=Q; end for;
function gplace(g, Q) kk,mk:=ResidueClassField(Q); E,P3,P15:=curve_qt(mk(q7),mk(t7)); P3n,P15n:=act(g,E,P3,P15);
  qn,tn:=qt_of(E,P3n,P15n); key,dd:=ptkey(qn,tn); return keyof7[Sprint(key)]; end function;
function gdiv(g, D) return &+[ Valuation(D,Q)*gplace(g,Q) : Q in Support(D)]; end function;
pure:=[D: D in R7reps | forall{Q: Q in Support(D) | Q notin cusps7}];
note(Sprintf("%o of the surviving F7 divisors are purely noncuspidal", #pure));
xs7:={psi(x): x in xs};
R0:=pure[1]; orbcl:={cl7(gdiv(g,R0)): g in AUT};
note(Sprintf("ORBIT: the 96 (=48 up to sign) automorphisms applied to one surviving divisor give %o classes; equal to the %o surviving classes: %o",
   #orbcl, #xs7, orbcl eq xs7));

// ---------------- the known point (Prop. 3.84) ----------------
QQz<Zz>:=PolynomialRing(Rationals()); Fz<zz>:=NumberField(Zz^2+Zz+1); Fw<W>:=PolynomialRing(Fz);
Krel<w>:=NumberField(W^4-5*W^3+8*W^2-6*W+3); assert IsIrreducible(W^4-5*W^3+8*W^2-6*W+3);
a1:=((-7*zz-8)*w^3+(23*zz+34)*w^2+(-8*zz-43)*w)/6 + (-3*zz+10)/2;
a3:=((17*zz+13)*w^3+(-82*zz-71)*w^2+(106*zz+119)*w)/12 - 2*zz - 5;
xQ:=((-zz+1)*w^3+(8*zz-2)*w^2+(-17*zz-1)*w)/6 + 3*zz/2;
yQ:=((zz+2)*w^3+(-4*zz-8)*w^2+(3*zz+6)*w+2*zz+3)/4;
K:=AbsoluteField(Krel); zK:=K!zz;
EK:=EllipticCurve([K!a1,0,K!a3,0,0]); PK:=EK![0,0]; QK:=EK![K!xQ,K!yQ];
assert Order(PK) eq 3 and Order(QK) eq 15;
note(Sprintf("known point: P3 order 3, P15 order 15, j = %o, Weil pairing e3(P3,5P15) = zeta^%o", jInvariant(EK), [i: i in [1..2] | WeilPairing(PK,5*QK,3) eq zK^i]));
// which sign of P3 puts the triple on the model with zeta_3 -> z
cand:=[]; for sgn in [1,-1] do qK,tK:=qt_of(EK, sgn*PK, QK); uK:=qK-1; vK:=tK*qK/(qK-1);
  if Evaluate(f3,[uK,vK,zK]) eq 0 then Append(~cand,<sgn,qK,tK>); end if; end for;
assert #cand eq 1; sgn,qK,tK:=Explode(cand[1]);
note(Sprintf("the triple (E, %oP3, P15) lies on the model; F(q,t) has degree %o over Q (so the point is a closed degree-4 point over F)", sgn, Degree(sub<K|qK,tK,zK>)));
OK:=MaximalOrder(K);
function redknown(p, zval, keyfun)   // specialisation at the primes of K above (p, zeta - zval) (zval=0: all primes above p); cuspidal primes reported separately
  out:=[]; cus:=[]; for pp in Factorization(p*OK) do kk,red:=ResidueClassField(pp[1]);
    if zval cmpne 0 and red(zK) ne zval then continue; end if;
    if Valuation(qK,pp[1]) lt 0 or Valuation(tK,pp[1]) lt 0 then Append(~cus,<pp[2],Degree(kk)>); continue; end if;
    qb:=red(qK); tb:=red(tK); Eb,_,_:=curve_qt(qb,tb); if Discriminant(Eb) eq 0 then Append(~cus,<pp[2],Degree(kk)>); continue; end if;
    key,dd:=ptkey(qb,tb); Append(~out,<keyfun(key), pp[2]*Degree(kk) div dd>); end for; return out, cus; end function;
sp7,cus7:=redknown(7, GF(7)!2, func<key | keyof7[Sprint(key)]>);
assert #cus7 eq 0; Dk7:=&+[t[2]*t[1]: t in sp7];
jr:=Rationals()!jInvariant(EK); note(Sprintf("known point: j = %o, valuations at 2 and 7: %o, %o", jr, Valuation(jr,2), Valuation(jr,7)));
xk7:=cl7(Dk7);
note(Sprintf("KNOWN: reduction at the prime (7, zeta-2): %o (noncuspidal); its class [D-4c0] is one of the twelve (F7 classes): %o",
  [<Degree(Q),Valuation(Dk7,Q)>: Q in Support(Dk7)], xk7 in xs7));
keyof4:=AssociativeArray(); for Q in NC4 do kk,mk:=ResidueClassField(Q); key,dd:=ptkey(mk(q4),mk(t4)); assert dd eq 2*Degree(Q); keyof4[Sprint(key)]:=Q; end for;
sp4,cus4:=redknown(2, 0, func<key | keyof4[Sprint(key)]>);
note(Sprintf("known point at 2: noncuspidal parts %o, primes (e,f) reducing to cusps %o", [<Degree(t[1]),t[2]>: t in sp4], cus4));
if #cus4 eq 0 then Dk4:=&+[t[2]*t[1]: t in sp4]; note(Sprintf("F4 class in the twelve: %o", cl4(Dk4) in xs)); end if;
// full torsion of the known curve
g:=0; used:=[]; for pq in PrimesInInterval(7,300) do if #used ge 8 then break; end if;
  for pp in Factorization(pq*OK) do if pp[2] eq 1 and Valuation(Discriminant(EK),pp[1]) eq 0 and Valuation(K!a1,pp[1]) ge 0 and Valuation(K!a3,pp[1]) ge 0 then
     kk,red:=ResidueClassField(pp[1]); Eb:=EllipticCurve([red(K!a1),0,red(K!a3),0,0]); h:=#Eb; g:=Gcd(g,h);
     Append(~used,<pq,Degree(kk),h,Invariants(AbelianGroup(Eb))>); break; end if; end for; end for;
note(Sprintf("full torsion: #E(k_v) at good unramified odd v: %o; gcd %o", used, g));
printf "X315_DONE (cpu %o s)\n", Cputime(t0); quit;
