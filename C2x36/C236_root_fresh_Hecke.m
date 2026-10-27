/* Exact C4 quotient of X1(2,36), paired with its9-point. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
A<t>:=FunctionField(GF(5));r:=t^2-t+1;c:=t*(r-1);b:=r*c;
E9:=EllipticCurve([A|1-c,-b,-b,0,0]);P9:=E9![0,0];assert 9*P9 eq E9!0 and 3*P9 ne E9!0;
j9:=jInvariant(E9);R<W>:=PolynomialRing(A);f:=16*(W^2+16)^3-j9*W^4;
assert IsIrreducible(f);K<w>:=FunctionField(f);assert Genus(K) eq 11;
print "C236_C4_QUOTIENT_GENUS11_VERIFIED",Genus(K),Cputime(t0);
a:=eval Read("C2x36/data/M2N36_D8_direct_full8_arithmetic.txt");
Rp:=PolynomialRing(GF(5));L<z>:=ext<GF(5)|Rp!a[1]>;ky:=a[8][1606][2];
b0:=Evaluate(Rp!ky[1..8],z);c0:=Evaluate(Rp!ky[9..16],z);
E:=EllipticCurve([L|0,c0,0,b0*(1-b0-c0),0]);P2:=E![0,0];Q:=E![b0,b0];
B,C:=Explode(TateNormalForm_bc(E,4*Q));tt:=C^2/(B-C);
T:=18*Q;U:=9*Q;third:=-c0-T[1];ww:=-4*(U[1]-T[1])/third;
assert 16*(ww^2+16)^3 eq jInvariant(E)*ww^4;
assert tt^25 eq tt and ww^25 eq ww and (tt^5 ne tt or ww^5 ne ww);
assert Evaluate(Numerator(j9),tt)/Evaluate(Denominator(j9),tt) eq jInvariant(E);
// g and g^2 preserve the invariant independently of any finite field hom.
for i in [1..4] do
 P2n:=P2+18*Q;Qn:=Q+P2;P2:=P2n;Q:=Qn;
 B,C:=Explode(TateNormalForm_bc(E,4*Q));tx:=C^2/(B-C);
 U:=9*Q;T:=18*Q;other:=P2+T;
 wx:=4*(U[1]-T[1])/(P2[1]-other[1]);assert tx eq tt and wx eq ww;
end for;
print "C236_C4_QUOTIENT_INVARIANT_MODULI_POINT",Eltseq(tt),Eltseq(ww),Cputime(t0);
cs:=Poles(K!j9);print "C236_C4_QUOTIENT_CUSPS",[Degree(p):p in cs],Cputime(t0);

Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
load "common/pair_af_place_linear.m";
function ZValues(E,P2,Q)
 B,C:=Explode(TateNormalForm_bc(E,4*Q));tt:=C^2/(B-C);
 U:=9*Q;T:=18*Q;other:=P2+T;ww:=4*(U[1]-T[1])/(P2[1]-other[1]);
 return [tt,ww];
end function;
vals:=ZValues(E,P2,Q);p:=C48PairAFPlaceLinear([K!t,w],vals);assert Degree(p) eq 2;
print "C236_C4_TARGET_PLACE_VERIFIED",Degree(p),Cputime(t0);
for ell in [7] do
 dl:=C48PairAFPlaceLinear([K!t,w],ZValues(E,P2,ell*Q));assert Degree(dl) eq 2;
 isogs:=MDIsogeniesFastPoints(E,[P2,Q],ell);TD:=DivisorGroup(K)!0;rows:=[* *];
 for item in isogs do
  EG,pts,k:=Explode(item);vals:=ZValues(EG,pts[1],pts[2]);yy:=C48PairAFPlaceLinear([K!t,w],vals);ee:=Degree(yy);
  assert (8*k) mod ee eq 0;mult:=(8*k) div ee;TD+:=mult*Divisor(yy);Append(~rows,[*k,ee,mult*]);
  print "C236_C4_HECKE_IMAGE",ell,k,ee,mult,Cputime(t0);
 end for;
 assert Degree(TD) eq 8*(ell+1);ss,mm:=Support(TD);assert &and[m mod 4 eq 0:m in mm];
 T2:=&+[(mm[i] div 4)*Divisor(ss[i]):i in [1..#ss]];assert Degree(T2) eq 2*(ell+1);
 A:=T2-Divisor(p)-ell*Divisor(dl);yes:=IsPrincipal(A);assert not yes;
 print "C236_C4_QUOTIENT_HECKE_NECESSARY",ell,yes,Cputime(t0);
 Write(Sprintf("C2x36/data/C236_final_Hecke_%o_arithmetic.txt",ell),Sprint([*ell,rows,yes*]):Overwrite:=true);
 if ell eq 7 then
  native:=MDIsogenies(E,ell);ND:=DivisorGroup(K)!0;assert #native eq ell+1;
  for phi in native do
   EE:=Domain(phi);pp:=phi(EE![BaseRing(EE)!a:a in Eltseq(P2)]);qq:=phi(EE![BaseRing(EE)!a:a in Eltseq(Q)]);
   yy:=C48PairAFPlaceLinear([K!t,w],ZValues(Codomain(phi),pp,qq));ee:=Degree(yy);assert 8 mod ee eq 0;
   ND+:=(8 div ee)*Divisor(yy);
  end for;
  assert ND eq TD;print "C236_C4_NATIVE_MD_ISOGENIES_REFERENCE_EQUAL",Degree(ND),Cputime(t0);
  assert #Basis(A:Reduction:=false,Simplification:="None") eq 0;
  assert Degree(K) eq 6 and Max(Degree(Numerator(j9)),Degree(Denominator(j9))) eq 36;
  print "C236_C4_RR_ZERO_INDEPENDENTLY_CONFIRMED",Cputime(t0);
 end if;
 if not yes then break;end if;
end for;
print "C236_ROOT_FRESH_CERTIFICATE_DONE complete=true",Cputime(t0);quit;
