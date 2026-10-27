/* Standalone global completion of the H4-stabilized cases.
   This source needs only Magma, no earlier-degree files or retained state.
   A map alone does not force an octic point's image to be Qi-rational:
   the separately certified H4-invariance is essential. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
K<i>:=QuadraticField(-1);
function CartanW(E,P,Q)
 assert 4*P eq E!0 and 4*Q eq E!0 and 2*P ne E!0 and 2*Q ne E!0 and 2*P ne 2*Q;
 R:=-Q;U:=P-Q;T:=2*U;P2:=2*R;assert 2*P2 eq E!0 and P2 ne T;
 return 4*(U[1]-T[1])/(P2[1]-(P2+T)[1]);
end function;
A<s>:=FunctionField(K);aa:=2*s/(1+s^2);bb:=(1-s^2)/(1+s^2);assert aa^2+bb^2 eq 1;
E4:=EllipticCurve([A|0,-(1+aa^2),0,aa^2,0]);
P4:=E4![aa,i*aa*(aa-1)];Q4:=E4![1+bb,bb*(1+bb)];
assert 2*P4 eq E4![0,0] and 2*Q4 eq E4![1,0];
ww:=CartanW(E4,P4,Q4);
assert jInvariant(E4)*ww^4 eq 16*(ww^2+16)^3;
assert Maximum(Degree(Numerator(ww)),Degree(Denominator(ww))) eq 4;
for M in [[3,2,2,3],[0,3,1,0]] do
 a,b,c,d:=Explode(M);PP:=a*P4+b*Q4;QQ:=c*P4+d*Q4;
 assert CartanW(E4,PP,QQ) eq ww;
end for;
print "C420_CARTAN_INVARIANT_PASS degree_on_X4=4 jdegree=6 w=",ww;

// Universal5-point family and its cyclic5 invariant, checked physically.
function TateBC(E,P)
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));x:=P[1];y:=P[2];
 aa1:=a1;aa3:=2*y+a3+a1*x;aa2:=3*x+a2;aa4:=3*x^2+2*x*a2+a4-a1*y;
 aaa1:=2*aa4/aa3+aa1;aaa2:=aa2-(aa4/aa3)^2-aa1*aa4/aa3;
 b:=-aaa2^3/aa3^2;c:=-(aaa1*aaa2-aa3)/aa3;return b,c;
end function;
B<b>:=FunctionField(K);E5:=EllipticCurve([B|1-b,-b,-b,0,0]);P5:=E5![0,0];assert 5*P5 eq E5!0;
tt5:=(b^2-11*b-1)/b;assert jInvariant(E5) eq (tt5^2+10*tt5+5)^3/tt5;
for a in [1..4] do
 ba,ca:=TateBC(E5,a*P5);assert ba eq ca;
 assert (ba^2-11*ba-1)/ba eq tt5;
end for;
print "C420_CYCLIC5_INVARIANT_PASS t=b-11-1/b";

// Compact genus1 fiber product and exact inverse birational maps.
At<t>:=FunctionField(K);R<W>:=PolynomialRing(At);j5:=(t^2+10*t+5)^3/t;
f:=16*(W^2+16)^3-j5*W^4;assert IsIrreducible(f);F<w>:=FunctionField(f);assert Genus(F) eq 1;
u:=w^2*(t^2+10*t+5)/(4*(w^2+16));
ss:=(u^2*(2*t+10)-16)/(4*(u^2+2*u-4));X:=2*ss+1;Y:=4*w*(X+2)/(X+1)^2;
assert Y^2 eq (X-3)*(X+1)*(X+2);
assert (X-3)^2*(X+2)/(X+1)^2 eq t;
assert Y*(X+1)^2/(4*(X+2)) eq w;
print "C420_COMPACT_BIRATIONAL_PASS genus=1 model_y2=x3-7x-6";

EQ:=EllipticCurve([Rationals()|0,0,0,-7,-6]);EQtw:=QuadraticTwist(EQ,-1);
lo,hi:=RankBounds(EQ);lo2,hi2:=RankBounds(EQtw);assert lo eq 0 and hi eq 0 and lo2 eq 0 and hi2 eq 0;
E:=BaseChange(EQ,K);T,tm:=TorsionSubgroup(E);assert Invariants(T) eq [2,4];
E9:=EllipticCurve([GF(9)|0,0,0,-7,-6]);G9,m9:=AbelianGroup(E9);assert Invariants(G9) eq [4,4];
assert not IsSquare(K!5) and not IsSquare(K!(-5));
print "C420_ELLIPTIC_RANK_TORSION_PASS ranks=[0,0] Qi_torsion=[2,4] F9=[4,4]";
pts:=[E!0,E![-2,0],E![-1,0],E![3,0]] cat [E![-1+eps*2*i,eta*(4-eps*2*i)]:eps,eta in [-1,1]];
assert #Seqset(pts) eq 8 and Seqset(pts) eq Seqset([tm(a):a in T]);
log:=Open("C4x20/logs/C420_compact_completion_root.log","w");
for P in pts do
 if P eq E!0 or P[2] eq 0 then
  fprintf log,"C420_GLOBAL_POINT cusp=true ellipticpoint=%o\n",P;
 else
  x:=P[1];y:=P[2];tv:=(x-3)^2*(x+2)/(x+1)^2;wv:=y*(x+1)^2/(4*(x+2));
  assert tv in {-11-2*i,-11+2*i} and wv in {-2*i,2*i};
  assert (tv^2+10*tv+5)^3/tv eq 1728 and 16*(wv^2+16)^3/wv^4 eq 1728;
  fprintf log,"C420_GLOBAL_POINT cusp=false modularj=1728 ellipticpoint=%o t=%o w=%o\n",P,tv,wv;
 end if;
end for;
fprintf log,"C420_COMPACT_COMPLETION_DONE globalpoints=8 cusps=4 CM1728=4 complete=true time=%o\n",Cputime(tt);delete log;
print "C420_COMPACT_COMPLETION_DONE globalpoints=8 cusps=4 CM1728=4 complete=true",Cputime(tt);quit;
