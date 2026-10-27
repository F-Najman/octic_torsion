/* Q-defined degree2 modular isogeny map, with exact order and independence guards. */
M248:=recformat<P:PtEll,Q:PtEll>;
function C48TwoIsogenousModuli(E,P)
 assert 48*P eq E!0 and 24*P ne E!0 and 16*P ne E!0;
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));T:=24*P;r:=T[1];
 A:=3*r+a2+a1^2/4;
 B:=3*r^2+2*r*(a2+a1^2/4)+a4+a1*a3/2;
 assert T[2]+(a1*r+a3)/2 eq 0;
 assert r^3+(a2+a1^2/4)*r^2+(a4+a1*a3/2)*r+a6+a3^2/4 eq 0;
 xx:=P[1]-r;yy:=P[2]+(a1*P[1]+a3)/2;
 assert xx ne 0 and yy^2 eq xx^3+A*xx^2+B*xx;
 F:=BaseRing(E);E2:=EllipticCurve([F!0,-2*A,F!0,A^2-4*B,F!0]);
 Q:=E2![yy^2/xx^2,yy*(B-xx^2)/xx^2];T2:=E2![0,0];
 assert 24*Q eq E2!0 and 12*Q ne E2!0 and 8*Q ne E2!0;
 assert 2*T2 eq E2!0 and T2 ne E2!0 and T2 ne 12*Q;
 assert aInvariants(E2)[1] eq 0 and aInvariants(E2)[3] eq 0;
 return E2,rec<M248|P:=T2,Q:=Q>;
end function;
