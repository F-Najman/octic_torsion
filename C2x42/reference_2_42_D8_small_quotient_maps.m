SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "common/pair_place_linear.m";
Y:=MDX11(2,14,GF(5));Z:=MDX1(21,GF(5));R:=PolynomialRing(GF(5));
for D in [4..8] do
 a:=eval Read(Sprintf("C2x42/data/M2N42_D8_direct_full%o_arithmetic.txt",D));
 K<z>:=ext<GF(5)|R!a[2]>;row:=[row:row in a[9]|row[1] eq D][1];key:=row[2];
 b:=Evaluate(R!key[1..D],z);c:=Evaluate(R!key[D+1..2*D],z);
 E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];
 assert 42*Q eq E!0 and 21*Q ne P and 21*Q ne E!0 and 14*Q ne E!0 and 6*Q ne E!0;
 A:=3*Q;bc:=[A[1]^3/A[2]^2,c*A[1]^2/A[2]^2];
 y:=C48PairPlaceLinear(Y`_coordinates,bc);
 yr:=ModuliPoint(Y,E,rec<recformat<P:PtEll,Q:PtEll>|P:=P,Q:=A>);
 assert y eq yr;xy:=TateNormalForm_xy(E,2*Q);
 zz:=C48PairPlaceLinear([FunctionField(Curve(Z)).1,FunctionField(Curve(Z)).2],xy);
 zr:=ModuliPoint(Z,E,rec<recformat<P:PtEll>|P:=2*Q>);
 assert zz eq zr;print "M2N42_NATIVE_SMALL_MAP_REFERENCES",D,Degree(y),Degree(zz),Cputime(t0);
end for;
print "M2N42_NATIVE_SMALL_MAP_REFERENCES_DONE";quit;
