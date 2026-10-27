/* Fresh degree8 X1(40)/F3 setup; retained continuation follows later.
   All finite cusp/Q-orbit data are recorded, including non-Q lifts. */
SetColumns(0);SetMemoryLimit(16*10^9);t40:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");
out40:=Open("C40/logs/classify_octic_40_setup.log","w");
procedure Note40(s)
 fprintf out40,"%o\n",s;Flush(out40);print s;
end procedure;
fp:=FrobeniusPolynomial(JOne(40),3);jf3:=Integers()!Evaluate(fp,1);
X:=MDX1(40,GF(3));C:=Curve(X);cs:=Cusps(X);osraw:=CuspOrbitsQ(X);
os:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in osraw];
assert Genus(X) eq 25 and jf3 eq 31579088486400;
assert #os eq 21 and #Seqset([CuspSignature(X,Support(o)[1]):o in os]) eq 21;
assert &+[Degree(c):c in cs] eq 48;
c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
Note40(Sprintf("C40_D8_MODEL_READY genus=25 finite_cusp_degrees=%o rational_orbit_degrees=%o JF3_order=%o time=%o",[Degree(c):c in cs],[Degree(o):o in os],jf3,Cputime(t40)));
SetVerbose("ClassGroup",1);
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
assert &*[n:n in Invariants(G)|n ne 0] eq jf3;
qc:=[m2(o-Degree(o)*Divisor(c0)):o in os];HQ:=sub<G|qc>;
assert Invariants(HQ) eq [2,2,2,8,120,895440] and 4592*#HQ eq jf3;
Note40(Sprintf("C40_D8_GROUP_READY invariants=%o rational_image=%o time=%o",Invariants(G),Invariants(HQ),Cputime(t40)));
places:=[NoncuspidalPlaces(X,d):d in [1..8]];
assert [#p:p in places] eq [0,0,0,32,72,80,408,708];
Note40(Sprintf("C40_D8_SETUP_DONE noncuspidal_counts=%o continuation_pending=true time=%o",[#p:p in places],Cputime(t40)));
