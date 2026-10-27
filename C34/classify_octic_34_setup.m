/* Fresh degree8 X1(34)/F3 setup. It retains the exact class group for a
   separately frozen continuation; no classification is asserted here. */
SetColumns(0);SetMemoryLimit(16*10^9);t34:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");
out34:=Open("C34/logs/classify_octic_34_setup.log","w");
procedure Note34(s)
 fprintf out34,"%o\n",s;Flush(out34);print s;
end procedure;
fp:=FrobeniusPolynomial(JOne(34),3);jf3:=Integers()!Evaluate(fp,1);
X:=MDX1(34,GF(3));C:=Curve(X);cs:=Cusps(X);osraw:=CuspOrbitsQ(X);
os:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in osraw];
assert Genus(X) eq 21 and jf3 eq 667924070400;
assert {*Degree(c):c in cs*} eq {*1^^16,8^^2*};
assert {*Degree(o):o in os*} eq {*1^^16,8^^2*};
assert #Seqset([CuspSignature(X,Support(o)[1]):o in os]) eq 18;
assert &and[#Support(o) eq 1 and Degree(Support(o)[1]) eq Degree(o):o in os];
c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
Note34(Sprintf("C34_D8_MODEL_READY genus=21 JF3_order=%o all_cusps_Q_liftable=true time=%o",jf3,Cputime(t34)));
SetVerbose("ClassGroup",1);
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
assert &*[n:n in Invariants(G)|n ne 0] eq jf3;
qc:=[m2(o-Degree(o)*Divisor(c0)):o in os];HQ:=sub<G|qc>;
assert Invariants(HQ) eq [8760,595680] and 128*#HQ eq jf3;
Note34(Sprintf("C34_D8_GROUP_READY invariants=%o rational_image=%o time=%o",Invariants(G),Invariants(HQ),Cputime(t34)));
places:=[NoncuspidalPlaces(X,d):d in [1..8]];
assert [#p:p in places] eq [0,0,0,26,32,80,256,734];
Note34(Sprintf("C34_D8_SETUP_DONE noncuspidal_counts=%o continuation_pending=true time=%o",[#p:p in places],Cputime(t34)));
