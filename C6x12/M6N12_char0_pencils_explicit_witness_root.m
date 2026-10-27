/* Fresh relative degree4 diagonal search for6x12. All468 quartic
   source places at7 and all noncuspidal degree4 patterns at13 are included.
   Cubic+cusp and entirelycuspidal p7 branches need the separate arguments. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();SetLogFile("C6x12/logs/M6N12_char0_pencils_explicit_witness_root.log");
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");
Attach("common/fast_hecke_11.m");
Attach("C6x12/model_6_12_root.m");
Attach("common/class_rational_reduction.m");
Gs:=[* *];Hs:=[* *];hs:=[* *];Xs:=[* *];ms:=[* *];css:=[* *];c0s:=[* *];sigs0:=[];ds0:=[];
data:=eval Read("C6x12/data/M6N12_global_L8_cusp_classes_root.m");results:=[* *];
R<T>:=PolynomialRing(Rationals());F<z>:=NumberField(T^2-T+1);
 A<u,v>:=AffineSpace(F,2);
 rawmodel:=Split(Read("mdmagma/models_X1_m_n/X1_6_12.txt"));
 for line in rawmodel do
  val:=Split(Split(line,"=")[2],";")[1];
  if line[1] eq "X" then equation:=eval(val);end if;
  if line[1] eq "r" then r:=eval(val);end if;
  if line[1] eq "s" then s:=eval(val);end if;
  if line[1] eq "q" then q:=eval(val);end if;
  if line[1] eq "E" then es:=eval(val);end if;
 end for;
 C:=ProjectiveClosure(Curve(A,equation));FF:=FunctionField(C);
 X:=New(MDCrvMod11);X`M:=6;X`N:=12;
 X`_P:=[FF!q];X`_Q:=[FF!0,FF!0];X`_coordinates:=[FF.1,FF.2];
 _InitMDCrvMod(X,12,C,[FF!a:a in es]);

C:=Curve(X);FF:=FunctionField(C);print "M6N12_CHAR0_MODEL",Genus(X),Cputime(t0);
 E:=EllipticCurve([FF!a:a in X`_E]);Q:=E![0,0];assert 12*Q eq E!0;
 px:=FF!X`_P[1];ai:=aInvariants(E);R<Y>:=PolynomialRing(FF);
 ys:=Roots(Y^2+(ai[1]*px+ai[3])*Y-px^3-ai[2]*px^2-ai[4]*px-ai[5]);
 assert #ys eq 2;P:=E![px,ys[1][1]];

 if WeilPairing(P,2*Q,6) ne z then P:=-P;end if;
 assert WeilPairing(P,2*Q,6) eq z;
 print "M6N12_CHAR0_WEIL",Cputime(t0);
 fs:=[(P+5*Q)[1]-P[1],(2*P+5*Q)[1]-P[1],(P+2*Q)[1]];
 assert &and[f ne 0:f in fs];print "M6N12_CHAR0_THREE_UNITS",Cputime(t0);
 cs:=Cusps(X);groups:=AssociativeArray();
 for c in cs do
  s:=[Valuation(f,c):f in fs];
  if IsDefined(groups,s) then Append(~groups[s],c);else groups[s]:=[c];end if;
 end for;
 cols:=[34,58,27];profiles:=[[sig[j]:j in cols]:sig in data[1]];
 assert #Seqset(profiles) eq 24 and Seqset(profiles) eq Keys(groups);
 blocks:=[&+[Divisor(c):c in groups[sig]]:sig in profiles];ds:=[Degree(D):D in blocks];
 print "M6N12_CHAR0_CUSP_BLOCKS",#blocks,ds,Cputime(t0);
 sigs:=data[1];
 c0:=Support(blocks[Index(ds,1)])[1];

 assert sigs eq data[1] and ds eq data[2];
 reduced:=eval Read("C6x12/data/M6N12_reduced_L8_cusp_coefficients_magma_root.m");

 out:=[* *];
 for row in [reduced[k]:k in [1..6]] do
  index:=row[1];vec:=row[2];D:=&+[vec[i]*blocks[i]:i in [1..#blocks]];
  assert Degree(D) eq 8;V,mp:=RiemannRochSpace(D);assert Dimension(V) eq 2;
  u:=FF.1;v:=FF.2;ff:=eval Read(Sprintf("C6x12/data/M6N12_char0_pencil_%o_expression_root.txt",index));
  ps,es:=Support(Divisor(ff));poles:=&+[-es[k]*Divisor(ps[k]):k in [1..#ps]|es[k] lt 0];
  assert Degree(poles) eq 8;yes,h:=IsPrincipal(poles-D);assert yes and D+Divisor(h) eq poles;
  print "M6N12_EXPLICIT_POLE_CLASS_IDENTITY",index,h,Cputime(t0);
  ram:=Divisor(Differential(ff))+2*poles;assert IsEffective(ram) and Degree(ram) eq 32;
  rp,rm:=Support(ram);degs:=[Degree(x):x in rp];assert degs eq [8,8,8,8] and rm eq [1,1,1,1];
  generators:=[* *];
  for x in rp do
   g1,g2:=TwoGenerators(x);Append(~generators,[*Sprint(g1),Sprint(g2)*]);
   print "M6N12_EXPLICIT_RAMIFICATION_PRIME",index,Degree(x),g1,g2;
  end for;
  Append(~out,[*index,2,8,degs,rm,generators*]);
  print "M6N12_EXPLICIT_PENCIL_WITNESS",index,degs,rm,Cputime(t0);
 end for;
 Write("C6x12/data/M6N12_explicit_pencil_witness_arithmetic_root.txt",Sprint(out):Overwrite:=true);
 print "M6N12_EXPLICIT_PENCIL_WITNESS_COMPLETE",Cputime(t0);quit;
