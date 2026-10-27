/* Exact forward map from an already PROVEN native class-group presentation.
   A rough degree-zero D is replaced by E-R-div(b), where deg R=g,
   b in L(D+R), and E is effective. All identities are checked as finite
   AND infinite fractional ideals, before any unknown-prime decomposition.
   The native abstract-generator basis uses REVERSED factorbase order. */
function MDClassGroupDirectRRCheckedV4(F0,certpath)
 F:=RationalExtensionRepresentation(F0);
 assert assigned F`cld and F`cld`proven and assigned F`cld`Cl;
 assert assigned F`cld`ZnToCl and assigned F`cld`ClToDivHigh;
 assert assigned F`cld`A and Degree(F`cld`A) eq 1;
 Cl:=F`cld`Cl;ff:=F`cld`ZnToCl;Zn:=Domain(ff);A:=F`cld`A;
 fb:=Reverse([v`place:v in F`cld`fb]);n:=#fb;
 assert Ngens(Zn) eq n;
 OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
 fi:=[*Ideal(P):P in fb*];fd:=[Degree(P):P in fb];bound:=Maximum(fd);
 isfin:=[Order(J) eq OF:J in fi];
 assert &and[isfin[j] or Order(fi[j]) eq OI:j in [1..n]];
 cert:=Open(certpath,"w");
 witnesses:=Open(certpath cat ".witnesses.py","w");
 fprintf witnesses,"witnesses = [\n";Flush(witnesses);
 function SparseCoordinates(v)
  e:=Eltseq(v);return [[j,e[j]]:j in [1..#e]|e[j] ne 0];
 end function;
 function SerializeFunction(a)
  assert #ConstantField(F) eq Characteristic(F);
  return [[ [Integers()!x:x in Coefficients(Numerator(z))],
             [Integers()!x:x in Coefficients(Denominator(z))] ]:z in Eltseq(a)];
 end function;
 function DeserializeFunction(a)
  B:=BaseField(F);P:=PolynomialRing(ConstantField(F));
  return F![B!(P!z[1])/B!(P!z[2]):z in a];
 end function;
 function SmoothSupport(D)
  ps,es:=Support(D);v:=Zn!0;
  for i in [1..#ps] do
   j:=Index(fb,ps[i]);if j eq 0 then return false,Zn!0;end if;
   v+:=es[i]*Zn.j;
  end for;
  return true,v;
 end function;
 oka,va:=SmoothSupport(A);assert oka;
 function SmoothIdeals(IE,JE)
  assert IsIntegral(IE) and IsIntegral(JE);
  mm:=Minimum(IE);fs:=Factorisation(mm);
  if &or[Degree(x[1]) gt bound:x in fs] then return false,Zn!0;end if;
  v:=Zn!0;Ileft:=IE;Jleft:=JE;
  for j in [1..n] do
   if isfin[j] then
    if not IsDivisibleBy(mm,Minimum(fi[j])) then continue;end if;
    e:=Valuation(IE,fi[j]);assert e ge 0;
    if e gt 0 then v+:=e*Zn.j;Ileft/:=fi[j]^e;end if;
   else
    e:=Valuation(JE,fi[j]);assert e ge 0;
    if e gt 0 then v+:=e*Zn.j;Jleft/:=fi[j]^e;end if;
   end if;
  end for;
  return IsOne(Ileft) and IsOne(Jleft),v;
 end function;
 function MapClass(D)
  ok,v:=SmoothSupport(D);
  if ok then return ff(v);end if;
  d:=Degree(D);D0:=D-d*A;g:=Genus(F);
  ID,JD:=Ideals(D0);
  for trial in [1..1000] do
   left:=g;rv:=Zn!0;R:=DivisorGroup(F)!0;
   while left gt 0 do
    j:=Random([i:i in [1..n]|fd[i] le left]);
    rv+:=Zn.j;R+:=Divisor(fb[j]);left-:=fd[j];
   end while;
   assert Degree(R) eq g;
   T:=D0+R;assert Degree(T) eq g;
   bb:=Basis(T:Reduction:=false,Simplification:="None");assert #bb gt 0;
   b:=bb[1];IR,JR:=Ideals(R);IT,JT:=Ideals(T);
   BF:=b*OF;BI:=b*OI;IE:=IT*BF;JE:=JT*BI;
   assert IsIntegral(IE) and IsIntegral(JE);
   // Degree(E)=g follows from E=T+div(b), degree(T)=g.
   // Avoid the map-degree convention of Degree(Norm(infinite_ideal)).
   // This is the full divisor identity D0=E-R-div(b), with no support
   // factorization of the principal divisor and with both orders checked.
   assert IE/(IR*BF) eq ID and JE/(JR*BI) eq JD;
   smooth,ev:=SmoothIdeals(IE,JE);
   if smooth then
    // Check the actual factorbase representation of E in both orders.
    EF:=ideal<OF|1>;EI:=ideal<OI|1>;ec:=Eltseq(ev);
    for j in [1..n] do
     if ec[j] ne 0 then
      if isfin[j] then EF*:=fi[j]^ec[j];else EI*:=fi[j]^ec[j];end if;
     end if;
    end for;
    assert EF eq IE and EI eq JE;
    assert &+[ec[j]*fd[j]:j in [1..n]] eq g;
    result:=ff(ev-rv+d*va);
    // Preserve exact proof objects, not only successful assertion markers.
    bdata:=SerializeFunction(b);assert DeserializeFunction(bdata) eq b;
    fprintf witnesses,"[%o,%o,%o,%o,%o,%o,%o],\n",trial,d,Eltseq(result),SparseCoordinates(ev),SparseCoordinates(rv),SparseCoordinates(va),bdata;Flush(witnesses);
    fprintf cert,"DIRECT_RR_CLASS_IDENTITY trials=%o degree=%o class=%o Edegree=%o Rdegree=%o both_orders_exact=true\n",trial,d,Eltseq(result),g,g;Flush(cert);
    return result;
   end if;
   if trial mod 10 eq 0 then
    fprintf cert,"DIRECT_RR_RETRY trials=%o both_orders_exact=true\n",trial;Flush(cert);
   end if;
  end for;
  error "Direct exact RR class descent exhausted its honest1000-trial bound";
 end function;
 high:=F`cld`ClToDivHigh;
 phi:=hom<DivisorGroup(F)->Cl|D:->MapClass(D),a:->high(a)>;
 if F0 cmpeq F then return Cl,Inverse(phi),phi,witnesses;end if;
 psi:=hom<DivisorGroup(F0)->Cl|D:->phi(DivisorGroup(F)!D),a:->DivisorGroup(F0)!(a@@phi)>;
 return Cl,Inverse(psi),psi,witnesses;
end function;
