/* Fresh D8 necessary list: every F7-rational member of the five
   exceptional quotient pencils, lifted through the oriented A3 cover.
   Cuspidal preimages are deliberately unrestricted (safe overestimate).
   This is an exploratory certificate until independently audited. */
if not assigned C230Out then C230Out:="C2x30/logs/certified_run_v2";end if;
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "C2x30/src/a3_model.m";
// Exceptional inputs are explicitly loaded by the literal-path wrapper.
assert A3signatures eq C230ExpectedCuspSignatures;
assert [Degree(o):o in A3orbits] eq C230ExpectedCuspDegrees;
log:=Open(C230Out cat "/exceptional_fibres.log","w");
ar:=Open(C230Out cat "/exceptional_fibres_data.m","w");

// Intrinsic exact pair description, independent of ambient field bases:
// chosen linear combination, its minimal polynomial, and the two
// coordinates in its power basis. No generic finite-field Embed is used.
function PairDescription(values)
 K:=Parent(values[1]);assert Parent(values[2]) cmpeq K;
 k:=PrimeField(K);p:=#k;e:=LCM([Degree(MinimalPolynomial(v,k)):v in values]);
 found:=false;
 for bits in [[1,0],[0,1]] cat [[1,a]:a in [1..p-1]] do
  gen:=&+[K!bits[i]*values[i]:i in [1,2]];pol:=MinimalPolynomial(gen,k);
  if Degree(pol) eq e then chosen:=bits;found:=true;break;end if;
 end for;
 assert found;M:=Matrix(k,[Eltseq(gen^i,k):i in [0..e-1]]);assert Rank(M) eq e;
 cs:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in values];
 assert &and[&+[K!cs[j][i+1]*gen^i:i in [0..e-1]] eq values[j]:j in [1,2]];
 return [*e,chosen,[Integers()!a:a in Coefficients(pol)],
  [[Integers()!a:a in row]:row in cs]*];
end function;

function OrientedValues(E,P2,Q30)
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));assert a1 eq 0 and a3 eq 0;
 assert Order(P2) eq 2 and Order(Q30) eq 30 and P2 ne 15*Q30;
 P:=2*Q30;e1:=P2[1];e2:=(15*Q30)[1];e3:=-a2-e1-e2;
 assert e1 ne e2 and e2 ne e3 and e3 ne e1;
 rr:=4*(e1-e2)*(e2-e3)*(e3-e1);assert rr^2 eq Discriminant(E);
 xx:=P[1];yy:=P[2];slope:=(3*xx^2+2*a2*xx+a4)/(2*yy);
 aa:=3*xx+a2-slope^2;assert aa ne 0;u:=2*yy/aa;
 bc:=TateNormalForm_bc(E,P);ET:=EllipticCurve([1-bc[2],-bc[1],-bc[1],0,0]);
 assert Discriminant(ET)*u^12 eq Discriminant(E);
 vv:=TateNormalForm_xy(E,P) cat [rr/u^6];assert vv[3]^2 eq Discriminant(ET);
 return vv;
end function;

pointdict:=AssociativeArray();pointdata:=[* *];
procedure PointIndex(E,P2,Q30,expected,~idx,~pointdict,~pointdata)
 assert OrientedValues(E,P2,Q30) eq expected;
 assert aInvariants(E)[1] eq 0 and aInvariants(E)[3] eq 0;
 bc:=DerickxNormalForm_bc(E,P2,Q30);desc:=PairDescription(bc);key:=Sprint(desc);
 L:=Parent(bc[1]);b,c:=Explode(bc);U:=EllipticCurve([L!0,c,L!0,b*(1-b-c),L!0]);
 UP:=U![0,0];UQ:=U![b,b];assert Order(UP) eq 2 and Order(UQ) eq 30 and UP ne 15*UQ;
 assert OrientedValues(U,UP,UQ) eq expected;
 if not IsDefined(pointdict,key) then
  Append(~pointdata,desc);pointdict[key]:=#pointdata;
 end if;idx:=pointdict[key];
end procedure;

procedure Compositions3(m,ids,~out)
 for a in [0..m] do for b in [0..m-a] do
  vv:=[a,b,m-a-b];Append(~out,[[ids[j],vv[j]]:j in [1..3]|vv[j] gt 0]);
 end for;end for;
end procedure;

candidatekeys:={};candidates:=[* *];fibres:=[* *];membercount:=0;
for ex in C230ExceptionalCases do
 idx,v,hh:=Explode(ex);D:=&+[v[i]*A3orbits[i]:i in [1..#v]];
 assert Degree(D) eq 8 and IsEffective(D);
 B:=Basis(D:Reduction:=false,Simplification:="None");assert #B eq hh and hh eq 2;
 sections:=[B[1]+a*B[2]:a in GF(7)] cat [B[2]];
 for sid in [1..#sections] do
  membercount+:=1;DD:=D+Divisor(sections[sid]);assert Degree(DD) eq 8 and IsEffective(DD);
  ss,mm:=Support(DD);current:=[*[]*];ncdegree:=&+[Integers()|Degree(ss[j])*mm[j]:j in [1..#ss]|Valuation(A3K!A3j15,ss[j]) ge 0];localparts:=[* *];valid:=true;
  for j in [1..#ss] do
   pp:=ss[j];mul:=mm[j];de:=Degree(pp);
   if Valuation(A3K!A3j15,pp) lt 0 then continue;end if;
   K:=ResidueClassField(pp);av:=[Evaluate(A3K!a,pp):a in aInvariants(A3E15)];
   ET:=EllipticCurve(av);P15:=ET![Evaluate(A3K!a,pp):a in Eltseq(A3P15)];
   assert Order(P15) eq 15;expected:=[Evaluate(f,pp):f in A3coordinates];
   assert TateNormalForm_xy(ET,P15) eq expected[1..2] and expected[3]^2 eq Discriminant(ET);
   a1,a2,a3,a4,a6:=Explode(av);cs:=[K!0,a2+a1^2/4,K!0,a4+a1*a3/2,a6+a3^2/4];
   E:=EllipticCurve(cs);PP:=E![P15[1],P15[2]+(a1*P15[1]+a3)/2];assert Order(PP) eq 15;
   R<T>:=PolynomialRing(K);f:=T^3+cs[2]*T^2+cs[4]*T+cs[5];fa:=Factorization(f);
   dd:=Sort([Degree(a[1]):a in fa]);assert dd eq [1,1,1] or dd eq [3];
   localoptions:=[* *];ids:=[];
   if dd eq [1,1,1] then
    L:=K;EL:=E;PL:=PP;roots:=[a[1]:a in Roots(f)];
   else
    if mul mod 3 ne 0 then
     valid:=false;Append(~localparts,[*de,mul,dd,[],0*]);current:=[* *];continue;
    end if;
    L<root>:=ext<K|f>;EL:=BaseChange(E,L);PL:=EL![L!a:a in Eltseq(PP)];
    roots:=[a[1]:a in Roots(PolynomialRing(L)!f)];
   end if;
   assert #roots eq 3 and #Seqset(roots) eq 3;expL:=[L!a:a in expected];
   for a in [1..3] do for b in [1..3] do
    if a eq b then continue;end if;cc:=6-a-b;
    orient:=4*(roots[a]-roots[b])*(roots[b]-roots[cc])*(roots[cc]-roots[a]);
    // The completed-square universal Tate curve has normalization scale1.
    xx:=PL[1];yy:=PL[2];sl:=(3*xx^2+2*cs[2]*xx+cs[4])/(2*yy);
    u:=2*yy/(3*xx+cs[2]-sl^2);assert u eq 1;
    if orient ne expL[3] then continue;end if;
    P2:=EL![roots[a],0];R2:=EL![roots[b],0];Q30:=R2+8*PL;
    assert 2*Q30 eq PL and 15*Q30 eq R2;
    pid:=0;PointIndex(EL,P2,Q30,expL,~pid,~pointdict,~pointdata);Append(~ids,pid);
   end for;end for;
   assert #ids eq 3;
   if dd eq [1,1,1] then
    assert #Seqset(ids) eq 3 and &and[pointdata[i][1] eq de:i in ids];
    Compositions3(mul,ids,~localoptions);
   else
    assert #Seqset(ids) eq 1 and pointdata[ids[1]][1] eq 3*de;
    Append(~localoptions,[[ids[1],mul div 3]]);
   end if;
   Append(~localparts,[*de,mul,dd,ids,#localoptions*]);new:=[* *];
   for old in current do for opt in localoptions do Append(~new,old cat opt);end for;end for;
   current:=new;
  end for;
  if not valid then current:=[* *];end if;
  if ncdegree eq 0 then current:=[* *];end if;
  exported:=[];
  for cand in current do
   combined:=AssociativeArray();
   for row in cand do if IsDefined(combined,row[1]) then combined[row[1]]+:=row[2];else combined[row[1]]:=row[2];end if;end for;
   clean:=Sort([[i,combined[i]]:i in Keys(combined)]);
   assert &+[pointdata[row[1]][1]*row[2]:row in clean] eq ncdegree;
   key:=Sprint(clean);
   if key notin candidatekeys then Include(~candidatekeys,key);Append(~candidates,clean);end if;
   Append(~exported,Index(candidates,clean));
  end for;
  Append(~fibres,[*idx,sid,ncdegree,valid,localparts,exported*]);
  fprintf log,"EXCEPTIONAL_FIBRE class=%o member=%o ncdegree=%o valid=%o lifts=%o points=%o totalcandidates=%o time=%o\n",idx,sid,ncdegree,valid,#current,#pointdata,#candidates,Cputime(tt);Flush(log);
 end for;
end for;
assert membercount eq 40;
fprintf ar,"C230Points := %o;\nC230Candidates := %o;\nC230Fibres := %o;\nC230FibresComplete := true;\nC230ExclusionClaimed := false;\n",pointdata,candidates,fibres;delete ar;
fprintf log,"EXCEPTIONAL_FIBRES_DONE members=40 points=%o candidates=%o complete=true exclusion_claimed=false time=%o\n",#pointdata,#candidates,Cputime(tt);delete log;
print "EXCEPTIONAL_FIBRES_DONE",#pointdata,#candidates,Cputime(tt);quit;
