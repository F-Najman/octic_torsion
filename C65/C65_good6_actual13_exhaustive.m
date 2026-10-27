AttachSpec("mdmagma/v2/mdmagma.spec");
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(3*10^9);st65:=Cputime();
X13:=MDX1(13,GF(2));C13:=Curve(X13);F13:=FunctionField(C13);
EE:=EllipticCurve(X13`_E);PP:=EE!X13`_P;
xy3:=TateNormalForm_xy(EE,3*PP);xy9:=TateNormalForm_xy(EE,9*PP);
phi3:=iso<C13->C13|[xy3[1],xy3[2],F13!1],[xy9[1],xy9[2],F13!1]:Check:=false>;
GG:=AutomorphismGroup(C13,[phi3]);assert #GG eq 3;
Q,pi:=CurveQuotient(GG);assert Genus(Q) eq 0;
function ZProjective65(E,R13)
 K:=BaseRing(E);xy:=TateNormalForm_xy(E,R13);pt:=C13(K)![xy[1],xy[2],K!1];qq:=pi(pt);
 if qq[1] eq 0 then assert qq[2] eq 0;return true,K!0;end if;
 zz:=qq[2]/qq[1]+1;assert zz^4+zz+1 eq 0;return false,zz;
end function;

function BCKey65(bc)
 L:=Parent(bc[1]);d:=LCM([Degree(MinimalPolynomial(v,GF(2))):v in bc]);
 found:=false;
 for bits in [[1,0],[0,1],[1,1]] do
  gen:=bits[1]*bc[1]+bits[2]*bc[2];mp:=MinimalPolynomial(gen,GF(2));
  if Degree(mp) eq d then found:=true;bb:=bits;break;end if;
 end for;assert found;
 MB:=Matrix(GF(2),[Eltseq(gen^i,GF(2)):i in [0..d-1]]);assert Rank(MB) eq d;
 rows:=[Eltseq(Solution(MB,Vector(GF(2),Eltseq(v,GF(2))))):v in bc];
 return [d] cat bb cat [Integers()!v:v in Eltseq(mp)] cat [Integers()!v:v in &cat rows];
end function;

FF<a>:=GF(64);curves65:=[];for jj in FF do for E0 in Twists(EllipticCurveFromjInvariant(jj)) do if #E0 eq 65 then Append(~curves65,E0);end if;end for;end for;assert #curves65 eq 1 and jInvariant(curves65[1]) eq 0;
E:=curves65[1];AG,ag:=AbelianGroup(E);assert Invariants(AG) eq [65];PN:=ag(AG.1);
maps65:=MDIsogenies(E,3);assert #maps65 eq 4;
total65:=0;infty65:=0;outfinite65:=0;outinf65:=0;degrees65:={};
for k in [k:k in [1..64]|GCD(k,65) eq 1] do
 R65:=k*PN;R13:=5*R65;assert BCKey65(TateNormalForm_bc(E,R65))[1] eq 6;imagekeys65:=[];inf,zz:=ZProjective65(E,R13);assert inf;infty65+:=1;
 for ph in maps65 do
  Et:=Codomain(ph);Rt:=ph(R65);inf,zz:=ZProjective65(Et,5*Rt);
  if inf then outinf65+:=1;else outfinite65+:=1;assert Degree(MinimalPolynomial(zz,GF(2))) eq 4;end if;
  bc:=TateNormalForm_bc(Et,Rt);dd:=LCM([Degree(MinimalPolynomial(v,GF(2))):v in bc]);Include(~degrees65,dd);assert dd eq 12;Append(~imagekeys65,BCKey65(bc));
  total65+:=1;
 end for;
 assert #Seqset(imagekeys65) eq 2 and &and[#[v:v in imagekeys65|v eq w] eq 2:w in imagekeys65];
end for;
printf "C65_GOOD6_ACTUAL13_CLOSED_COVERAGE_DONE generators=%o input_z_infinity=%o isogeny_images=%o finite_z=%o infinity_z=%o full_source_degrees=%o time=%o\n",infty65,infty65,total65,outfinite65,outinf65,degrees65,Cputime(st65);
quit;
