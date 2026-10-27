load "common/pair_af_place_linear.m";
/* Exact extension of the checked two-coordinate AF linear inverse.
   Input three generators of a fine affine moduli chart; assert unique place. */
function MDPlaceAFThreeViaBase(coords,values,basecoords)
 assert #coords eq 3 and #values eq 3;
 K:=Parent(values[1]);assert &and[Parent(a) cmpeq K:a in values];
 k:=PrimeField(K);p:=#k;e:=LCM([Degree(MinimalPolynomial(a,k)):a in values]);
 bitses:=[[1,0,0],[0,1,0],[0,0,1]] cat [[1,a,b]:a,b in [0..p-1]] cat [[0,1,a]:a in [0..p-1]];
 found:=false;
 for bits in bitses do
  gen:=&+[K!bits[i]*values[i]:i in [1..3]];pol:=MinimalPolynomial(gen,k);
  if Degree(pol) eq e then found:=true;chosen:=bits;break;end if;
 end for;assert found;
 M:=Matrix(k,[Eltseq(gen^i,k):i in [0..e-1]]);assert Rank(M) eq e;
 cs:=[Eltseq(Solution(M,Vector(k,Eltseq(a,k)))):a in values];
 assert &and[&+[K!cs[j][i+1]*gen^i:i in [0..e-1]] eq values[j]:j in [1..3]];
 F:=Parent(coords[1]);assert &and[Parent(a) cmpeq F:a in coords];
 gg:=&+[F!chosen[i]*coords[i]:i in [1..3]];
 fs:=[Evaluate(pol,gg)] cat [coords[j]-&+[F!cs[j][i+1]*gg^i:i in [0..e-1]]:j in [1..3]];
 bp:=C48PairAFPlaceLinear(basecoords,values[1..2]);aa,bb:=TwoGenerators(bp);
 raw:=Setseq(Seqset(Decomposition(F,bp)));
 ps:=[q:q in raw|Valuation(F!aa,q) gt 0 and Valuation(F!bb,q) gt 0];
 ps:=[q:q in ps|&and[f eq 0 or Valuation(f,q) gt 0:f in fs]];
 assert #ps eq 1 and Degree(ps[1]) eq e;
 assert &and[f eq 0 or Valuation(f,ps[1]) gt 0:f in fs];
 return ps[1];
end function;
