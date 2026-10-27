// helpers8.m -- shared helpers for the independent degree-8 classification checks.
// Written for this review (adapted from ~/degree6_torsion/scratch/degree7_review/helpers7.m).
// Models: Sutherland's raw X_1(N) models FFFcN.txt (van Hoeij's (x,y)-coordinates), read from the mdmagma copy.
MODELDIR := "mdmagma/models_X1_n/";

// Tate normal form parameters (b,c) of (E,Q), Q of order >= 4, by direct substitution.
function tate(E,Q)
  k:=BaseRing(E); a1,a2,a3,a4,a6:=Explode(aInvariants(E)); x0:=Q[1]; y0:=Q[2];
  R<X,Y>:=PolynomialRing(k,2); f:=Y^2+a1*X*Y+a3*Y-(X^3+a2*X^2+a4*X+a6);
  g:=Evaluate(f,[X+x0,Y+y0]); assert MonomialCoefficient(g,R!1) eq 0;
  A3:=MonomialCoefficient(g,Y); A4:=-MonomialCoefficient(g,X); assert A3 ne 0;
  g:=Evaluate(g,[X,Y+(A4/A3)*X]); assert MonomialCoefficient(g,X) eq 0;
  A1:=MonomialCoefficient(g,X*Y); A2:=-MonomialCoefficient(g,X^2); A3:=MonomialCoefficient(g,Y);
  assert MonomialCoefficient(g,Y^2) eq 1 and MonomialCoefficient(g,X^3) eq -1;
  return -A2^3/A3^2, 1-A1*A2/A3;
end function;

// raw (x,y) -> (b,c), van Hoeij / Sutherland convention
function xy2bc(xx,yy)
  rr:=(xx^2*yy - xx*yy + yy - 1)/(xx^2*yy - xx); ss:=(xx*yy - yy + 1)/(xx*yy);
  return rr*ss*(rr-1), ss*(rr-1);
end function;

TateCurve := func<b,c | EllipticCurve([1-c,-b,-b,0,0])>;

// Canonical identity of the closed point of Y_1(N)_{F_p} with Tate coordinates (bb,cc) (any finite field containing them).
// Depends only on the Frobenius orbit of (bb,cc); also returns its degree [F_p(bb,cc):F_p].
function ptkey(bb, cc)
  k:=Parent(bb); F0:=PrimeField(k); cc:=k!cc;
  db:=Degree(MinimalPolynomial(bb,F0)); dc:=Degree(MinimalPolynomial(cc,F0)); d:=LCM(db,dc);
  cands:=[bb, cc] cat [bb + l*cc : l in F0 | l ne 0];
  for i in [1..#cands] do g:=cands[i];
    if Degree(MinimalPolynomial(g,F0)) eq d then
      M:=Matrix(F0,[Eltseq(g^e,F0): e in [0..d-1]]);
      ok1,s1:=IsConsistent(M,Vector(F0,Eltseq(bb,F0))); ok2,s2:=IsConsistent(M,Vector(F0,Eltseq(cc,F0)));
      assert ok1 and ok2;
      return <d, i, Eltseq(MinimalPolynomial(g,F0)), Eltseq(s1), Eltseq(s2)>, d;
    end if; end for;
  error "ptkey: no generator";
end function;

// key of a noncuspidal place of the function field, given b,c in the function field
function placekey(Q, b, c) k,mk:=ResidueClassField(Q); key,d:=ptkey(mk(b),mk(c)); assert d eq Degree(Q); return key; end function;

// all effective divisors of degree rem supported on the places cusps[i..]
procedure effdivs(cusps, i, rem, cur, ~L)
  if rem eq 0 then Append(~L,cur); return; end if; if i gt #cusps then return; end if;
  for k in [0..(rem div Degree(cusps[i]))] do effdivs(cusps, i+1, rem-k*Degree(cusps[i]), cur+k*cusps[i], ~L); end for;
end procedure;
// same, but as multiplicity vectors (sequences of integers indexed like cusps)
procedure effvecs(degs, i, rem, cur, ~L)
  if rem eq 0 then Append(~L,cur); return; end if; if i gt #degs then return; end if;
  for k in [0..(rem div degs[i])] do c2:=cur; c2[i]:=k; effvecs(degs, i+1, rem-k*degs[i], c2, ~L); end for;
end procedure;

// Specialisation at p of the closed point given by Tate coordinates (bK,cK) in a number field K.
// Returns a list of <key, multiplicity, degree> (multiplicity e*f/d at each prime above p), and
// a flag "cuspidal" if some prime above p gives bad reduction.
function specialise(K, bK, cK, N, p : Diamond:=1)
  EK:=TateCurve(bK,cK); PK:=EK![0,0]; assert Order(PK) eq N;
  if Diamond ne 1 then bK,cK:=tate(EK, Diamond*PK); end if;
  O:=MaximalOrder(K); out:=[]; cusp:=false;
  for pp in Factorization(p*O) do e:=pp[2]; kk,red:=ResidueClassField(pp[1]); f:=Degree(kk);
    if Valuation(bK,pp[1]) lt 0 or Valuation(cK,pp[1]) lt 0 then cusp:=true; continue; end if;
    b1:=red(bK); c1:=red(cK); Ek:=TateCurve(b1,c1);
    if Discriminant(Ek) eq 0 then cusp:=true; continue; end if;
    assert Order(Ek![0,0]) eq N;
    key,d:=ptkey(b1,c1); Append(~out, <key, e*f div d, d>);
  end for;
  return out, cusp;
end function;

// Sutherland's raw model of X_1(N) over the finite field k, as a function field; returns F, b, c.
function rawmodel(N, k)
  s:=Read(MODELDIR cat "FFFc" cat IntegerToString(N) cat ".txt");
  R<x,y>:=PolynomialRing(Integers(),2); f:=eval s;
  Fx<T>:=RationalFunctionField(k); PY<Y>:=PolynomialRing(Fx);
  if Degree(f,y) le Degree(f,x) then
    g:=Evaluate(f,[PY!T, Y]); assert IsIrreducible(g); F:=FunctionField(g); xx:=F!T; yy:=F.1;
  else
    g:=Evaluate(f,[Y, PY!T]); assert IsIrreducible(g); F:=FunctionField(g); xx:=F.1; yy:=F!T;
  end if;
  b,c:=xy2bc(xx,yy); return F,b,c;
end function;

// multiset of degrees of the Gal(Qbar/Q)-orbits of cusps of X_1(N) (P a point of order N):
// cusp <-> (bb mod N up to sign, a in (Z/g)^*), g=gcd(bb,N); field Q(zeta_g) (or its real subfield if 2bb = 0 mod N).
function cusporbitdegrees(N)
  out:=[];
  for g in Divisors(N) do
    if N div g le 2 then Append(~out, g le 2 select 1 else EulerPhi(g) div 2);
    else out cat:= [EulerPhi(g) : i in [1..EulerPhi(N div g) div 2]]; end if;
  end for;
  return Sort(out);
end function;

// index of Gamma_1(N).{+-1} in PSL_2(Z) = degree of j on X_1(N)
function jdegree(N) d:=N^2/2; for pr in PrimeDivisors(N) do d*:=(1-1/pr^2); end for; return Integers()!d; end function;

// #J_1(N)(F_p) = prod over even Dirichlet characters eps mod N of det(1 - T_p + p*eps(p)) on S_2(N,eps)
// (Eichler--Shimura), via sign +1 cuspidal modular symbols, one character per Galois orbit, norms to Q.
function JacOrderMS(N, p)
  G:=DirichletGroup(N, CyclotomicField(EulerPhi(N))); tot:=1; dim:=0;
  for eps0 in GaloisConjugacyRepresentatives(G) do
    if Evaluate(eps0,-1) ne 1 then continue; end if; eps:=MinimalBaseRingCharacter(eps0);
    M:=CuspidalSubspace(ModularSymbols(eps,2,+1)); if Dimension(M) eq 0 then continue; end if;
    Tp:=HeckeOperator(M,p); K:=BaseRing(Tp); n:=Nrows(Tp);
    d:=Determinant(IdentityMatrix(K,n) - Tp + p*(K!Evaluate(eps,p))*IdentityMatrix(K,n));
    tot*:=Norm(d); dim+:=n*Degree(K);
  end for;
  return Integers()!tot, dim;
end function;
