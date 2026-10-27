// formal_common.m -- shared routines for the independent degree-8 review of the
// formal-immersion (Cartier) certificates at 2.  Written from scratch for this review;
// nothing here is taken from the Codex scripts.
//
// Conventions.  F is a global function field over a finite field K, g its genus,
// V, h := DifferentialSpace(DivisorGroup(F)!0) a basis of H^0(Omega).  Everything is
// expressed in coordinates with respect to this one basis.
// Magma's DifferentialSpace(D) is {omega : div(omega) >= D}; so for a place P,
// DifferentialSpace(m*P) is the kernel of restriction of H^0(Omega) to the jets mod m*P.

// Magma -b buffers stdout, so every progress line is also appended to a progress file
// (Write opens and closes the file each time, so it is visible immediately).
procedure Say(pf, s)
  print s;
  if pf ne "" then Write(pf, s); end if;
end procedure;

function DedupPlaces(L)
  out := [];
  for P in L do
    if not P in out then Append(~out, P); end if;
  end for;
  return out;
end function;

// coefficient space of DifferentialSpace(m*P) inside K^g
function KernelCoords(P, m, h, K, g)
  Vm, hm := DifferentialSpace(m*P);
  W := sub< VectorSpace(K, g) | [Vector(K, Eltseq(hm(Vm.k) @@ h)) : k in [1..Dimension(Vm)]] >;
  assert Dimension(W) eq Dimension(Vm);
  return W;
end function;

// all kernel spaces W[<i,m>] for i in idx, 1 <= m <= mmax(i); checks nesting and the
// trivial lower bound dim >= g - m deg P.  Prints a greppable line per space.
function AllKernels(cusps, idx, mmax, h, K, g, t0, pf)
  W := AssociativeArray();
  for i in idx do
    P := cusps[i];
    for m in [1..mmax[i]] do
      t1 := Cputime();
      W[<i,m>] := KernelCoords(P, m, h, K, g);
      dm := Dimension(W[<i,m>]);
      assert dm ge g - m*Degree(P);
      if m gt 1 then assert W[<i,m>] subset W[<i,m-1>]; end if;
      Say(pf, Sprintf("KERNEL cusp %o deg %o m %o dim %o (codim %o) [%o s, total %o s]",
             i, Degree(P), m, dm, g - dm, Cputime(t1), Cputime(t0)));
    end for;
  end for;
  return W;
end function;

// same, for an explicit increasing list mlist[i] of multiplicities at cusp i
function AllKernelsList(cusps, idx, mlist, h, K, g, t0, pf)
  W := AssociativeArray();
  for i in idx do
    P := cusps[i]; prev := 0;
    for m in mlist[i] do
      t1 := Cputime();
      W[<i,m>] := KernelCoords(P, m, h, K, g);
      dm := Dimension(W[<i,m>]);
      assert dm ge g - m*Degree(P);
      if prev gt 0 then assert W[<i,m>] subset W[<i,prev>]; end if;
      prev := m;
      Say(pf, Sprintf("KERNEL cusp %o deg %o m %o dim %o (codim %o) [%o s, total %o s]",
             i, Degree(P), m, dm, g - dm, Cputime(t1), Cputime(t0)));
    end for;
  end for;
  return W;
end function;

// Cartier matrix in the basis h(V.1..V.g): row k = coordinates of C(omega_k).
// Over F_2 the Cartier operator is linear, so U = rowspace(MC - I) = (C-1)H^0(Omega).
function CartierU(V, h, g)
  K := BaseRing(V);
  B := [h(V.i) : i in [1..g]];
  MC := Matrix(K, g, g, &cat[Eltseq(Cartier(B[k]) @@ h) : k in [1..g]]);
  I := IdentityMatrix(K, g);
  return RowSpace(MC - I), MC;
end function;

// number of effective divisors of degree d on closed points of the given degrees
function CountDivisors(degs, d)
  R<x> := PowerSeriesRing(Integers(), d+1);
  s := &*[ R | 1/(1 - x^e) : e in degs | e le d ];
  return Coefficient(s, d);
end function;

// Depth-first enumeration of all effective divisors sum_j m_j P_{idx[j]} of degree rem
// (weights degs).  Kc is the running intersection (inside the ambient space) of the
// spaces WU[<i,m>].  At each leaf the rank dimU - dim Kc is compared with target.
procedure EnumRanks(j, rem, Kc, idx, degs, WU, dimU, target, cur, ~cnt, ~bad, ~minr, ~maxr, pf)
  if rem eq 0 then
    cnt +:= 1;
    r := dimU - Dimension(Kc);
    if r lt minr then minr := r; end if;
    if r gt maxr then maxr := r; end if;
    if r ne target then
      bad +:= 1;
      if bad le 30 then Say(pf, Sprintf("FAILURE divisor %o rank %o", cur, r)); end if;
    end if;
    if cnt mod 100000 eq 0 then Say(pf, Sprintf("PROGRESS enumerated %o failures %o", cnt, bad)); end if;
    return;
  end if;
  if j gt #idx then return; end if;
  i := idx[j];
  EnumRanks(j+1, rem, Kc, idx, degs, WU, dimU, target, cur, ~cnt, ~bad, ~minr, ~maxr, pf);
  for m in [1..(rem div degs[i])] do
    EnumRanks(j+1, rem - m*degs[i], Kc meet WU[<i,m>], idx, degs, WU, dimU, target,
              Append(cur, <i,m>), ~cnt, ~bad, ~minr, ~maxr, pf);
  end for;
end procedure;

// numerator of x([n]Q) as a polynomial in X = x(Q), for odd n, over the base of E.
// Uses Magma's convention (checked in check_divpoly.m): for even k, DivisionPolynomial(E,k)
// is psi_k * psi_2, so x([n]Q) = X - f_{n-1} f_{n+1} / (psi_2^2 f_n^2).
function XMultNumerator(E, n)
  assert IsOdd(n);
  f := [DivisionPolynomial(E, k) : k in [n-1, n, n+1]];
  R := Parent(f[1]); X := R.1;
  b2, b4, b6 := Explode(bInvariants(E)[1..3]);
  ps2 := 4*X^3 + b2*X^2 + 2*b4*X + b6;
  num := X*ps2*f[2]^2 - f[1]*f[3];
  q, r := Quotrem(num, ps2);
  assert r eq 0;
  return q;
end function;

// ---------------------------------------------------------------------------------------
// Kernel spaces by successive residues (no Riemann-Roch spaces).  Added 2026-10-27 because
// DifferentialSpace(m*P) costs ~2300-3800 s per (P, odd m) on the degree-49 model of X1(49).
// For a place P with residue field k (over the constant field K), choose w0 among the basis
// differentials with v_P(w0) = 0 and a local uniformizer pi.  With f_a = (sum a_k w_k)/w0,
//   S_0 = K^g,   S_{n+1} = { a in S_n : residue at P of f_a / pi^n is 0 },
// so S_m = { a : v_P(sum a_k w_k) >= m } = coefficient space of DifferentialSpace(m*P).
// Returns the sequence [S_1, ..., S_mmax] (subspaces of K^g) and the number of residue
// evaluations.  NOT USED in the end: Evaluate(f, P) turned out to be very slow at infinite
// places (minutes per call); KernelsByExpansion below is used instead.
function KernelsByResidues(P, mmax, B, K, g)
  k := ResidueClassField(P);
  pi := LocalUniformizer(P);
  assert Valuation(pi, P) eq 1;
  i0 := 0;
  for i in [1..g] do
    if Valuation(B[i], P) eq 0 then i0 := i; break; end if;
  end for;
  assert i0 gt 0;
  w0 := B[i0];
  Vg := VectorSpace(K, g);
  Sb := [Vg.i : i in [1..g]];            // basis of S_n
  fn := [B[i]/w0 : i in [1..g]];         // f_a / pi^n for a in Sb
  e := Degree(k, K);
  Ws := [* *];                            // a list: sequences of subspaces do not coerce
  nev := 0;
  for n in [0..mmax-1] do
    rows := [];
    for f in fn do
      r := Evaluate(f, P); nev +:= 1;
      Append(~rows, Eltseq(k!r, K));
    end for;
    M := Matrix(K, #fn, e, &cat rows);
    N := Kernel(M);                       // row vectors c with c*M = 0
    newSb := []; newfn := [];
    for c in Basis(N) do
      cs := Eltseq(c);
      Append(~newSb, &+[cs[i]*Sb[i] : i in [1..#Sb] | cs[i] ne 0]);
      Append(~newfn, (&+[cs[i]*fn[i] : i in [1..#Sb] | cs[i] ne 0]) / pi);
    end for;
    Sb := newSb; fn := newfn;
    Append(~Ws, sub< Vg | Sb >);
    assert Dimension(Ws[#Ws]) eq #Sb;
  end for;
  return Ws, nev;
end function;

// all kernel spaces for cusps idx (mmax[i] each) by residues; one progress line per cusp;
// each space is also written to the file kf (one line per space: i m dim, then bit rows).
function AllKernelsByResidues(cusps, idx, mmax, h, K, g, t0, pf, kf)
  V := Domain(h);
  B := [h(V.i) : i in [1..g]];
  W := AssociativeArray();
  for i in idx do
    t1 := Cputime();
    Ws, nev := KernelsByResidues(cusps[i], mmax[i], B, K, g);
    for m in [1..mmax[i]] do
      W[<i,m>] := Ws[m];
      assert Dimension(Ws[m]) ge g - m*Degree(cusps[i]);
      if kf ne "" then
        Write(kf, Sprintf("SPACE %o %o %o", i, m, Dimension(Ws[m])));
        for v in Basis(Ws[m]) do Write(kf, &cat[Sprint(x) : x in Eltseq(v)]); end for;
      end if;
    end for;
    Say(pf, Sprintf("RESKERNEL cusp %o deg %o mmax %o codims %o residues %o [%o s, total %o s]",
        i, Degree(cusps[i]), mmax[i], [g - Dimension(Ws[m]) : m in [1..mmax[i]]], nev, Cputime(t1), Cputime(t0)));
  end for;
  return W;
end function;

// ---------------------------------------------------------------------------------------
// Kernel spaces from local expansions (used for X1(49) and X1(51); validated against
// DifferentialSpace(m*P) on all 124 resp. 103 spaces of X1(45) and X1(63)/<8> in
// validate_residue_kernels.m, and on two spaces of each of X1(49), X1(51) themselves).
// Choose w0 among the basis differentials with v_P(w0) = 0; then f_i = w_i/w0 is regular at P
// and v_P(sum a_i w_i) = v_P(sum a_i f_i).  With Expand(f_i, P : AbsPrec := mmax) (Laurent series
// in Magma's local uniformizer, coefficients in the residue field k), the coefficient space of
// DifferentialSpace(m*P) is the kernel of the g x (m deg P) matrix of the coefficients of
// pi^0, ..., pi^(m-1), written over K.
function KernelsByExpansion(P, mmax, B, K, g)
  k := ResidueClassField(P);
  e := Degree(k, K);
  i0 := 0;
  for i in [1..g] do
    if Valuation(B[i], P) eq 0 then i0 := i; break; end if;
  end for;
  assert i0 gt 0;
  w0 := B[i0];
  ex := [Expand(B[i]/w0, P : AbsPrec := mmax) : i in [1..g]];
  assert &and[ Valuation(sr) ge 0 : sr in ex ];
  rows := [ &cat[ Eltseq(k!Coefficient(sr, j), K) : j in [0..mmax-1] ] : sr in ex ];
  M := Matrix(K, g, e*mmax, &cat rows);
  Vg := VectorSpace(K, g);
  Ws := [* *];
  for m in [1..mmax] do
    Append(~Ws, sub< Vg | Basis(Kernel(Submatrix(M, 1, 1, g, e*m))) >);
  end for;
  return Ws;
end function;

function AllKernelsByExpansion(cusps, idx, mmax, h, K, g, t0, pf, kf)
  V := Domain(h);
  B := [h(V.i) : i in [1..g]];
  W := AssociativeArray();
  for i in idx do
    t1 := Cputime();
    Ws := KernelsByExpansion(cusps[i], mmax[i], B, K, g);
    for m in [1..mmax[i]] do
      W[<i,m>] := Ws[m];
      assert Dimension(Ws[m]) ge g - m*Degree(cusps[i]);
      if m gt 1 then assert Ws[m] subset Ws[m-1]; end if;
      if kf ne "" then
        Write(kf, Sprintf("SPACE %o %o %o", i, m, Dimension(Ws[m])));
        for v in Basis(Ws[m]) do Write(kf, &cat[Sprint(x) : x in Eltseq(v)]); end for;
      end if;
    end for;
    Say(pf, Sprintf("EXPKERNEL cusp %o deg %o mmax %o codims %o [%o s, total %o s]",
        i, Degree(cusps[i]), mmax[i], [g - Dimension(Ws[m]) : m in [1..mmax[i]]], Cputime(t1), Cputime(t0)));
  end for;
  return W;
end function;
