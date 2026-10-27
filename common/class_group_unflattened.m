ProjectiveNumberingMap := function(B, k)
//
//  A set {1..n} and map: f: {1..n} -> Universe(B) with image
//  " ProjectiveHull( Span_k(B) ) "
//
   k := SetToSequence(Set(k));
   len := #k;
   NumToElt := function(n)
      if n eq -1 then
         return B[1];
      end if;
      p := len;
      i := 2;
      while n - p ge 0 do
         n -:= p;
         p *:= len;
         i +:= 1; 
      end while;
      b := Zero(Universe(B));
      j := 1;
      repeat
         n, r := Quotrem(n, len);
         b := b + k[r+1]*B[j];
         j +:= 1;
      until n eq 0;
      b := b + B[i];
      return b;
   end function;
   return map< { 1..(len^#B-1) div (len-1) } -> Universe(B) 
                  | n :-> NumToElt(n-2) >;
end function;
   
   
/* Exact forward class map from Magma2.30-1 package/Ring/FldFun/classgroup.m.
   Bypasses ONLY the inverse-generator flattening step; all forward logic below
   is copied byte-for-byte. Requires an already PROVEN class group and its
   evaluated abstract group data. No relation search or unproven group used. */
function MDClassGroupUnflattened(F0)
 F:=RationalExtensionRepresentation(F0);
 assert assigned F`cld and assigned F`cld`proven and F`cld`proven;
 assert assigned F`cld`Cl and assigned F`cld`ZnToCl and assigned F`cld`ClToDivHigh;
   f := F`cld`ZnToCl;
   Zn := Domain(f);             
   Cl := Codomain(f);
   ClToDivHigh := F`cld`ClToDivHigh;
   fb := Reverse( [ fbp`place : fbp in F`cld`fb] );
   A1 := F`cld`A1;             
   
   ClgDivisorFactorFb := function(F, D)
      // only apply this to reduced (small) divisors!
      I, _ := Ideals(D); // rem: does a superflous ideal inversion
      N := Factorisation(Minimum(I));
      ecfdim := DimensionOfExactConstantField(F);
      if #N gt 0 and
         Maximum([ Degree(x[1]) : x in N ]) gt F`cld`fbmaxdeg * ecfdim then
         return false;
      end if;
      places, expos := Support(D);
      v := Universe(F`cld`rels)!0;
      for i in [1..#places] do
         j := Index(fb, places[i]);
         if j eq 0 then
            return false;
         end if;
         v[F`cld`fb[j]`i] := expos[i];
      end for;
      return v;
   end function;
      
   DivToCl := function(D)
      
      F := FunctionField(D);
       
      if assigned F`cld`A and not assigned F`cld`Avec then
         F`cld`Avec := ClgDivisorFactorFb(F, F`cld`A);
         assert F`cld`Avec cmpne false;
      end if;
      
      places, expos := Support(D);
      v := Universe(F`cld`rels)!0;
      for i in [1..#places] do
         j := Index(fb, places[i]);
         if j eq 0 then
            v := false;
            break;
         end if;
         v[F`cld`fb[j]`i] := expos[i];
      end for;
      
      if v cmpne false then
         vprint ClassGroup :
           "is smooth";
      else
         
         Dred1, r1, xinfty := Reduction(D);
         v := ClgDivisorFactorFb(F, Dred1);
      
         if not assigned F`cld`xinftyvec or F`cld`xinfty ne xinfty then
            F`cld`xinfty := xinfty;
            // is smooth or something is wrong
            F`cld`xinftyvec := ClgDivisorFactorFb(F, xinfty);
         end if;
      
         if v cmpne false then
            vprint ClassGroup : 
              "first reduction ok";
            v := v + r1*F`cld`xinftyvec;
         else
            vprint ClassGroup : 
              "try several random reductions ...";
            tries := 0;
            n := 1;
            while v cmpeq false do
               tries +:= 1;
               n := Minimum( [n, #F`cld`fb, 10] );
               Drand := &+ [ DivisorGroup(F) | 
                             Random([1..n]) * Random(F`cld`fb)`place 
                           : i in [1..n] ];
               if assigned F`cld`A then
                  Dred2, r2 := Reduction(Dred1 + Drand, F`cld`A);
               else
                  Dred2, r2, xinfty2 := Reduction(Dred1 + Drand);
                  // should be the case otherwise we have to factor xinfty2
                  assert xinfty2 eq F`cld`xinfty;
               end if;
               N := ProjectiveNumberingMap(Basis(Dred2), ConstantField(F));
               s := Domain(N);
               vprint ClassGroup :
                 "enum number", Maximum(s), "at weight", n;
               i := 1;
               while v cmpeq false and i in s do
                  v := ClgDivisorFactorFb(F, Dred2 + Divisor(N(i)));
                  i := i + 1;
               end while;
               if tries gt Maximum([ #F`cld`fb / 10, 10 ]) then
                  tries := 0;
                  n +:= 1;
               end if;
            end while;
            Drandvec := ClgDivisorFactorFb(F, Drand);
            if assigned F`cld`A then
               v := v - Drandvec + r2*F`cld`Avec + r1*F`cld`xinftyvec;
            else
               v := v - Drandvec + (r1 + r2)*F`cld`xinftyvec;
            end if;
            
         end if;
      end if;
      
      return f(Zn!Eltseq(v));
      
   end function;
      

 phi:=hom<DivisorGroup(F)->Cl|x:->DivToCl(x),y:->ClToDivHigh(y)>;
 if F0 cmpeq F then return Cl,Inverse(phi),phi;end if;
 psi:=hom<DivisorGroup(F0)->Cl|x:->phi(DivisorGroup(F)!x),y:->DivisorGroup(F0)!(y@@phi)>;
 return Cl,Inverse(psi),psi;
end function;
