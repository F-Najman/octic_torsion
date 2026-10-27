/* Native exact abstract generator construction, with profiling call omitted.
   Works after ClassGroupAbelianInvariants(...Proof:=true), and deliberately
   never invokes the expensive inverse-generator divisor reductions. */
procedure MDProvenClassGroupGenerators(F)
 assert assigned F`cld and F`cld`proven;

                      
   // The cyclic generators should
   // 1) consist of small degree places and         
   // 2) have small exponents.
   // 3) The infinite cyclic generator should be the standard 
   //    divisor of degree one.
            
   timeon := Cputime();                      
   vprint ClassGroup : "evaluating abelian group data ...";

   // reverse the fb and rels because then the generators
   // get nicer in the following (is done below as well)

   fb := Reverse( [ fbp`place : fbp in F`cld`fb] );
   Zn := FreeAbelianGroup(#fb);
   N :=  sub< Zn | [ Zn!Reverse(Eltseq(v)) : v in F`cld`rels ] >;
   Cl, g := quo< Zn | N >;
   A1rep := g( Zn![ Valuation(F`cld`A1, P) : P in fb ] );
   h := hom< Cl -> Cl | [ Cl.i : i in [1..NumberOfGenerators(Cl)-1] ] 
                        cat [ A1rep ] >;
   f := hom< Zn -> Cl | g( [ Zn.i : i in [1..NumberOfGenerators(Zn)] ] )@@h >;
   
   genshigh := [ &+ [ DivisorGroup(F) | v[j]*fb[j] : j in [1..#fb] ] 
                       where v is Eltseq(Cl.i@@f) 
               : i in [1..NumberOfGenerators(Cl)-1] ];
   Append(~genshigh, F`cld`A1);                       
                       
   ClToDivHigh := function(a)
      return &+ [ DivisorGroup(F) | v[j] * genshigh[j] : j in [1..#genshigh] ] 
                     where v is Eltseq(a);
   end function;
      
   F`cld`Cl := Cl;   
   F`cld`ZnToCl := f;
   F`cld`ClToDivHigh := ClToDivHigh;
   F`cld`Clabelinv := AbelianInvariants(Cl);
   
   vprint ClassGroup : "done.";

   
end procedure;
