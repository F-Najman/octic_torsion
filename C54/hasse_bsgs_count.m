/* Exact point count by all annihilator candidates in the Hasse interval.
   No SEA call, extension fields, group structure or point Order calls. */
intrinsic EllipticCountHasseBSGS(E::CrvEll[FldFin] : TrialLimit:=100) -> RngIntElt,RngIntElt
{ Determine #E by exact Hasse-interval BSGS and intersection of random-point annihilators }
 F:=BaseRing(E);q:=#F;RF:=RealField(100);
 bound:=Floor(Sqrt(RF!(4*q)));assert bound^2 le 4*q and 4*q lt (bound+1)^2;
 lo:=q+1-bound;hi:=q+1+bound;width:=hi-lo+1;
 m:=Ceiling(Sqrt(RF!width));assert (m-1)^2 lt width and width le m^2;
 function PointKey(P)
  if P[3] eq 0 then return [F!0,F!1,F!0];end if;
  assert P[3] eq 1;return [P[1],P[2],F!1];
 end function;
 remaining:={Integers()|n:n in [lo..hi]};zero:=E!0;
 for trial in [1..TrialLimit] do
  P:=Random(E);babies:=AssociativeArray();B:=zero;
  for j in [0..m-1] do
   key:=PointKey(B);
   if IsDefined(babies,key) then Append(~babies[key],j);else babies[key]:=[j];end if;
   B+:=P;
  end for;
  A:=-lo*P;step:=-m*P;candidate:={Integers()|};
  for k in [0..(hi-lo) div m] do
   key:=PointKey(A);
   if IsDefined(babies,key) then
    for j in babies[key] do
     n:=lo+k*m+j;
     if n le hi then assert n*P eq zero;Include(~candidate,n);end if;
    end for;
   end if;
   A+:=step;
  end for;
  assert not IsEmpty(candidate);remaining meet:=candidate;
  assert not IsEmpty(remaining);
  if #remaining eq 1 then return Rep(remaining),trial;end if;
 end for;
 error "Hasse BSGS trial cap reached without a unique count";
end intrinsic;
