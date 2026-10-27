SetColumns(0);SetMemoryLimit(2*10^9);t0:=Cputime();
AttachSpec("gl2/magma.spec");
for tup in [<33,2>,<34,3>,<35,2>,<38,3>,<39,2>,<40,3>,<42,5>,<44,3>,<45,2>,<46,3>,<48,5>,<49,2>,<50,3>,<51,2>,<52,3>,<54,5>,<55,2>,<56,3>,<57,2>,<60,7>,<63,2>,<64,3>,<65,2>,<69,2>,<72,5>,<75,2>,<77,2>,<85,2>,<91,2>,<95,2>] do
 n,p:=Explode(tup);GG:=GL(2,Integers(n));
 HH:=sub<GG|[GG![1,0,1,1]] cat [GG![a,0,0,1]:a in [1..n]|GCD(a,n) eq 1 ]>;
 zs:=GL2PointCounts(HH,p,8);
 ys:=[zs[d]-GL2RationalCuspCount(HH,p^d):d in [1..8]];
 assert &+GL2jCounts(HH,p) eq ys[1];closed:=[];
 for d in [1..8] do
  z:=ys[d]-&+([0] cat [e*closed[e]:e in [1..d-1]|d mod e eq 0]);
  assert z mod d eq 0;Append(~closed,z div d);
 end for;
 printf "DEGREE8_FAST_COUNT m=1 n=%o p=%o genus=%o rational_point_counts=%o closed_place_counts=%o time=%o\n",n,p,GL2Genus(HH),ys,closed,Cputime(t0);
end for;
print "DEGREE8_FAST_COUNTS_DONE complete=true",Cputime(t0);
quit;
