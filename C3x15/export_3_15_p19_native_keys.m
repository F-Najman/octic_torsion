/* Native degree4 key export, one fresh bounded process per shard.
   Every key is an actual normalized-curve place, and the complete
   source count is checked independently from Frobenius. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(6*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_full_aut_hash_m3_15.m");
Attach("C3x15/fast_full_aut_hash_m3_15_native.m");
// The shard is the parameter Shard (magma -b Shard:=k); the original runs set the environment variable M3N15_SHARD.
if assigned Shard then shard:=StringToInteger(Shard); else shard:=StringToInteger(GetEnv("M3N15_SHARD")); end if;assert shard in [1..6];
X:=MDX11(3,15,GF(19));C:=Curve(X);AF,mm:=AlgorithmicFunctionField(FunctionField(C));
coords:=[mm(f):f in X`_coordinates];assert #[c:c in Cusps(X)|Degree(c) eq 4] eq 0;
xs:=Places(AF,4);assert #xs eq 32076;
J:=JH(45,[16]);fp:=FrobeniusPolynomial(J,19);assert Dimension(J) eq 9;
tr:=[];cl:=[];
for n in [1..4] do
 sn:=-n*Coefficient(fp,18-n);
 for j in [1..n-1] do sn-:=Coefficient(fp,18-j)*tr[n-j];end for;
 Append(~tr,Integers()!sn);num:=19^n+1-Integers()!sn;
 for e in Divisors(n) do if e lt n then num-:=e*cl[e];end if;end for;
 assert num mod n eq 0;Append(~cl,num div n);
end for;
assert cl[4] eq #xs;
lo:=1+(shard-1)*5500;hi:=Minimum(shard*5500,#xs);records:=[* *];
print "M3N15_P19_NATIVE_SHARD_START",shard,lo,hi,#xs,cl,Cputime(t0);
for i in [lo..hi] do
 values:=[Evaluate(f,xs[i]):f in coords];key:=M3FifteenNativeKey(values,19,4);
 Append(~records,[*i,key*]);
 if (i-lo+1) mod 500 eq 0 then print "M3N15_P19_NATIVE_SHARD_PROGRESS",shard,i,Cputime(t0);end if;
end for;
Write("C3x15/data/M3N15_p19_native_keys_shard" cat IntegerToString(shard) cat ".txt",
 Sprint([*shard,lo,hi,#xs,Coefficients(fp),records*]):Overwrite:=true);
print "M3N15_P19_NATIVE_SHARD_DONE",shard,#records,"complete=true",Cputime(t0);quit;
