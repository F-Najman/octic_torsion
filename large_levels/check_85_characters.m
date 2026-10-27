SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);
AttachSpec("mdmagma/Magma/magma.spec");
N:=85;a:=2;
idx:=[j:j in [1..N-1]|GCD(j,N) eq 1 and ConreyCharacterOrbitLabel(N,j) eq "85.j"];
assert #idx gt 0 and &and[ConreyCharacterValue(N,j,a) eq 1:j in idx];
printf "OCTIC85_CHARACTER_KERNEL indices=%o value_at2=1 verified=true\n",idx;
print "OCTIC85_CHARACTER_CERTIFICATE_DONE";
quit;
