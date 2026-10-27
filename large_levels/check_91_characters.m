SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);
AttachSpec("mdmagma/Magma/magma.spec");
N:=91;a:=3;
for lab in ["g","u"] do
 idx:=[j:j in [1..N-1]|GCD(j,N) eq 1 and ConreyCharacterOrbitLabel(N,j) eq Sprintf("91.%o",lab)];
 assert #idx gt 0 and &and[ConreyCharacterValue(N,j,a) eq 1:j in idx];
 printf "OCTIC91_CHARACTER_KERNEL orbit=%o indices=%o value_at3=1 verified=true\n",lab,idx;
end for;
print "OCTIC91_CHARACTER_CERTIFICATE_DONE";
quit;
