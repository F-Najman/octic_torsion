import sys;from pathlib import Path;sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'common'))
from finite_field5_checked import FF,EC,pack
import ast,json
from pathlib import Path
p=Path(__file__).resolve().parent/'data'
a=ast.literal_eval((p/'M2N36_D8_direct_full8_arithmetic.txt').read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
f=a[0];d,ky=a[7][1605];assert d==8
F=FF(f);b,c=pack(ky[:8]),pack(ky[8:]);E=EC(F,b,c);Q=(b,b);P=(0,0);qs=[None]
for _ in range(36):qs.append(E.add(qs[-1],Q))
st=[];ims=[]
for bb in [0,1]:
 for cc in [0,1]:
  for dd in [1,5,7,11,13,17]:
   PP=E.add(P,qs[18] if bb else None);QQ=E.add(P if cc else None,qs[dd]);key=F.key(E.bc(PP,QQ));ims.append(key)
   if key==(b,c):st.append([bb,cc,dd])
assert len(set(ims))==6
assert st==[[0,0,1],[0,0,17],[1,1,1],[1,1,17]]
out={'source_record':1606,'field_polynomial':f,'key':ky,'orbit_size':len(set(ims)),'stabilizer_parameters':st}
(p/'M2N36_survivor_Aut24_stabilizer.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));print('M2N36_SURVIVOR_AUT24_STABILIZER_DONE')
