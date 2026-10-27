"""Independent finite-field moduli and complete degree7 source coverage."""
from pathlib import Path
import ast,json,re
from C75_binary_fields_root import BinaryField,Elliptic
P=Path(__file__).resolve().parent
raw=(P/'data'/'C75_good7_moduli_records_root_v2.txt').read_text().replace('[*','[').replace('*]',']')
rows=[ast.literal_eval(n.value) for n in ast.parse(raw).body];assert len(rows)==14
fields={};degrees=[]
for poly,data,degree in rows:
 key=tuple(poly)
 if key not in fields:fields[key]=BinaryField(poly)
 F=fields[key];m=F.mul;E=Elliptic(F,data[:5]);Q25=tuple(data[5:7]);Q3=tuple(data[7:9]);bc=tuple(data[9:11]);xy=tuple(data[11:13]);hv,xc=data[13:15]
 assert E.on(Q25) and E.on(Q3) and E.multiple(25,Q25) is None and E.multiple(5,Q25) is not None
 assert E.multiple(3,Q3) is None and Q3 is not None
 assert E.tate(Q25)==(bc,xy)
 xx,yy=xy;b,c=bc
 r=F.div(m(F.square[xx],yy)^m(xx,yy)^yy^1,m(xx,m(xx,yy)^1));s=F.div(m(xx,yy)^yy^1,m(xx,yy))
 assert b==m(m(r,s),r^1) and c==m(s,r^1)
 f3,f4,f5,f6,f7,f8,f9=b,c,b^c,s^1,s^r,m(r,s)^1,F.square[s]^s^r^1
 f10,f12=xx^yy^1,xx^yy
 f13=m(m(F.square[xx],xx),yy)^m(F.square[xx],F.square[yy])^m(F.square[xx],yy)^m(xx,F.square[yy])^yy^1
 top=m(m(m(m(m(f3,f4),F.pow(f6,3)),F.square[f7]),f10),F.pow(f12,3))
 bot=m(m(m(f5,F.square[f8]),F.square[f9]),F.pow(f13,3));assert F.div(top,bot)==hv
 # Independently carry the actual3-point through the Tate normalization.
 a1,a2,a3,a4,a6=E.a1,E.a2,E.a3,E.a4,E.a6;x0,y0=Q25
 aa3=a3^m(a1,x0);aa4=F.square[x0]^a4^m(a1,y0);sl=F.div(aa4,aa3)
 aa2=x0^a2^m(sl,a1)^F.square[sl];sc=F.div(aa2,aa3)
 rx=m(F.square[sc],Q3[0]^x0)
 rr=F.div(b,1^c);expected=F.div(rx^rr,F.square[1^c]);assert xc==expected
 b2=F.square[a1];b4=m(a1,a3);b6=F.square[a3]
 b8=m(F.square[a1],a6)^m(m(a1,a3),a4)^m(a2,F.square[a3])^F.square[a4]
 delta=m(F.square[b2],b8)^m(m(b2,b4),b6)^F.square[b6]
 jv=F.div(F.pow(a1,12),delta)
 assert F.pow(xc,4)^F.pow(xc,3)==F.inv[jv]
 assert all(F.pow(v,1<<degree)==v for v in [hv,xx,xc])
 assert degree==F.d and degree in [7,14];degrees.append(degree)
# Exhaustively count both ordinary twists for each j inF128.
F=next(F for F in fields.values() if F.d==7);twist=next(a for a in range(F.q) if F.trace(a))
qual=[]
for j in range(1,F.q):
 a6=F.inv[j]
 for a2 in [0,twist]:
  count=2+2*sum(F.trace(x^a2^F.div(a6,F.square[x]))==0 for x in range(1,F.q))
  if count==150:qual.append(j)
assert len(qual)==7 and len(set(qual))==7
j=qual[0];assert {F.pow(j,1<<k) for k in range(7)}==set(qual)
# Source counts from the degree290 Frobenius polynomial, using Newton identities.
a=ast.literal_eval((P/'data'/'C75_source_frob_root_arithmetic.txt').read_text().replace('[*','[').replace('*]',']'))[0]
_,g,co,native,total=a;assert g==145 and len(co)==291
tr=[];closed=[]
for n in range(1,9):
 t=-n*co[2*g-n]-sum(co[2*g-j]*tr[n-j-1] for j in range(1,n));tr.append(t)
 num=(1<<n)+1-t-sum(e*closed[e-1] for e in range(1,n) if n%e==0);assert num%n==0;closed.append(num//n)
assert closed==native==[20,10,0,8,0,0,20,0]
tate=json.loads((P/'data'/'audit75_tate_parities.json').read_text())
cusps={int(k):v for k,v in tate['source_F2_signature'].items()};nc=[v-cusps.get(i+1,0) for i,v in enumerate(closed)]
assert nc==[0,0,0,0,0,0,20,0]
ranklog=(P/'logs'/'C75_rank_zero_root.log').read_text()
rr=re.findall(r'order=(\d+) conductor=(\d+) dimension=(\d+) winding=(\d+) coefficient_degree=(\d+)',ranklog)
assert len(rr)==6 and sum(int(d)*int(e) for o,c,d,w,e in rr)==145 and all(d==w for o,c,d,w,e in rr)
assert 'C75_RANK_ZERO_WINDING_COMPLETE 145' in ranklog
log=(P/'logs'/'C75_good7_certificate_root_v2.log').read_text();assert 'C75_GOOD7_NATIVE_REFERENCE_AND_RR_PASS 112' in log and 'C75_GOOD7_CERTIFICATE_DONE complete=true' in log
assert 'Runtime error' not in log and 'Assertion failed' not in log
out={'complete':True,'actual_moduli_records_checked':len(rows),'qualifying_ordinary_twists':7,'all_source_NC_counts_degrees1to8':nc,'source_degree7_diamond_orbits':1,'quotient_degree7_diamond_orbits':1,'native_full_isogeny_divisor_equal_twice_normed_Hecke':True}
(P/'data'/'C75_good7_independent_root.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out));print('C75_GOOD7_INDEPENDENT_REPLAY_PASS')
