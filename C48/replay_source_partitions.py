"""Independent F5 arithmetic and group-law replay of all level48 inputs.

No Magma code or finite-field implementation is imported. Distinct exact
Tate moduli and Newton's independent count certify exhaustive partitions.
"""
from collections import Counter
import hashlib
import json
from math import gcd
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent / 'data'

class Field:
    def __init__(self, polynomial):
        self.poly = polynomial
        self.p = 5
        self.d = len(polynomial) - 1
        self.q = 5**self.d
        assert polynomial[-1] == 1
        self.powers = [5**i for i in range(self.d)]
        self.digits = [tuple((a // w) % 5 for w in self.powers) for a in range(self.q)]
        self.neg = [sum((-c % 5) * w for c, w in zip(v, self.powers)) for v in self.digits]
        primes, remaining, divisor = [], self.q-1, 2
        while divisor * divisor <= remaining:
            if remaining % divisor == 0:
                primes.append(divisor)
                while remaining % divisor == 0:
                    remaining //= divisor
            divisor += 1
        if remaining > 1:
            primes.append(remaining)
        for primitive in range(5, self.q):
            if all(self.raw_power(primitive, (self.q-1)//r) != 1 for r in primes):
                break
        else:
            raise AssertionError("No primitive element")
        self.logs, self.exps = [-1]*self.q, []
        a = 1
        for i in range(self.q-1):
            assert self.logs[a] == -1
            self.logs[a] = i
            self.exps.append(a)
            a = self.raw_product(a, primitive)
        assert a == 1 and all(x >= 0 for x in self.logs[1:])
        self.exps *= 2

    def encode(self, v):
        assert len(v) == self.d and all(0 <= c < 5 for c in v)
        return sum(c*w for c,w in zip(v, self.powers))

    def add(self, a, b):
        return sum(((x+y) % 5)*w for x,y,w in zip(self.digits[a],self.digits[b],self.powers))

    def sub(self, a, b):
        return self.add(a, self.neg[b])

    def raw_product(self, a, b):
        co = [0]*(2*self.d-1)
        for i,x in enumerate(self.digits[a]):
            for j,y in enumerate(self.digits[b]):
                co[i+j] = (co[i+j] + x*y) % 5
        for i in range(2*self.d-2, self.d-1, -1):
            if co[i]:
                c = co[i]
                for j,v in enumerate(self.poly[:-1]):
                    co[i-self.d+j] = (co[i-self.d+j] - c*v) % 5
        return sum(c*w for c,w in zip(co[:self.d], self.powers))

    def raw_power(self, a, n):
        ans = 1
        while n:
            if n & 1:
                ans = self.raw_product(ans, a)
            a = self.raw_product(a, a)
            n //= 2
        return ans

    def mul(self, a, b):
        return 0 if not a or not b else self.exps[self.logs[a] + self.logs[b]]

    def div(self, a, b):
        assert b
        return 0 if not a else self.exps[(self.logs[a] - self.logs[b]) % (self.q-1)]

    def frob(self, a):
        return 0 if not a else self.exps[(5*self.logs[a]) % (self.q-1)]

    def key(self, b, c):
        out = []
        for _ in range(self.d):
            out.append((b,c))
            b,c = self.frob(b), self.frob(c)
        assert len(set(out)) == self.d
        return min(out)

class Curve:
    def __init__(self, F, b, c):
        self.F, self.b, self.c = F,b,c
        a1, a2, a3 = F.sub(1,c), F.neg[b], F.neg[b]
        self.A2 = F.add(a2, F.div(F.mul(a1,a1),4))
        self.A4 = F.div(F.mul(a1,a3),2)
        self.A6 = F.div(F.mul(a3,a3),4)
        self.P = (0,F.div(a3,2))
        disc_factor = F.add(F.sub(F.mul(1,F.mul(b,b)), F.mul(b,F.sub(F.mul(3,F.mul(c,c)),1))),
                            F.mul(c,F.mul(F.sub(c,1),F.mul(F.sub(c,1),F.sub(c,1)))))
        # This is 16b^2-b(8c^2+20c-1)+c(c-1)^3, reduced modulo5.
        assert b and disc_factor and self.on_curve(self.P)

    def on_curve(self, P):
        if P is None:
            return True
        F=self.F
        x,y=P
        return F.mul(y,y) == F.add(F.add(F.mul(x,F.mul(x,x)),F.mul(self.A2,F.mul(x,x))),F.add(F.mul(self.A4,x),self.A6))

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        F=self.F
        x,y=P
        u,v=Q
        if x == u:
            if y != v or not y:
                assert F.add(y,v) == 0
                return None
            slope = F.div(F.add(F.add(F.mul(3,F.mul(x,x)),F.mul(2,F.mul(self.A2,x))),self.A4),F.mul(2,y))
        else:
            slope=F.div(F.sub(v,y),F.sub(u,x))
        z=F.sub(F.sub(F.sub(F.mul(slope,slope),self.A2),x),u)
        w=F.sub(F.mul(slope,F.sub(x,z)),y)
        assert self.on_curve((z,w))
        return z,w

    def multiple(self, n):
        ans,Q=None,self.P
        while n:
            if n & 1:
                ans=self.add(ans,Q)
            Q=self.add(Q,Q)
            n//=2
        return ans

    def tate_pair(self, Q):
        F=self.F
        assert Q is not None
        x,y=Q
        assert y
        slope=F.div(F.add(F.add(F.mul(3,F.mul(x,x)),F.mul(2,F.mul(self.A2,x))),self.A4),F.mul(2,y))
        tmp=F.sub(F.add(self.A2,F.mul(3,x)),F.mul(slope,slope))
        assert tmp
        scale=F.div(F.mul(2,y),tmp)
        b=F.neg[F.div(tmp,F.mul(scale,scale))]
        c=F.sub(1,F.div(F.mul(2,slope),scale))
        return b,c

def main():
    src=HERE/'source_partition_inputs_complete.py'
    a=runpy.run_path(str(src))
    assert a['complete'] is True and a['prime'] == 5 and a['level'] == 48
    report={}
    for d in range(3,9):
        F=Field(a['field_polynomials'][d])
        assert F.d == d
        seen=set()
        sizes=[]
        for raw in a['representative_pairs'][d]:
            b,c=map(F.encode,raw)
            E=Curve(F,b,c)
            assert E.multiple(48) is None and E.multiple(24) is not None and E.multiple(16) is not None
            assert E.tate_pair(E.P) == (b,c)
            orbit={F.key(*E.tate_pair(E.multiple(u))) for u in range(1,25) if gcd(u,48)==1}
            assert seen.isdisjoint(orbit)
            seen.update(orbit)
            sizes.append(len(orbit))
        assert sizes == a['partition_orbit_sizes'][d]
        assert len(seen) == a['noncuspidal_counts'][d] == [16,64,656,2712,10784,48992][d-3]
        report[d]={'representatives':len(sizes),'places':len(seen),'orbit_size_distribution':dict(Counter(sizes))}
        print('C48_SOURCE_PARTITION_INDEPENDENT_PASS',d,report[d],flush=True)
    result={'complete':True,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
            'all_points_exact_order48':True,'all_diamond_orbits_disjoint':True,
            'all_newton_counts_exhausted':True,'degrees':report}
    (HERE/'source_partition_replay.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    main()
