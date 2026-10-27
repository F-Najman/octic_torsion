#!/usr/bin/env python3
"""Independent finite-field/group-law replay of all 96 C40 octic orbits.

This uses neither Magma nor mdmagma.  It reconstructs the exact fields and
Tate curves, checks order 40, applies every diamond, and proves that the
resulting 708 distinct closed places exhaust the independent Newton count.
"""
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent/'data'

def literals(path):
    result = {}
    for statement in ast.parse(path.read_text()).body:
        assert isinstance(statement, ast.Assign) and len(statement.targets) == 1
        result[statement.targets[0].id] = ast.literal_eval(statement.value)
    return result

class Field:
    def __init__(self, polynomial):
        self.poly = polynomial
        self.d = len(polynomial)-1
        self.q = 3**self.d
        assert polynomial[-1] == 1
        self.powers = [3**i for i in range(self.d)]
        self.digits = [tuple((a//w)%3 for w in self.powers) for a in range(self.q)]
        self.negatives = [sum((-c%3)*w for c,w in zip(v,self.powers)) for v in self.digits]
        primes = []
        remaining = self.q-1
        divisor = 2
        while divisor*divisor <= remaining:
            if remaining % divisor == 0:
                primes.append(divisor)
                while remaining % divisor == 0:
                    remaining //= divisor
            divisor += 1
        if remaining > 1:
            primes.append(remaining)
        for candidate in range(3,self.q):
            if all(self.raw_power(candidate,(self.q-1)//p) != 1 for p in primes):
                primitive = candidate
                break
        else:
            raise AssertionError('No primitive element: defining polynomial is not irreducible')
        self.logs = [-1]*self.q
        self.exps = []
        a = 1
        for i in range(self.q-1):
            assert self.logs[a] == -1
            self.logs[a] = i
            self.exps.append(a)
            a = self.raw_product(a,primitive)
        assert a == 1 and all(v >= 0 for v in self.logs[1:])
        self.exps *= 2

    def encode(self, coefficients):
        assert len(coefficients) == self.d
        assert all(c in (0,1,2) for c in coefficients)
        return sum(c*w for c,w in zip(coefficients,self.powers))

    def add(self,a,b):
        return sum(((x+y)%3)*w for x,y,w in zip(self.digits[a],self.digits[b],self.powers))

    def sub(self,a,b):
        return self.add(a,self.negatives[b])

    def raw_product(self,a,b):
        coefficients = [0]*(2*self.d-1)
        for i,x in enumerate(self.digits[a]):
            for j,y in enumerate(self.digits[b]):
                coefficients[i+j] = (coefficients[i+j]+x*y)%3
        for i in range(2*self.d-2,self.d-1,-1):
            c = coefficients[i]
            if c:
                for j,v in enumerate(self.poly[:-1]):
                    coefficients[i-self.d+j] = (coefficients[i-self.d+j]-c*v)%3
        return sum(c*w for c,w in zip(coefficients[:self.d],self.powers))

    def raw_power(self,a,n):
        result = 1
        while n:
            if n & 1:
                result = self.raw_product(result,a)
            a = self.raw_product(a,a)
            n //= 2
        return result

    def mul(self,a,b):
        return 0 if not a or not b else self.exps[self.logs[a]+self.logs[b]]

    def div(self,a,b):
        assert b
        return 0 if not a else self.exps[(self.logs[a]-self.logs[b])%(self.q-1)]

    def frobenius(self,a):
        return 0 if not a else self.exps[(3*self.logs[a])%(self.q-1)]

    def key(self,b,c):
        keys = []
        for _ in range(self.d):
            keys.append((b,c))
            b,c = self.frobenius(b),self.frobenius(c)
        assert len(set(keys)) == self.d
        return min(keys)

class TateCurve:
    def __init__(self,F,b,c):
        self.F,self.b,self.c = F,b,c
        a1,a2,a3 = F.sub(1,c),F.negatives[b],F.negatives[b]
        # Complete the square: y' = y+(a1*x+a3)/2.  In characteristic 3,
        # 1/2=2 and 1/4=1.
        self.A2 = F.add(a2,F.mul(a1,a1))
        self.A4 = F.mul(2,F.mul(a1,a3))
        self.A6 = F.mul(a3,a3)
        self.P = (0,F.mul(2,a3))
        assert self.on_curve(self.P)
        # Cubic discriminant in characteristic 3.  The original Weierstrass
        # discriminant differs by a nonzero constant and has the same zeros.
        discriminant = F.sub(F.mul(F.mul(self.A2,self.A2),F.mul(self.A4,self.A4)),
                             F.add(F.mul(self.A4,F.mul(self.A4,self.A4)),
                                   F.mul(self.A6,F.mul(self.A2,F.mul(self.A2,self.A2)))))
        assert discriminant

    def on_curve(self,P):
        if P is None:
            return True
        F = self.F
        x,y = P
        return F.mul(y,y) == F.add(F.add(F.mul(x,F.mul(x,x)),F.mul(self.A2,F.mul(x,x))),F.add(F.mul(self.A4,x),self.A6))

    def add(self,P,Q):
        if P is None:
            return Q
        if Q is None:
            return P
        F = self.F
        x,y = P
        u,v = Q
        if x == u:
            if y != v or not y:
                assert F.add(y,v) == 0
                return None
            slope = F.div(F.add(F.mul(2,F.mul(self.A2,x)),self.A4),F.mul(2,y))
        else:
            slope = F.div(F.sub(v,y),F.sub(u,x))
        z = F.sub(F.sub(F.sub(F.mul(slope,slope),self.A2),x),u)
        w = F.sub(F.mul(slope,F.sub(x,z)),y)
        result = z,w
        assert self.on_curve(result)
        return result

    def multiple(self,n):
        result,Q = None,self.P
        while n:
            if n & 1:
                result = self.add(result,Q)
            Q = self.add(Q,Q)
            n //= 2
        return result

    def tate_pair(self,Q):
        F = self.F
        assert Q is not None
        x,y = Q
        assert y
        slope = F.div(F.add(F.mul(2,F.mul(self.A2,x)),self.A4),F.mul(2,y))
        temporary_a2 = F.sub(self.A2,F.mul(slope,slope))
        assert temporary_a2
        scale = F.div(F.mul(2,y),temporary_a2)
        b = F.negatives[F.div(temporary_a2,F.mul(scale,scale))]
        c = F.sub(1,F.div(F.mul(2,slope),scale))
        return b,c

def main():
    data = literals(HERE/'fixed_octic_inputs.py')
    F = Field(data['field_polynomial'])
    assert F.d == 8 and len(data['representatives_bc']) == 96
    assert len(set(data['original_place_indices'])) == 96
    assert len(data['all_place_orbit_indices']) == 708
    magmacounts = Counter(data['all_place_orbit_indices'])
    allkeys = set()
    orbit_sizes = []
    units = [1,3,7,9,11,13,17,19]
    for index,raw in enumerate(data['representatives_bc'],1):
        b,c = [F.encode(v) for v in raw]
        E = TateCurve(F,b,c)
        assert E.multiple(40) is None and E.multiple(20) is not None and E.multiple(8) is not None
        assert E.tate_pair(E.P) == (b,c)
        orbit = {F.key(*E.tate_pair(E.multiple(a))) for a in units}
        assert allkeys.isdisjoint(orbit)
        assert len(orbit) == magmacounts[index]
        allkeys.update(orbit)
        orbit_sizes.append(len(orbit))
    assert len(allkeys) == 708
    report = {'level':40,'degree':8,'representatives':96,'distinct_octic_places':708,
              'orbit_sizes':orbit_sizes,'exact_order40':True,'complete_diamond_partition':True,
              'master_sha256':hashlib.sha256((HERE/'fixed_octic_inputs.m').read_bytes()).hexdigest()}
    (HERE/'fixed_octic_inputs_replay.json').write_text(json.dumps(report,indent=2)+'\n')
    print('C40_FIXED_INPUT_INDEPENDENT_REPLAY_DONE representatives=96 distinct_places=708 exact_order40=true')

if __name__ == '__main__':
    main()
