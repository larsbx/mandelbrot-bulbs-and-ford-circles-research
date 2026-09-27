"""Exact ι_{p/q} ∈ ℚ(ζ_q): series arithmetic over ℚ[x]/Φ_q(x). Prints basis coordinates + float."""
from __future__ import annotations
from fractions import Fraction as Fr
from functools import reduce
from itertools import product
import cmath, sys
from math import pi
from sympy import cyclotomic_poly, Poly, symbols

x = symbols('x')

def cyclo(q):
    return tuple(Fr(int(c)) for c in reversed(Poly(cyclotomic_poly(q, x), x).all_coeffs()))  # low→high

def field(q):
    phi = cyclo(q); n = len(phi) - 1                       # monic degree n = φ(q)
    def reduce_(v):                                          # reduce polynomial (low→high) mod Φ_q
        v = list(v) + [Fr(0)] * max(0, n - len(v))
        for k in range(len(v) - 1, n - 1, -1):
            if v[k]:
                c = v[k]
                for i in range(n + 1):
                    v[k - n + i] -= c * phi[i]
        return tuple(v[:n])
    def mul(a, b):
        out = [Fr(0)] * (2 * n - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj: out[i + j] += ai * bj
        return reduce_(out)
    add = lambda a, b: tuple(u + v for u, v in zip(a, b))
    zero = tuple(Fr(0) for _ in range(n)); one = reduce_([Fr(1)])
    zeta = reduce_([Fr(0), Fr(1)])
    return n, zero, one, zeta, add, mul, reduce_

def exact_index(p, q):
    n, zero, one, zeta, add, mul, red = field(q)
    scal = lambda c, a: tuple(c * t for t in a)
    lam = reduce(lambda a, _: mul(a, zeta), range(p), one)  # ζ^p
    N = 2 * q + 2
    smul = lambda A, B: tuple(reduce(add, (mul(A[i], B[k - i]) for i in range(k + 1)), zero) for k in range(N))
    s = (zero, one) + (zero,) * (N - 2)
    for _ in range(q):
        s = tuple(add(mul(lam, a), b) for a, b in zip(s, smul(s, s)))
    zminus = tuple(add(scal(Fr(1 if k == 1 else 0), one) if k == 1 else zero, scal(Fr(-1), a)) for k, a in enumerate(s))
    assert all(v == zero for v in zminus[:q + 1]), "z^2..z^q coefficients must vanish"
    P = zminus[q + 1:]                                       # length q+1
    # inverse series of P to order q: I_0 = 1/P_0 ... field inverse via linear algebra is heavy; use
    # Newton-free approach: solve P * I = 1 coefficientwise needs 1/P_0 in the field -> use sympy.
    from sympy import Rational, Symbol, rem, invert
    P0 = sum(Rational(c.numerator, c.denominator) * x ** i for i, c in enumerate(P[0]))
    inv0 = Poly(invert(P0, cyclotomic_poly(q, x)), x).all_coeffs()[::-1]
    inv0 = red([Fr(int(c.p), int(c.q)) for c in inv0])
    I = [inv0]
    for k in range(1, q + 1):
        acc = reduce(add, (mul(P[i], I[k - i]) for i in range(1, k + 1)), zero)
        I.append(mul(scal(Fr(-1), inv0), acc))
    return I[q]

def to_complex(v, q):
    return sum(complex(c) * cmath.exp(2j * pi * i / q) for i, c in enumerate(v))

if __name__ == "__main__":
    qs = [int(a) for a in sys.argv[1:]] or list(range(1, 13))
    for q in qs:
        for p in [1] + ([q - 1] if q > 2 else []):
            if q == 1: continue
            v = exact_index(p, q)
            print(f"q={q} p={p} ι = {to_complex(v, q):.12g}   coords(1,ζ,…) = {[str(c) for c in v]}", flush=True)
