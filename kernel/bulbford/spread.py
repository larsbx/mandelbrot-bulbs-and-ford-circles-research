"""The bulb-Ford quantities in exact rational trigonometry: no angle, no π, no complex numbers.

With λ₀ = ζ_q^p the multiplier at the root of B_{p/q}, the first-order prediction is
pred = 2|c'(λ₀)| q⁻² = |1 − λ₀| q⁻², so the exact object is the quadrance

    Qd(1 − λ₀) = |1 − ζ_q^p|² = 4 s,    s = sin²(πp/q) the spread of the half-turn p/q,

and pred² = Qd(1 − λ₀)/q⁴ is Qd(1 − λ₀) times the squared Ford diameter (1/q²)².
Since sin(qθ) = sin θ · U_{q−1}(cos θ), the numbers c_p = cos(πp/q), 0 < p < q, are exactly
the q − 1 simple roots of the integer Chebyshev polynomial U_{q−1}, in decreasing order of
p; and Qd(1 − λ₀) = 4(1 − c_p²). So a Sturm sequence of U_{q−1} over Q isolates c_p by
counting, bisection with dyadic midpoints narrows it, and every bracket returned is a pair
of Fractions. The same polynomials give the spread polynomials of rational trigonometry,
S_n(1 − c²) = (1 − c²) U_{n−1}(c)², with S_n ∘ S_m = S_{nm}.

The Sturm chain is kept in integer polynomials (pseudo-remainders with the sign of the
multiplier fixed, divided by their content), which is practical to q of a few hundred.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache, reduce
from math import gcd

Poly = tuple[int, ...]  # integer coefficients, low -> high
Bracket = tuple[Fraction, Fraction]


# --- integer polynomials -------------------------------------------------------------------------


def trim(p: Poly) -> Poly:
    end = len(p)
    while end > 1 and p[end - 1] == 0:
        end -= 1
    return tuple(p[:end])


def poly_add(a: Poly, b: Poly) -> Poly:
    n = max(len(a), len(b))
    return trim(tuple((a[k] if k < len(a) else 0) + (b[k] if k < len(b) else 0) for k in range(n)))


def poly_scale(a: Poly, c: int) -> Poly:
    return trim(tuple(c * x for x in a))


def poly_mul(a: Poly, b: Poly) -> Poly:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(tuple(out))


def compose(f: Poly, g: Poly) -> Poly:
    """f(g(x)) by Horner over polynomials."""
    return reduce(lambda acc, c: poly_add(poly_mul(acc, g), (c,)), reversed(f), (0,))


@lru_cache(maxsize=None)
def chebyshev_u(n: int) -> Poly:
    """U_n(c): U_0 = 1, U_1 = 2c, U_{k+1} = 2c U_k − U_{k−1}."""
    if n < 0:
        raise ValueError("U_n needs n >= 0")
    if n <= 1:
        return (1,) if n == 0 else (0, 2)
    return poly_add(poly_mul((0, 2), chebyshev_u(n - 1)), poly_scale(chebyshev_u(n - 2), -1))


@lru_cache(maxsize=None)
def spread_polynomial(n: int) -> Poly:
    """S_n(s): S_0 = 0, S_1 = s, S_{k+1} = 2(1 − 2s) S_k − S_{k−1} + 2s; S_{−n} = S_n."""
    n = abs(n)
    if n <= 1:
        return (0,) if n == 0 else (0, 1)
    step = poly_mul((2, -4), spread_polynomial(n - 1))
    return poly_add(poly_add(step, poly_scale(spread_polynomial(n - 2), -1)), (0, 2))


# --- exact signs and Sturm counting --------------------------------------------------------------


def sign_at(p: Poly, x: Fraction) -> int:
    """sign p(a/b) = sign Σ c_k a^k b^(n−k) for b > 0: integer arithmetic only."""
    a, b = x.numerator, x.denominator
    v, b_power = p[-1], 1
    for c in reversed(p[:-1]):                             # homogeneous Horner: O(deg) products
        b_power *= b
        v = v * a + c * b_power
    return (v > 0) - (v < 0)


def _derivative(p: Poly) -> Poly:
    return trim(tuple(k * c for k, c in enumerate(p))[1:]) if len(p) > 1 else (0,)


def _primitive(p: Poly) -> Poly:
    g = reduce(gcd, p, 0)
    return tuple(c // g for c in p) if g > 1 else p


def _negated_pseudo_remainder(a: Poly, b: Poly) -> Poly:
    """A positive multiple of −rem(a, b): from lc(b)^steps a = Q b + r, return −sign(lc^steps) r / content."""
    r, lc, steps = list(a), b[-1], 0
    while len(r) >= len(b) and any(r):
        c, shift = r[-1], len(r) - len(b)
        r = [lc * x for x in r]
        for i, y in enumerate(b):
            r[shift + i] -= c * y
        r = list(trim(tuple(r[:-1]) or (0,)))
        steps += 1
    sign = 1 if lc > 0 or steps % 2 == 0 else -1
    return _primitive(trim(tuple(-sign * x for x in r)))


@lru_cache(maxsize=None)
def sturm_chain(p: Poly) -> tuple[Poly, ...]:
    chain = [p, _derivative(p)]
    while len(chain[-1]) > 1 or chain[-1][0] != 0:
        nxt = _negated_pseudo_remainder(chain[-2], chain[-1])
        if nxt == (0,):
            break
        chain.append(nxt)
    return tuple(chain)


def _variations(chain: tuple[Poly, ...], x: Fraction) -> int:
    signs = [s for s in (sign_at(f, x) for f in chain) if s]
    return sum(1 for u, v in zip(signs, signs[1:]) if u != v)


def roots_in(p: Poly, lo: Fraction, hi: Fraction) -> int:
    """Distinct real roots of a squarefree p in (lo, hi] (Sturm's theorem)."""
    chain = sturm_chain(p)
    return _variations(chain, lo) - _variations(chain, hi)


# --- brackets ------------------------------------------------------------------------------------


def _reduce(p: int, q: int) -> int:
    if q < 1:
        raise ValueError("denominator must be positive")
    r = p % q
    if r == 0:
        raise ValueError("p/q is an integer: the multiplier is 1 and the bulb is the cardioid itself")
    return r


@lru_cache(maxsize=None)
def cosine_isolation(q: int) -> tuple[Bracket, ...]:
    """Intervals (lo, hi], each holding exactly one root of U_{q−1}, in decreasing order of root:
    entry p − 1 holds cos(πp/q). One bisection tree for all roots; each probe is counted once."""
    if q < 2:
        return ()
    chain = sturm_chain(chebyshev_u(q - 1))
    count = lru_cache(maxsize=None)(lambda x: _variations(chain, x))

    def split(lo: Fraction, hi: Fraction) -> list[Bracket]:
        n = count(lo) - count(hi)
        if n <= 1:
            return [(lo, hi)] * n
        mid = (lo + hi) / 2
        return split(mid, hi) + split(lo, mid)

    return tuple(split(Fraction(-1), Fraction(1)))         # U_{q−1}(±1) = (±1)^(q−1) q ≠ 0


@lru_cache(maxsize=None)
def turn_cosine_bracket(p: int, q: int, bits: int = 64) -> Bracket:
    """[lo, hi] ∋ cos(πp/q) for any integer p, of width ≤ 2^−bits.

    cos(πp/q) has period 2q and is even, so r = p mod 2q reduces it to k = min(r, 2q − r); the
    endpoints k = 0, q are ±1 exactly, and otherwise it is the k-th largest root of U_{q−1}.
    """
    if q < 1:
        raise ValueError("denominator must be positive")
    r = p % (2 * q)
    if r % q == 0:
        one = Fraction(1 if r == 0 else -1)
        return one, one
    lo, hi = cosine_isolation(q)[min(r, 2 * q - r) - 1]
    u = chebyshev_u(q - 1)
    # One simple root r in (lo, hi]: u has the sign of u(hi) on (r, hi] and the other sign
    # just below r, so the sign at a midpoint says which half holds r.
    s_hi = sign_at(u, hi)
    if s_hi == 0:
        return hi, hi
    width = Fraction(1, 2**bits)
    while hi - lo > width:
        mid = (lo + hi) / 2
        s = sign_at(u, mid)
        if s == 0:
            return mid, mid
        lo, hi = (lo, mid) if s == s_hi else (mid, hi)
    return lo, hi


def multiplier_quadrance_bracket(p: int, q: int, bits: int = 64) -> Bracket:
    """[lo, hi] ∋ Qd(1 − ζ_q^p) = 4 sin²(πp/q) = 4(1 − c²), of width ≤ 8·2^−bits."""
    lo, hi = turn_cosine_bracket(_reduce(p, q), q, bits)
    sq_lo = Fraction(0) if lo <= 0 <= hi else min(lo * lo, hi * hi)
    sq_hi = max(lo * lo, hi * hi)
    return 4 * (1 - sq_hi), 4 * (1 - sq_lo)


def ford_prediction_squared_bracket(p: int, q: int, bits: int = 64) -> Bracket:
    """[lo, hi] ∋ pred² = (2|c'(λ₀)| q⁻²)² = Qd(1 − λ₀)/q⁴ on the main cardioid."""
    lo, hi = multiplier_quadrance_bracket(p, q, bits)
    return lo / q**4, hi / q**4
