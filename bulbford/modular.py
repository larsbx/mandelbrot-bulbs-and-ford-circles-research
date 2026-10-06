"""Exact normal-form invariants of f(z) = ζz + z² in ℚ(ζ_q) by multimodular arithmetic (V57).

For a prime p ≡ 1 (mod q) the cyclotomic integers map to 𝔽_p by ζ ↦ ω, one ring homomorphism for each primitive q-th
root ω ∈ 𝔽_p; together they are the φ(q) Galois embeddings.  The homological recursion of normal_form.py has only the
divisions λ₀ − λ₀^k = ζ(1 − ζ^{k−1}), k ≢ 1 (mod q), which are units at p.  So, word-size and exact:
    a = [z^{q+1}] f^q = q·b₁ ∈ ℤ[ζ],      X = q·ι·a² = q²((q² − 1)b₁²/2 + b₂)        (normal_form: qι = (q²−1)/2 + b₂/b₁²).
The φ(q) values of each at the embeddings give its power-basis coordinates (1, ζ, …, ζ^{φ(q)−1}) mod p by a Vandermonde
solve; primes are combined by CRT, a by symmetric lifting (it is integral), X by rational reconstruction (its
denominator is what conjecture C22 is about, so it is not assumed).  Stops when a further batch of primes changes
nothing.  The Galois conjugates σ_j, i.e. the values at λ₀ = e^{2πij/q}, are embed(·, q, j)."""
from __future__ import annotations
from fractions import Fraction
from math import gcd, isqrt
from typing import NamedTuple
import flint
import mpmath as mp
import numpy as np
from .dynamics import njit

P_MAX = 2 ** 31                      # residues < 2^31, so products fit in int64


@njit(cache=True)
def _powmod(x, e, p):
    r, x = 1, x % p
    while e:
        if e & 1:
            r = r * x % p
        x, e = x * x % p, e >> 1
    return r


@njit(cache=True)
def _b12(lam, q, p):
    """(b₁, b₂) mod p of the normal form of λz + z², λ ∈ 𝔽_p a primitive q-th root of unity."""
    H = np.zeros(2 * q + 2, np.int64)
    H[1] = 1
    lp = np.ones(2 * q + 2, np.int64)                       # λ^k
    for k in range(1, 2 * q + 2):
        lp[k] = lp[k - 1] * lam % p
    beta1 = 0
    for k in range(2, 2 * q + 2):
        c = 0                                                # (H²)_k
        for i in range(1, k):
            c = (c + H[i] * H[k - i]) % p
        if k == q + 1:
            beta1 = c
        elif k == 2 * q + 1:
            linv = _powmod(lam, p - 2, p)
            return beta1 * linv % p, c * linv % p
        else:
            rhs = (p - c) % p
            if k > q + 1:
                rhs = (rhs + (k - q) % p * lp[k - q - 1] % p * beta1 % p * H[k - q]) % p
            H[k] = rhs * _powmod((lp[1] - lp[k]) % p, p - 2, p) % p
    return 0, 0


def _primes(q):
    """Primes p ≡ 1 (mod q) below P_MAX, descending."""
    p = (P_MAX - 2) // q * q + 1
    while p > q:
        if flint.fmpz(p).is_prime():
            yield p
        p -= q


def _roots(q, p):
    """The primitive q-th roots of unity ω^j (j coprime to q, ascending) in 𝔽_p."""
    ells = [int(l) for l, _ in flint.fmpz(q).factor()]
    w = next(w for w in (pow(g, (p - 1) // q, p) for g in range(2, p))
             if all(pow(w, q // l, p) != 1 for l in ells))
    return [pow(w, j, p) for j in range(1, q + 1) if gcd(j, q) == 1]


def _coords_mod(q, p):
    """Power-basis coordinates of (a, X) mod p."""
    roots = _roots(q, p)
    half = pow(2, -1, p)
    vals = []
    for w in roots:
        b1, b2 = _b12(w, q, p)
        vals.append((q * b1 % p, q * q % p * (((q * q - 1) % p * half % p * b1 % p * b1 + b2) % p) % p))
    V = flint.nmod_mat([[pow(w, i, p) for i in range(len(roots))] for w in roots], p)
    solve = lambda col: tuple(int(x) for x in V.solve(flint.nmod_mat([[v] for v in col], p)).entries())
    return solve([v[0] for v in vals]), solve([v[1] for v in vals])


def _crt(c, m, r, p):
    """x ≡ c (mod m), x ≡ r (mod p) → x mod mp."""
    return c + m * ((r - c) * pow(m, -1, p) % p)


def _symmetric(c, m):
    return c - m if c > m // 2 else c


def _ratrecon(u, m):
    """n/d ≡ u (mod m) with |n|, d ≤ √(m/2), or None."""
    bound = isqrt(m // 2)
    r0, r1, t0, t1 = m, u % m, 0, 1
    while r1 > bound:
        k = r0 // r1
        r0, r1, t0, t1 = r1, r0 - k * r1, t1, t0 - k * t1
    if t1 == 0 or abs(t1) > bound or gcd(r1, abs(t1)) != 1:
        return None
    return Fraction(r1, t1)


class Exact(NamedTuple):
    q: int
    a: tuple            # integer coordinates of a in 1, ζ, …, ζ^{φ(q)−1}
    X: tuple            # rational coordinates of X = q·ι·a²
    primes: int         # primes used


def exact_invariants(q: int, batch: int = 4) -> Exact:
    """a and X = q·ι·a² of e^{2πi/q}z + z² (λ₀ = ζ) exactly; stable under one further batch of primes."""
    m, ca, cx, last, n = 1, None, None, None, 0
    for p in _primes(q):
        ra, rx = _coords_mod(q, p)
        ca = ra if ca is None else tuple(_crt(c, m, r, p) for c, r in zip(ca, ra))
        cx = rx if cx is None else tuple(_crt(c, m, r, p) for c, r in zip(cx, rx))
        m, n = m * p, n + 1
        if n % batch:
            continue
        X = tuple(_ratrecon(c, m) for c in cx)
        now = (tuple(_symmetric(c, m) for c in ca), X)
        if None not in X and now == last:
            return Exact(q, now[0], X, n)
        last = now
    raise RuntimeError("ran out of primes")


def embed(coords, q: int, j: int = 1):
    """σ_j of the element with these power-basis coordinates: ζ ↦ e^{2πij/q} (current mpmath precision)."""
    z = mp.expjpi(mp.mpf(2 * j) / q)
    return mp.fsum(mp.mpf(c.numerator) / c.denominator * z ** i if isinstance(c, Fraction) else c * z ** i
                   for i, c in enumerate(coords))
