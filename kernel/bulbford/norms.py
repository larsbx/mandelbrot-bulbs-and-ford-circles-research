"""The parabolic coefficient a_q = [z^{q+1}] f^q(z), f(z) = λz + z², and its norm from ℚ(ζ_q) to ℚ.

At λ = ζ_q, f^q(z) = z + a_q z^{q+1} + …, and ι = b/a_q² in normal form, so the primes of
N(a_q) bound the denominators of ι. At three multipliers the map is exactly solvable:

- λ = 2:  2z + z² = (1 + z)² − 1 is conjugate to w ↦ w², so a_q(2) = C(2^q, q + 1);
- λ = 4:  4z + z² is conjugate to the Chebyshev map w ↦ w² − 2 at its fixed point w = 2, and
          2T_N(1 + z/2) gives a_q(4) = 2N (N + q)! / ((N − q − 1)! (2q + 2)!), N = 2^q;
- λ = −2: the same Chebyshev map at its interior fixed point w = −1 (no closed form used).

Hence (bridges spike, B3): if ℓ > 2q + 2 is prime and ord_ℓ(2) = q, a prime 𝔩 | ℓ of ℤ[ζ_q] with
ζ_q ≡ 2 (mod 𝔩) contains a_q(ζ_q) ≡ C(2^q, q + 1) ≡ 0, because ℓ | 2^q − 1 divides the numerator and
not (q + 1)!; likewise for ord_ℓ(4) = q through the factor N² − 1 = 4^q − 1.
"""
from __future__ import annotations

from math import comb, factorial, gcd

import mpmath as mp


def a_coefficient(lam, q: int):
    """[z^{q+1}] f^q(z) by truncated series; `lam` may be an int, a residue, or an mp number."""
    zero = lam * 0
    s = [zero, zero + 1] + [zero] * q
    for _ in range(q):
        s = [lam * s[n] + sum(s[i] * s[n - i] for i in range(n + 1)) for n in range(q + 2)]
    return s[q + 1]


def a_at_2(q: int) -> int:
    """a_q(2) = C(2^q, q + 1) (power-map conjugacy)."""
    return comb(2**q, q + 1)


def a_at_4(q: int) -> int:
    """a_q(4) = 2N (N + q)! / ((N − q − 1)! (2q + 2)!), N = 2^q (Chebyshev conjugacy)."""
    n, k = 2**q, q + 1
    return 2 * n * factorial(n + k - 1) // (factorial(n - k) * factorial(2 * k))


def norm_a(q: int) -> int:
    """|N(a_q)| = |Π_{p ∈ (ℤ/q)^×} a_q(ζ_q^p)|, an integer, from a product at 40 + 8q digits."""
    with mp.workdps(40 + 8 * q):
        prod = mp.fprod(a_coefficient(mp.expj(2 * mp.pi * p / q), q) for p in range(1, max(q, 2)) if gcd(p, q) == 1)
        n = int(mp.nint(mp.re(prod)))
        if abs(prod - n) > mp.mpf(10) ** -10:
            raise ArithmeticError(f"N(a_{q}) is not resolved at this precision")
    return abs(n)
