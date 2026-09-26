"""Horn map of the ½-parabolic point of the quadratic family (V49): f(z) = −z + z² (c = −3/4, the root shared by the
cardioid's ½-satellite and the period-2 disc).  g = f∘f = z − 2z³ + z⁴ has a triple fixed point at 0 with attracting
petals along ±ℝ and repelling petals along ±iℝ; f exchanges them.

Exact formal Fatou coordinate of g (Lemma 7.2 of the paper):  Φ(z) = ¼z⁻² + d₋₁z⁻¹ + β log z + Σ_{k≥1} d_k z^k,
Φ∘g = Φ + 1, in rationals (the equation at order z^{k+2} fixes d_k, the one at z² fixes β = 11/8).  The horn map is the
case p/q = 1/2 of horn_rational (upper map: orbits from +iℝ above the escape band reach −ℝ, and f carries them to +ℝ)."""
from __future__ import annotations
from fractions import Fraction as Fr
from functools import lru_cache
import mpmath as mp
from .horn_rational import horn_coeffs_pq, root

HALF = root(1, 2)                              # germ tag for renorm.py

H = {2: Fr(-2), 3: Fr(1)}                      # g(z) = z(1 + h(z)), h = −2z² + z³


def _mul(a: list, b: list, L: int) -> list:
    return [sum((a[i] * b[k - i] for i in range(k + 1)), Fr(0)) for k in range(L + 1)]


def _pow1p(k: int, L: int) -> list:
    """(1 + h)^k as a power series to order L, any integer k (generalized binomial)."""
    h = [H.get(j, Fr(0)) for j in range(L + 1)]
    out, term, binom = [Fr(0)] * (L + 1), [Fr(1)] + [Fr(0)] * L, Fr(1)
    for j in range(L + 1):
        out = [o + binom * t for o, t in zip(out, term)]
        term, binom = _mul(term, h, L), binom * (k - j) / (j + 1)
    return out


def _log1p(L: int) -> list:
    h = [H.get(j, Fr(0)) for j in range(L + 1)]
    out, term = [Fr(0)] * (L + 1), [Fr(1)] + [Fr(0)] * L
    for j in range(1, L + 1):
        term = _mul(term, h, L)
        out = [o + Fr((-1) ** (j + 1), j) * t for o, t in zip(out, term)]
    return out


@lru_cache(maxsize=None)
def fatou_series_half(K: int = 30) -> tuple[dict, Fr]:
    """({k: d_k}, β) for k = −2, −1, 1, …, K."""
    L = K + 4
    P = {k: _pow1p(k, L) for k in range(-2, K + 1) if k}
    lg = _log1p(L)
    d, beta = {}, Fr(0)
    for m in range(K + 3):                                  # order z^m of Φ∘g − Φ − 1
        known = sum((c * P[k][m - k] for k, c in d.items() if 0 <= m - k <= L), Fr(0)) + beta * lg[m]
        rhs = (1 if m == 0 else 0) - known
        if m == 2:
            beta = rhs / lg[2]
        elif m - 2 <= K:
            k = m - 2
            d[k] = rhs / P[k][2]                            # [z²]((1+h)^k − 1) = −2k
    return d, beta


def horn_coeffs_half(M: int, h: float = 4.5, N: int = 64, n: int = 600, dps: int = 40, K: int = 30):
    """Upper horn map of the ½-root sampled at the absolute height h (above the escape band ≈ [1.5, 3]); the general
    construction is horn_rational.horn_coeffs_pq(1, 2, …)."""
    return horn_coeffs_pq(1, 2, M, upper=True, h=h, N=N, n=n, K=K, dps=dps)


def kappa_half(h: float = 4.5, dps: int = 60, **kw):
    """(a₁, a₂, a₂/(2πi a₁²), aliasing) for the ½-parabolic horn map."""
    a, al = horn_coeffs_half(2, h=h, dps=dps, **kw)
    with mp.workdps(dps):
        return a[1], a[2], a[2] / (2j * mp.pi * a[1] ** 2), al
