"""Horn map of the ½-parabolic point of the quadratic family (V49): f(z) = −z + z² (c = −3/4, the root shared by the
cardioid's ½-satellite and the period-2 disc).  g = f∘f = z − 2z³ + z⁴ has a triple fixed point at 0 with attracting
petals along ±ℝ and repelling petals along ±iℝ; f exchanges them.

Formal Fatou coordinate of g:  Φ(z) = ¼z⁻² + d₋₁z⁻¹ + β log z + Σ_{k≥1} d_k z^k,  Φ∘g = Φ + 1  (exact rationals; the
equation at order z^{k+2} fixes d_k, the one at z² fixes β).  Horn map from the repelling petal +iℝ: for Im Z above β·π/2 the
orbit of Ψ_rep(Z) arrives in the attracting petal −ℝ, and f (a half step of g) carries it to +ℝ, so
E(Z) = Φ_att(f(g^m(Ψ_rep(Z)))) − m − ½ = Z + a₀ + Σ_{n≥1} a_n e^{2πinZ} on the cylinder of f-orbits."""
from __future__ import annotations
from fractions import Fraction as Fr
from functools import lru_cache
import mpmath as mp

H = {2: Fr(-2), 3: Fr(1)}                      # g(z) = z(1 + h(z)), h = −2z² + z³
HALF = ("½",)                                  # germ tag for renorm.py (length 1: quadratic family, one antipode)


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
    """Fourier coefficients a_0..a_M of the horn map sampled on Im Z = h and the aliasing level |ĉ_{N/2}|.  Above the
    centre β·π/2 ≈ 2.16 (upper map) orbits reach the petal −ℝ and f carries them to +ℝ; below it (lower map) they
    reach +ℝ directly.  a_n carries the noise e^{2πnh}·10^{−dps}."""
    dq, bq = fatou_series_half(K)
    with mp.workdps(dps):
        num = lambda x: mp.mpf(x.numerator) / x.denominator
        d = {k: num(c) for k, c in dq.items()}
        beta = num(bq)
        Phi = lambda z: mp.fsum(c * z ** k for k, c in d.items()) + beta * mp.log(z)
        dPhi = lambda z: mp.fsum(k * c * z ** (k - 1) for k, c in d.items()) + beta / z
        f = lambda z: -z + z * z
        g = lambda z: f(f(z))
        centre = beta * mp.pi / 2                            # Im Φ on the repelling petal (arg z = π/2): the two
                                                             # horn maps live above and below this height
        def E(Z):
            Zp = Z - n
            z = mp.sqrt(1 / (4 * Zp))
            z = z if mp.im(z) > 0 else -z                    # repelling petal +iℝ
            for _ in range(80):
                s = (Phi(z) - Zp) / dPhi(z); z -= s
                if abs(s) < mp.mpf(10) ** (-dps + 5) * abs(z): break
            for _ in range(2 * n):
                z = g(z)
            if abs(z) > 1 or (mp.re(z) > 0) != (h < centre):
                raise ValueError(f"orbit of Z = {Z} did not reach the expected attracting petal")
            return Phi(z if mp.re(z) > 0 else f(z)) - 2 * n - (0 if mp.re(z) > 0 else mp.mpf(1) / 2)
        Zs = [mp.mpf(j) / N + 1j * mp.mpf(h) for j in range(N)]
        D = [E(Z) - Z for Z in Zs]
        coef = lambda k: mp.fsum(D[j] * mp.expjpi(-2 * mp.mpf(k) * j / N) for j in range(N)) / N
        return [coef(k) * mp.exp(2 * mp.pi * k * h) for k in range(M + 1)], abs(coef(N // 2))


def kappa_half(h: float = 4.5, dps: int = 60, **kw):
    """(a₁, a₂, a₂/(2πi a₁²), aliasing) for the ½-parabolic horn map."""
    a, al = horn_coeffs_half(2, h=h, dps=dps, **kw)
    with mp.workdps(dps):
        return a[1], a[2], a[2] / (2j * mp.pi * a[1] ** 2), al


def horn_coeffs_half_normalized(M: int, dps: int):
    """a_0..a_M of the upper horn map translated so that a₀ = 0 and a₁ = 1/(2πi) (Z ↦ Z + c with e^{2πic} =
    1/(2πi a₁): a_n ↦ a_n (2πi a₁)^{−n}), i.e. 𝒫_½(W) = W + W² + …, the form renorm.py expects.  Working precision
    grows with M (the n-th mode loses ≈ 3.8n digits at h = 4.5)."""
    a, _ = horn_coeffs_half(M, h=4.5, dps=dps + 4 * M, K=60, n=1000)
    with mp.workdps(dps + 4 * M):
        s = 2j * mp.pi * a[1]
        return [mp.mpc(0)] + [a[n] / s ** n for n in range(1, M + 1)]
