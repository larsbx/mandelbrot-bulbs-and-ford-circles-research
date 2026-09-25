"""Horn maps of parabolic polynomial germs g(v) = v + v² + g₃v³ + … and κ₀ = a₂/(2πi a₁²) (V38, V41).

In w = −1/v the global map is F(w) = −1/g(−1/w) = w + 1 + β/w + …, β = 1 − g₃ (résidu itératif).
Formal Fatou coordinate Φ(w) = w − β log w + Σ_{k≥1} c_k w^{−k}, exact rationals from Φ∘F = Φ + 1.
Φ_att = Φ(F^n w) − n (principal log); Ψ_rep(Z) = F^n(Φ⁻¹(Z − n)) (log with arg ∈ (0, 2π));
upper horn map E(Z) = Φ_att(Ψ_rep(Z)) = Z + a₀ + Σ_{n≥1} a_n e^{2πinZ}.  In W = e^{2πiZ} the phase-0 germ is
𝒫₀(W) = W exp(2πi Σ_{n≥1} a_n W^n), and ι(𝒫₀) − ½ = a₂/(2πi a₁²).

Germs are tuples of Fractions (g₂, g₃, …) = coefficients of v², v³, …; QUADRATIC = z + z² (the cusp of M),
CUBIC = v + v² + v³/3 (the cusp z* = 1/√3, c₀ = 2/(3√3) of the z³ + c main component, rescaled by v = √3·w)."""
from __future__ import annotations
from fractions import Fraction as Fr
from functools import lru_cache
import mpmath as mp

QUADRATIC = (Fr(1),)
CUBIC = (Fr(1), Fr(1, 3))


def _smul(a, b, L):
    return [sum((a[i] * b[k - i] for i in range(k + 1)), Fr(0)) for k in range(L + 1)]


def _sinv(a, L):
    """1/a for a power series with a[0] ≠ 0."""
    out = [Fr(1) / a[0]] + [Fr(0)] * L
    for k in range(1, L + 1):
        out[k] = -sum((a[j] * out[k - j] for j in range(1, k + 1)), Fr(0)) / a[0]
    return out


def _slog1p(x, L):
    """log(1 + x) for x[0] = 0."""
    out, pw = [Fr(0)] * (L + 1), [Fr(1)] + [Fr(0)] * L
    for j in range(1, L + 1):
        pw = _smul(pw, x, L)
        out = [o + Fr((-1) ** (j + 1), j) * t for o, t in zip(out, pw)]
    return out


@lru_cache(maxsize=None)
def fatou_series(germ: tuple = QUADRATIC, K: int = 30):
    """(β, [c_0 = 0, c_1, …, c_K]) with Φ(w) = w − β log w + Σ c_k w^{−k} solving Φ∘F = Φ + 1."""
    L = K + 2
    g = [Fr(0), Fr(1)] + list(germ) + [Fr(0)] * L
    U = [-g[j] * (-1) ** j for j in range(L + 1)]            # U(u) = 1/F = −g(−u)
    V = U[1:] + [Fr(0)]                                        # V = U/u = 1 − u + …
    beta = 1 - g[3]
    invV = _sinv(V, L)
    base = [a - (1 if k == 0 else 0) for k, a in enumerate(invV)]             # 1/V − 1 = u·(…)
    lhs = [Fr(0)] * (L + 1)
    for k in range(L):                                         # (1/U − 1/u) = (1/V − 1)/u
        lhs[k] = base[k + 1]
    lhs[0] -= 1                                                # − 1
    logV = _slog1p([a - (1 if k == 0 else 0) for k, a in enumerate(V)], L)
    lhs = [x + beta * y for x, y in zip(lhs, logV)]
    assert lhs[0] == 0 and lhs[1] == 0, "β must cancel the order-u term"
    c = [Fr(0)] * (K + 1)
    Upow = [Fr(1)] + [Fr(0)] * L
    powers = []
    for k in range(1, K + 1):
        Upow = _smul(Upow, U, L); powers.append(Upow)
    for m in range(2, K + 2):
        s = lhs[m] + sum((c[k] * (powers[k - 1][m] - (1 if m == k else 0)) for k in range(1, m - 1)), Fr(0))
        c[m - 1] = s / (m - 1)                                 # coefficient of c_{m−1} at u^m is −(m−1)
    return beta, c


def fatou_coeffs(K: int, germ: tuple = QUADRATIC) -> list:
    return fatou_series(germ, K)[1]


def horn_coeffs(M, h=0.25, N=64, R=400, n=500, dps=40, K=30, germ: tuple = QUADRATIC):
    """Upper horn-map Fourier coefficients a_0..a_M (mpmath) and the aliasing level |ĉ_{N/2}|."""
    beta_f, cf = fatou_series(germ, K)
    with mp.workdps(dps):
        beta = mp.mpf(beta_f.numerator) / beta_f.denominator
        cm = [mp.mpf(x.numerator) / x.denominator for x in cf]
        gm = [mp.mpf(x.numerator) / x.denominator for x in germ]
        g = lambda v: v + mp.fsum(gm[j] * v ** (j + 2) for j in range(len(gm)))
        Fm = lambda w: -1 / g(-1 / w)
        def ph(w, rep):
            lg = mp.log(w) if not rep else mp.log(abs(w)) + 1j * (mp.arg(w) % (2 * mp.pi))
            return w - beta * lg + mp.fsum(cm[k] * w ** (-k) for k in range(1, K + 1))
        dph = lambda w: 1 - beta / w - mp.fsum(k * cm[k] * w ** (-k - 1) for k in range(1, K + 1))
        def E(Z):
            Zp = Z - n; w = Zp + beta * mp.log(-Zp)
            for _ in range(80):
                s = (ph(w, True) - Zp) / dph(w); w -= s
                if abs(s) < mp.mpf(10) ** (-dps + 5): break
            for _ in range(n): w = Fm(w)
            m = 0
            while abs(w) < R or w.real < abs(w) / 2:
                w = Fm(w); m += 1
            return ph(w, False) - m
        Zs = [mp.mpf(j) / N + 1j * mp.mpf(h) for j in range(N)]
        D = [E(Z) - Z for Z in Zs]
        coef = lambda k: mp.fsum(D[j] * mp.expjpi(-2 * mp.mpf(k) * j / N) for j in range(N)) / N
        return [coef(k) * mp.exp(2 * mp.pi * k * h) for k in range(M + 1)], abs(coef(N // 2))


def kappa0_mp(h=0.25, N=48, R=400, n=500, dps=40, K=30, germ: tuple = QUADRATIC):
    """(a₁, a₂, κ₀ = a₂/(2πi a₁²), aliasing) for the upper horn map of `germ`."""
    a, al = horn_coeffs(2, h=h, N=N, R=R, n=n, dps=dps, K=K, germ=germ)
    with mp.workdps(dps):
        return a[1], a[2], a[2] / (2j * mp.pi * a[1] ** 2), al


if __name__ == "__main__":
    for name, gm in (("z + z²", QUADRATIC), ("v + v² + v³/3", CUBIC)):
        a, al = horn_coeffs(3, dps=30, germ=gm)
        with mp.workdps(30):
            print(f"{name}: a₀ = {mp.nstr(a[0], 6)}, a₁ = {mp.nstr(a[1], 15)}, a₂ = {mp.nstr(a[2], 15)}, "
                  f"κ₀ = {mp.nstr(a[2] / (2j * mp.pi * a[1] ** 2), 15)}, aliasing {mp.nstr(al, 2)}")
