"""Renormalization predictions from the upper horn map 𝒫₀ of z + z² (V38–V39, C23).

For fixed p and q = pN + r → ∞ (r coprime to p), the first near-parabolic renormalization of f_{p/q} tends to
G = μ·𝒫₀ with μ = e^{−2πi/α} = e^{−2πi r/p}, and κ(p/q) → (ι(G^p) − ½)/p.  p = 1 gives κ₀."""
from __future__ import annotations
from functools import lru_cache
import mpmath as mp
from .horn import horn_coeffs, QUADRATIC
from .germ import horn_germ, iterate_index, compose


@lru_cache(maxsize=None)
def _coeffs(M: int, dps: int, germ: tuple = QUADRATIC):
    return horn_coeffs(M, h=0.25, N=64, dps=dps, germ=germ)[0]


def bounded_p_limit(p: int, r: int, dps: int = 40, M: int = 14, germ: tuple = QUADRATIC):
    """Predicted lim κ(p/(pN + r)) as N → ∞."""
    with mp.workdps(dps):
        G = horn_germ(_coeffs(M, dps, germ), 2 * p + 1, mp.expjpi(-2 * mp.mpf(r) / p))
        return (iterate_index(G, p)[0] - mp.mpf(1) / 2) / p


def unfolding(p: int, r: int, M: int = 20, dps: int = 60, germ: tuple = QUADRATIC):
    """C23 unfolding M_u(W) = μ e^{u/p²} 𝒫₀(W), μ = e^{−2πir/p} (numpy, double precision).
    Returns (cycle, seed): cycle(u, W) → (W', ρ) refines a point of the p-cycle bifurcating from 0 and returns its
    multiplier; seed(u) is its leading-order position W^p = −u/(p·A_p), A_p = [W^{p+1}]((μ𝒫₀)^p)."""
    import numpy as np
    a = np.array([complex(x) for x in _coeffs(M, dps, germ)])   # dps 60: a_n noise·|W|^n < 1e-16 for |W| ≲ 10, n ≤ 20
    n = np.arange(1, M + 1)
    mu = np.exp(-2j * np.pi * r / p)
    P = lambda W: W * np.exp(2j * np.pi * np.sum(a[1:] * W ** n))
    dP = lambda W: np.exp(2j * np.pi * np.sum(a[1:] * W ** n)) * (1 + 2j * np.pi * np.sum(n * a[1:] * W ** n))
    with mp.workdps(30):
        G = horn_germ(_coeffs(M, dps, germ), p + 1, mp.expjpi(-2 * mp.mpf(r) / p))
        Gp = [mp.mpc(0), mp.mpc(1)] + [mp.mpc(0)] * p
        for _ in range(p):
            Gp = compose(G, Gp, p + 1)
        A_p = complex(Gp[p + 1])

    def orbit(u, W):
        c, z, d = mu * np.exp(u / p ** 2), W, 1 + 0j
        for _ in range(p):
            d *= c * dP(z); z = c * P(z)
        return z, d

    def cycle(u, W):
        for _ in range(100):                          # Newton on M^p(W)/W − 1 (root W = 0 deflated)
            z, d = orbit(u, W)
            s = (z / W - 1) / (d / W - z / W ** 2); W -= s
            if abs(s) < 1e-15 * abs(W): break
        return W, orbit(u, W)[1]

    seed = lambda u: (-u / (p * A_p)) ** (1 / p)
    return cycle, seed


def rho_path(p: int, r: int, us, germ: tuple = QUADRATIC):
    """ρ along a path us (starting near 0), by continuation of the bifurcating p-cycle."""
    import numpy as np
    cycle, seed = unfolding(p, r, germ=germ)
    W, out = seed(us[0]), []
    for u in us:
        W, d = cycle(u, W); out.append(d)
    return np.array(out), W


def unfolding_taylor(p: int, r: int, radius: float = 1.5, N: int = 128, germ: tuple = QUADRATIC):
    """Taylor coefficients ρ_0…ρ_{N−1} of the multiplier ρ_{p,r}(u) of M_u = μe^{u/p²}𝒫 (Cauchy/FFT on |u| = radius,
    reached by continuation along [radius/N, radius]; the analogue of taylor.taylor for the renormalized family)."""
    import numpy as np
    circle = radius * np.exp(2j * np.pi * np.arange(N) / N)
    v, _ = rho_path(p, r, np.concatenate([np.linspace(radius / N, radius, N), circle]), germ=germ)
    return np.fft.fft(v[N:]) / N / radius ** np.arange(N)


def lavaurs_compose(rho, p: int, q: int):
    """[u^k] ρ(ũ), ũ = u/(1 + u/(2πipq)): with λ = λ₀e^{u/q²}, α = p/q + u/(2πiq²), the Inou–Shishikura multiplier
    e^{−2πi/α} equals μ e^{ũ/p²} exactly, so R_{p/q}(u) = ρ_{p,r}(ũ(u)) up to the O(q⁻²) change of the germ.
    Uses ũ^m = Σ_j C(−m, j) u^{m+j} c^{−j}, c = 2πipq."""
    from math import comb
    c = 2j * mp.pi * p * q
    binom_neg = lambda m, j: (-1) ** j * comb(m + j - 1, j) if m else int(j == 0)
    return [sum(complex(rho[m] * binom_neg(m, k - m) / c ** (k - m)) for m in range(k + 1)) for k in range(len(rho))]


def G_limits(p: int, r: int, germ: tuple = QUADRATIC, steps: int = 400) -> tuple:
    """C23 predictions |u_a|/2 at every antipode: in the conformal coordinate s = ρ^{1/(d−1)} follow the segment
    from the root to each s with s^{d−1} = −1, continuing u by Newton on ρ(u) = s(t)^{d−1} (d − 1 = 2 for a germ
    with a cubic term: the centre is a double zero of ρ; d − 1 = 1 for z + z²); k = len(germ)."""
    import numpy as np
    k = len(germ)                                    # = deg g − 1 = order of the critical point = branching of ρ
    cycle, seed = unfolding(p, r, germ=germ)
    out = []
    for end in (np.exp(1j * np.pi * (2 * j + 1) / k) for j in range(k)):     # s with s^k = −1
        u = 0.02 + 0j
        W, rho = cycle(u, seed(u))
        s0 = rho ** (1 / k)
        for t in np.linspace(0, 1, steps + 1)[1:]:
            target = ((1 - t) * s0 + t * end) ** k
            for _ in range(30):
                W, d0 = cycle(u, W)
                h = 1e-7 * max(1.0, abs(u))
                _, d1 = cycle(u + h, W)
                du = (d0 - target) / ((d1 - d0) / h); u -= du
                if abs(du) < 1e-13 * max(1.0, abs(u)): break
        out.append(abs(u) / 2)
    return tuple(out)


def G_limit(p: int, r: int, germ: tuple = QUADRATIC) -> float:
    return G_limits(p, r, germ)[0]
