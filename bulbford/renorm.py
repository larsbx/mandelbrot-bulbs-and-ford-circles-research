"""Renormalization predictions from the upper horn map 𝒫₀ of z + z² (V38–V39, C23).

For fixed p and q = pN + r → ∞ (r coprime to p), the first near-parabolic renormalization of f_{p/q} tends to
G = μ·𝒫₀ with μ = e^{−2πi/α} = e^{−2πi r/p}, and κ(p/q) → (ι(G^p) − ½)/p.  p = 1 gives κ₀."""
from __future__ import annotations
from functools import lru_cache
import mpmath as mp
from .horn import horn_coeffs
from .germ import horn_germ, iterate_index, compose


@lru_cache(maxsize=None)
def _coeffs(M: int, dps: int):
    return horn_coeffs(M, h=0.25, N=64, dps=dps)[0]


def bounded_p_limit(p: int, r: int, dps: int = 40, M: int = 14):
    """Predicted lim κ(p/(pN + r)) as N → ∞."""
    with mp.workdps(dps):
        G = horn_germ(_coeffs(M, dps), 2 * p + 1, mp.expjpi(-2 * mp.mpf(r) / p))
        return (iterate_index(G, p)[0] - mp.mpf(1) / 2) / p


def unfolding(p: int, r: int, M: int = 20, dps: int = 60):
    """C23 unfolding M_u(W) = μ e^{u/p²} 𝒫₀(W), μ = e^{−2πir/p} (numpy, double precision).
    Returns (cycle, seed): cycle(u, W) → (W', ρ) refines a point of the p-cycle bifurcating from 0 and returns its
    multiplier; seed(u) is its leading-order position W^p = −u/(p·A_p), A_p = [W^{p+1}]((μ𝒫₀)^p)."""
    import numpy as np
    a = np.array([complex(x) for x in _coeffs(M, dps)])   # dps 60: a_n noise·|W|^n < 1e-16 for |W| ≲ 10, n ≤ 20
    n = np.arange(1, M + 1)
    mu = np.exp(-2j * np.pi * r / p)
    P = lambda W: W * np.exp(2j * np.pi * np.sum(a[1:] * W ** n))
    dP = lambda W: np.exp(2j * np.pi * np.sum(a[1:] * W ** n)) * (1 + 2j * np.pi * np.sum(n * a[1:] * W ** n))
    with mp.workdps(30):
        G = horn_germ(_coeffs(M, dps), p + 1, mp.expjpi(-2 * mp.mpf(r) / p))
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


def rho_path(p: int, r: int, us):
    """ρ along a path us (starting near 0), by continuation of the bifurcating p-cycle."""
    import numpy as np
    cycle, seed = unfolding(p, r)
    W, out = seed(us[0]), []
    for u in us:
        W, d = cycle(u, W); out.append(d)
    return np.array(out), W


def G_limit(p: int, r: int) -> float:
    """C23 prediction for lim G(p/(pN + r)) = |u_a|/2 with ρ(u_a) = −1 (continuation 0 → 2, then Newton in u)."""
    import numpy as np
    cycle, seed = unfolding(p, r)
    us = np.linspace(0.02, 2.0, 400) + 0j
    W = seed(us[0])
    for u in us:
        W, _ = cycle(u, W)
    u, h = us[-1], 1e-7
    for _ in range(50):
        W, d0 = cycle(u, W)
        _, d1 = cycle(u + h, W)
        du = (d0 + 1) / ((d1 - d0) / h); u -= du
        if abs(du) < 1e-13: break
    return abs(u) / 2
