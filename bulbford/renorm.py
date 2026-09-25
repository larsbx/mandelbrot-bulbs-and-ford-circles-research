"""Renormalization predictions from the upper horn map 𝒫₀ of z + z² (V38–V39, C23).

For fixed p and q = pN + r → ∞ (r coprime to p), the first near-parabolic renormalization of f_{p/q} tends to
G = μ·𝒫₀ with μ = e^{−2πi/α} = e^{−2πi r/p}, and κ(p/q) → (ι(G^p) − ½)/p.  p = 1 gives κ₀."""
from __future__ import annotations
from functools import lru_cache
import mpmath as mp
from .horn import horn_coeffs
from .germ import horn_germ, iterate_index


@lru_cache(maxsize=None)
def _coeffs(M: int, dps: int):
    return horn_coeffs(M, h=0.25, N=64, dps=dps)[0]


def bounded_p_limit(p: int, r: int, dps: int = 40, M: int = 14):
    """Predicted lim κ(p/(pN + r)) as N → ∞."""
    with mp.workdps(dps):
        G = horn_germ(_coeffs(M, dps), 2 * p + 1, mp.expjpi(-2 * mp.mpf(r) / p))
        return (iterate_index(G, p)[0] - mp.mpf(1) / 2) / p
