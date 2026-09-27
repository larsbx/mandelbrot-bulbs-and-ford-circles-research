"""Taylor coefficients of R_q(u) := ρ(c_of(λ₀ e^{u/q²})) at u = 0, by Cauchy/FFT on |u| = r."""
from __future__ import annotations
from typing import NamedTuple
import numpy as np
from .dynamics import Family, MAIN2, rho_on_path


class TaylorData(NamedTuple):
    p: int
    q: int
    r: float
    coeffs: np.ndarray        # r_k = [u^k] R_q(u), k = 0..N-1

    def R(self, u: complex, nterms: int | None = None) -> complex:
        return complex(np.polyval(self.coeffs[:nterms][::-1], u))

    def solve(self, target: complex, u0: complex, nterms: int | None = None) -> complex:
        """Newton on R(u) = target from u0 (truncated series)."""
        c = self.coeffs[:nterms]
        dc = c[1:] * np.arange(1, len(c))
        u = u0
        for _ in range(100):
            s = (np.polyval(c[::-1], u) - target) / np.polyval(dc[::-1], u)
            u -= s
            if abs(s) < 1e-14:
                break
        return complex(u)


def taylor(p: int, q: int, fam: Family = MAIN2, r: float = 2.5, N: int = 256) -> TaylorData:
    us = r * np.exp(2j * np.pi * np.arange(N) / N)
    vals = rho_on_path(fam, p, q, us)
    coeffs = np.fft.fft(vals) / N / r ** np.arange(N)
    return TaylorData(p, q, r, coeffs)


def kappa_fft(p: int, q: int, fam: Family = MAIN2, r: float = 2.5, N: int = 64) -> TaylorData:
    """Cheap Taylor data (N=64 suffices: |r_k| decays like q^{-k})."""
    return taylor(p, q, fam, r=r, N=N)
