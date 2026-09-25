import mpmath as mp
import numpy as np
import pytest
from bulbford.normal_form import normal_form
from bulbford.renorm import bounded_p_limit, G_limit, rho_path


@pytest.mark.parametrize("p,r", [(2, 1), (3, 1), (3, 2)])
def test_V39_bounded_p_kappa(p, r):
    """κ(p/(pN+r)) + i/(2πpq) → (ι((μ𝒫₀)^p) − ½)/p, μ = e^{−2πir/p}  (C23 (ii), V39)."""
    q = p * 256 + r
    with mp.workdps(40):
        nf = normal_form(p, q, dps=40)
        k = (nf.iota - mp.mpf(1) / 2) / q + 1j / (2 * mp.pi * p * q)
        assert abs(k - bounded_p_limit(p, r)) < 10.0 / q ** 2


def test_V40_unfolding_taylor_coefficients():
    """ρ(u) of the p-cycle of e^{u}𝒫₀ has r₀, r₁, r₂ = 1, −1, κ₀ (P0, V38)."""
    us = np.linspace(0.02, 0.3, 30)
    v, _ = rho_path(1, 1, us)
    c = np.polyfit(us, v, 8)[::-1]
    assert abs(c[0] - 1) < 1e-8 and abs(c[1] + 1) < 1e-6
    assert abs(c[2] - (0.0238258874 - 0.0523046591j)) < 1e-4


def test_V40_G_limit_p1():
    """lim G(1/q) = 1.026411… (V24 at q = 16001: 1.026411) from the horn-map family alone."""
    assert abs(G_limit(1, 1) - 1.026411) < 1e-5
