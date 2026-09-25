import mpmath as mp
import numpy as np
import pytest
from bulbford.normal_form import normal_form
from bulbford.renorm import bounded_p_limit, G_limit, G_limits, rho_path


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


def test_P6_degree3_bulbs_have_two_antipodes():
    from bulbford.dynamics import MAIN2, MAIN3, bulb
    assert len(bulb(MAIN2, 2, 7).ants) == 1
    b = bulb(MAIN3, 2, 7)
    assert len(b.ants) == 2 and all(abs(a.rho + 1) < 1e-10 for a in b.ants)
    assert abs(b.ants[0].c - b.ants[1].c) > 1e-3 * abs(b.cen.c - b.root)


def test_V41_degree3_kappa0_is_cubic_horn_index():
    """κ₃(1/q) + i/(2πq) → a₂/(2πi a₁²) of the horn map of v + v² + v³/3 (Cauchy radius 1, E-h)."""
    import math
    from bulbford.horn import CUBIC, kappa0_mp
    from bulbford.dynamics import MAIN3
    from bulbford.taylor import taylor
    k0 = complex(kappa0_mp(dps=30, germ=CUBIC)[2])
    t = taylor(1, 512, MAIN3, r=1.0, N=64)
    assert abs(t.coeffs[1] + 1) < 1e-8                       # P0 holds on the safe radius
    assert abs(t.coeffs[2] + 1j / (2 * math.pi * 512) - k0) < 2e-5


def test_V41_degree3_bounded_p_G_both_antipodes():
    from bulbford.horn import CUBIC
    from bulbford.dynamics import MAIN3, bulb
    pred = sorted(G_limits(3, 1, germ=CUBIC))
    got = sorted(bulb(MAIN3, 3, 3 * 512 + 1).G_ants)
    assert all(abs(a - b) < 2e-4 for a, b in zip(pred, got))
