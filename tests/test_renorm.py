import mpmath as mp
import numpy as np
import pytest
from bulbford.normal_form import normal_form
from bulbford.renorm import bounded_p_limit, G_limit, G_limits, rho_path, unfolding_taylor, lavaurs_compose


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


def test_lavaurs_compose_is_series_of_rho_after_moebius():
    """[u^k] ρ(u/(1 + u/(2πipq))) for a random polynomial ρ, against direct evaluation on a small circle."""
    rng = np.random.default_rng(1)
    rho = rng.normal(size=8) + 1j * rng.normal(size=8)
    p, q, R, n = 2, 5, 0.05, 64
    us = R * np.exp(2j * np.pi * np.arange(n) / n)
    vals = np.polyval(rho[::-1], us / (1 + us / (2j * np.pi * p * q)))
    direct = np.fft.fft(vals)[:8] / n / R ** np.arange(8)
    assert np.allclose(lavaurs_compose(rho, p, q), direct, atol=1e-10)


def test_unfolding_taylor_p1():
    """Cauchy coefficients of ρ for M_u = e^{u}𝒫₀: r₀ = 1, r₁ = −1 (P0, diagnostic of the radius, E-h), r₂ = κ₀."""
    c = unfolding_taylor(1, 1)
    assert abs(c[0] - 1) < 1e-12 and abs(c[1] + 1) < 1e-12
    assert abs(c[2] - (0.023825887402200569 - 0.052304659114003666j)) < 1e-12


@pytest.mark.parametrize("p,r", [(1, 1), (2, 1)])
def test_first_order_correction_is_kinematic(p, r):
    """R_{p/q}(u) = ρ_{p,r}(u/(1 + u/(2πipq))) + O(q⁻²) (register V47): the residual of r₂…r₄ falls by 4 per doubling
    of q, while the plain ρ leaves an O(1/q) error (factor ≤ 2 asymptotically)."""
    from bulbford.taylor import taylor
    rho = unfolding_taylor(p, r)[:5]
    res, naive = {}, {}
    for N in (64, 128):
        q = p * N + r
        c = taylor(p, q, r=1.5, N=128).coeffs[:5]
        res[N] = np.abs(c - np.array(lavaurs_compose(rho, p, q)))[2:]
        naive[N] = np.abs(c - rho)[2:]
    assert np.all((res[64] / res[128] > 3.7) & (res[64] / res[128] < 4.3))
    assert np.all(naive[64] / naive[128] < 3)                                   # not O(q⁻²)
