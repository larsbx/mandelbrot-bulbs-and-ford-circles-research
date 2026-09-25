import mpmath as mp
from bulbford.horn import fatou_coeffs, kappa0_mp
from bulbford.normal_form import normal_form
from fractions import Fraction as Fr

with mp.workdps(40):     # parse at full precision (the default 15 digits would round K0 at 2e-18)
    K0 = mp.mpc("0.023825887402200569000729189645402", "-0.052304659114003666316522486207571")


def test_fatou_coefficients():
    """Φ(w) = w − log w + Σ c_k w^{−k} for F(w) = w²/(w−1) (f = z + z² in w = −1/z)."""
    assert fatou_coeffs(4)[1:] == [Fr(1, 2), Fr(1, 3), Fr(13, 36), Fr(113, 240)]


def test_V38_kappa0_is_horn_map_index():
    """κ₀ := lim κ(1/q) = a₂/(2πi a₁²) for the upper horn map of z + z²: the horn map side to 25 digits."""
    _, _, k, _ = kappa0_mp(h=0.25, N=32, R=300, n=400, dps=32, K=30)
    with mp.workdps(32):
        assert abs(k - K0) < mp.mpf(10) ** -25


def test_V37_bulb_side_approaches_horn_value():
    """κ(1/q) + i/(2πq) = κ₀ + O(q⁻²) with |c₂| ≈ 2.93 (V37)."""
    q = 512
    with mp.workdps(40):
        nf = normal_form(1, q, dps=40)
        k = (nf.iota - mp.mpf(1) / 2) / q + 1j / (2 * mp.pi * q)
        assert abs(k - K0) < 3.0 / q ** 2
