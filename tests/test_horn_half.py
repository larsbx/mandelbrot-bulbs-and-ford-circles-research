from fractions import Fraction as Fr
import mpmath as mp
from bulbford.horn_half import fatou_series_half, horn_coeffs_half, kappa_half


def test_fatou_series_solves_abel_equation():
    """Φ(g(z)) − Φ(z) = 1 + O(z^{K}) for g = f∘f, f(z) = −z + z² (formal, exact rationals); Φ ~ 1/(4z²)."""
    d, beta = fatou_series_half(12)
    assert d[-2] == Fr(1, 4)
    with mp.workdps(50):
        z = mp.mpf("0.01") * mp.expjpi(mp.mpf("0.1"))
        g = lambda z: z - 2 * z ** 3 + z ** 4
        Phi = lambda z: sum(mp.mpf(c.numerator) / c.denominator * z ** k for k, c in d.items()) + \
            mp.mpf(beta.numerator) / beta.denominator * mp.log(z)
        assert abs(Phi(g(z)) - Phi(z) - 1) < mp.mpf(10) ** -24


def test_horn_map_half_is_stable_in_height():
    """a₂/a₁² of the upper horn map (above the centre β·π/2 of the repelling cylinder) does not depend on the height."""
    a, _ = horn_coeffs_half(2, h=4.5, dps=60)
    b, _ = horn_coeffs_half(2, h=5.0, dps=60)
    with mp.workdps(60):
        assert abs(a[2] / a[1] ** 2 - b[2] / b[1] ** 2) / abs(a[2] / a[1] ** 2) < 1e-20


def test_half_root_constant_matches_bulbs():
    """V49: ι(𝒫(f_½)) − ½ equals the bulb limit at the ½-root (disc 1/q and cardioid [0;1,1,N], V48) to < 1e-8."""
    k = complex(kappa_half(h=4.5, dps=60)[2])
    assert abs(k - (0.018405261733 - 0.042903723571j)) < 1e-8
    assert abs(k - (0.018405261699 - 0.042903722127j)) < 1e-8


def test_half_root_unfolding_predicts_disc_bulbs():
    """V49: the renormalized family e^{u}𝒫_½ gives the disc's bulb size G(1/q) → 1.0201625870 and, with the phase
    μ = −1, the p = 2 class κ (V48 table) — the bounded-p theory transfers to the ½-root germ."""
    from bulbford.horn_half import HALF
    from bulbford.renorm import bounded_p_limit, G_limit
    assert abs(G_limit(1, 1, germ=HALF) - 1.0201625870) < 1e-8
    assert abs(complex(bounded_p_limit(2, 1, dps=60, M=20, germ=HALF)) - (0.0475525586 - 0.0197495608j)) < 1e-8
