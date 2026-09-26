import mpmath as mp
import pytest
from bulbford.horn_rational import fatou_series_pq, horn_index_pq

K0 = complex(0.023825887402200569, -0.052304659114003666)
KHALF = complex(0.018405261617964781, -0.042903723880230475)


@pytest.mark.parametrize("p,q", [(0, 1), (1, 2), (1, 3), (2, 5)])
def test_fatou_series_solves_abel_equation(p, q):
    """Φ(f^q(z)) − Φ(z) = 1 + O(z^K) for the formal Fatou coordinate of the q-petal germ f = e^{2πip/q}z + z²."""
    with mp.workdps(40):
        s, beta, g = fatou_series_pq(p, q, K=14)
        z = mp.mpf("0.004") * mp.expjpi(mp.mpf("0.37"))
        Phi = lambda z: mp.fsum(c * z ** k for k, c in s.items()) + beta * mp.log(z)
        assert abs(Phi(g(z)) - Phi(z) - 1) < mp.mpf(10) ** -18


@pytest.mark.parametrize("p,q,upper,expected", [(0, 1, True, K0), (0, 1, False, K0.conjugate()), (1, 2, True, KHALF)])
def test_horn_index_reproduces_known_roots(p, q, upper, expected):
    """The general construction reproduces κ₀ (cusp, both horn maps) and κ_½ (V38, V49)."""
    assert abs(complex(horn_index_pq(p, q, upper=upper, dps=40)) - expected) < 1e-15


@pytest.mark.parametrize("upper,bulbs", [(True, 0.016581162601 - 0.039504714589j), (False, 0.015840773241 + 0.040175221243j)])
def test_one_sided_limits_at_one_third(upper, bulbs):
    """V50: the upper/lower horn index of e^{2πi/3}z + z² equals the limit of κ along [0;2,1,N] (from above) /
    [0;3,N] (from below), N → ∞ (bulb fits of degree 4, N ≤ 2048); the two sides are not conjugate."""
    assert abs(complex(horn_index_pq(1, 3, upper=upper, dps=30)) - bulbs) < 1e-9


def test_polynomial_germ_spec_matches_quadratic():
    """poly_germ([λ₀, 1]) is the same germ as the default e^{2πi/3}z + z²; two_cycle_germ has multiplier 4(c+1)."""
    from bulbford.horn_rational import poly_germ, two_cycle_germ, _poly_coeffs
    with mp.workdps(40):
        spec = poly_germ([mp.expjpi(mp.mpf(2) / 3), mp.mpc(1)])
        assert abs(horn_index_pq(1, 3, dps=20, poly=spec) - horn_index_pq(1, 3, dps=20)) < 1e-15
        c = -1 + mp.expjpi(mp.mpf(2) / 3) / 4
        assert abs(_poly_coeffs(1, 3, two_cycle_germ(c))[0] - 4 * (c + 1)) < 1e-30


def test_lower_tail_uses_conjugate_phase():
    """V55 (referee A1): along [0;2,N,3] → ½⁻ the tail phase is μ̄_t, not μ_t.  1 − [0;2,N,3] = [0;1,1,N,3] → ½⁺, so
    the lower prediction is the conjugate of the upper one; the bulbs give κ([0;2,256,3]) = 0.0534760 + 0.0026281i
    (q = 1541, O(1/q) away from the limit), the rule with μ_t would give 0.0537821 + 0.0209937i."""
    from bulbford.renorm import bounded_p_limit
    from bulbford.horn_rational import root
    lo, up = (complex(bounded_p_limit(3, 1, dps=60, M=20, germ=root(1, 2, u))) for u in (False, True))
    assert abs(lo - up.conjugate()) < 1e-12
    assert abs(lo - (0.0533070 + 0.0026060j)) < 1e-7


def test_deep_lower_cusp_germ_is_conjugate_of_cusp():
    """V55: at n = 3000 the lower map of z + z² (as the root 1/1) needs an escape radius independent of the depth; it
    then gives the conjugate of the cusp prediction for the tail 3 (1 − [0;1,N,3] = [0;N+1,3])."""
    from bulbford.renorm import bounded_p_limit
    from bulbford.horn_rational import root
    lo = complex(bounded_p_limit(3, 1, dps=60, M=20, germ=root(1, 1, False)))
    assert abs(lo - complex(bounded_p_limit(3, 1, dps=60, M=20)).conjugate()) < 1e-15
