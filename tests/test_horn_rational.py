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
