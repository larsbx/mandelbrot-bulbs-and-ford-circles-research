import numpy as np
import pytest
from bulbford.index import index, index_complex, kappa
from bulbford.taylor import taylor


def test_index_q1_q2_exact():
    assert abs(index_complex(1, 1) - 0.0) < 1e-20          # z+z²: no z³ term
    assert abs(index_complex(1, 2) - 0.125) < 1e-20        # z−2z³+z⁴ → 1/8


def test_index_symmetry():
    """ι(−p/q) = conj ι(p/q) (complex conjugation of the family)."""
    for p, q in [(1, 5), (2, 7), (3, 11)]:
        assert abs(index_complex(q - p, q) - np.conj(index_complex(p, q))) < 1e-12


@pytest.mark.parametrize("p,q", [(1, 3), (2, 5), (3, 7), (4, 13)])
def test_P2_second_order_coefficient(p, q):
    """[u^2] R_q = (ι − ½)/q;  [u^0] = 1;  [u^1] = −1  (P0, P2)."""
    t = taylor(p, q, r=1.0, N=64)
    assert abs(t.coeffs[0] - 1) < 1e-9
    assert abs(t.coeffs[1] + 1) < 1e-9
    assert abs(t.coeffs[2] - kappa(p, q)) < 1e-8
