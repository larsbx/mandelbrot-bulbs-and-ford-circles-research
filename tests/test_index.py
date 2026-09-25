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


def test_series_reconstructs_diameter_to_second_order():
    """|u_a|/2 from the Taylor data at u=0 equals G_ant up to O(q⁻²) (nonlinearity of c_of)."""
    from bulbford.dynamics import MAIN2, bulb
    for p, q in [(1, 23), (7, 23), (10, 23)]:
        t = taylor(p, q, r=2.5, N=128)
        assert abs(abs(t.solve(-1, 2.0)) / 2 - bulb(MAIN2, p, q).G_ant) < 5 / q ** 2


@pytest.mark.parametrize("fam_name,p,q", [("main2", 2, 7), ("disk2", 2, 7), ("disk2", 3, 11), ("main3", 2, 7), ("main3", 3, 11)])
def test_P2_other_families(fam_name, p, q):
    """P2 is family-generic: [u²]R_q = (ι − ½)/q with ι the index of the parabolic cycle point."""
    from bulbford.dynamics import FAMILIES
    from bulbford.index import index_at
    fam = FAMILIES[fam_name]
    t = taylor(p, q, fam, r=1.0, N=64)
    assert abs(t.coeffs[2] - (index_at(fam, p, q) - 0.5) / q) < 1e-7
