import mpmath as mp
import pytest
from flint import ctx
from bulbford.normal_form import normal_form
from bulbford.index import index_complex, _series_fq


@pytest.mark.parametrize("p,q", [(1, 2), (1, 3), (2, 5), (3, 7), (4, 13), (10, 23), (26, 59)])
def test_P4_index_from_normal_form(p, q):
    """qι = (q²−1)/2 + b₂/b₁²  (P4)."""
    assert abs(complex(normal_form(p, q).iota) - index_complex(p, q)) < 1e-12


@pytest.mark.parametrize("p,q", [(1, 3), (3, 7), (10, 23)])
def test_P4_leading_coefficient(p, q):
    """a = [z^{q+1}] f^q = q b₁ = q λ₀⁻¹ Σ_{i+j=q+1} H_i H_j  (P4)."""
    old, oc = ctx.dps, ctx.cap
    ctx.dps, ctx.cap = 2 * q + 40, 2 * q + 2
    try:
        a = complex(_series_fq(p, q, 2 * q + 2).coeffs()[q + 1])
    finally:
        ctx.dps, ctx.cap = old, oc
    assert abs(complex(normal_form(p, q).a) / a - 1) < 1e-12


def test_normal_form_is_precision_stable():
    """No cancellation: κ agrees to 30 digits at 50 and 120 working digits (q = 211)."""
    k = lambda dps: (lambda nf: (nf.iota - mp.mpf(1) / 2) / 211)(normal_form(5, 211, dps=dps))
    with mp.workdps(120):
        assert abs(k(50) - k(120)) < mp.mpf(10) ** -30
