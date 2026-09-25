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


def test_germ_normal_form_matches_quadratic_P4():
    """normal_form_index on the germ λ₀z + z² reproduces P4's ι for f_{λ₀}."""
    from bulbford.germ import normal_form_index
    with mp.workdps(40):
        for p, q in [(1, 3), (3, 7), (4, 13)]:
            lam = mp.expjpi(mp.mpf(2 * p) / q)
            g = [mp.mpc(0), lam, mp.mpc(1)]
            assert abs(normal_form_index(g, q) - normal_form(p, q, dps=40).iota) < mp.mpf(10) ** -30


@pytest.mark.parametrize("p,q,B", [(1, 7, 4), (3, 23, 8), (26, 59, 16), (5, 211, 7), (1, 257, None)])
def test_ball_normal_form_encloses_mpmath(p, q, B):
    """Certified κ (Arb balls, blocked convolution) contains the 80-digit mpmath value, for several block sizes."""
    from bulbford.normal_form_ball import normal_form_ball
    nb = normal_form_ball(p, q, B=B)
    with mp.workdps(80):
        km = (normal_form(p, q, dps=80).iota - mp.mpf(1) / 2) / q
        mid = mp.mpc(nb.kappa.real.mid().str(75, radius=False), nb.kappa.imag.mid().str(75, radius=False))
        assert float(nb.kappa.rad()) < 1e-40
        assert abs(mid - km) <= 2 * mp.mpf(float(nb.kappa.rad())) + mp.mpf(10) ** -70


def test_ball_normal_form_leading_coefficient():
    """a = q·b₁ from the ball normal form equals [z^{q+1}] f^q (series route)."""
    from bulbford.normal_form_ball import normal_form_ball
    old, oc = ctx.dps, ctx.cap
    ctx.dps, ctx.cap = 2 * 23 + 40, 2 * 23 + 2
    try:
        a = complex(_series_fq(10, 23, 2 * 23 + 2).coeffs()[23 + 1])
    finally:
        ctx.dps, ctx.cap = old, oc
    assert abs(complex(normal_form_ball(10, 23).a) / a - 1) < 1e-12
