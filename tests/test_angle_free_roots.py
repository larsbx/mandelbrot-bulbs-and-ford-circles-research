"""The certified lane builds roots of unity from their polynomial, never from an angle.

Two independent angle-free certificates of zeta_q must meet: the Arb ball that
`index.py` selects from the roots of Phi_q, and the Krawczyk box of `antipode.py`.
The untrusted float seeds of the Krawczyk boxes come from the exact cosine brackets
of `spread.py` and must land within reach of Newton's method.
"""
from __future__ import annotations

from math import gcd

import pytest
from flint import acb, arb, ctx, fmpq

from bulbford.antipode import _unity_seed, lambda_box, unity_boxes, zeta_box
from bulbford.index import lambda_ball, zeta_ball


def as_arb(i) -> arb:
    """The rational interval [lo, hi] as an Arb ball containing it."""
    lo, hi = arb(fmpq(i.lo.numerator, i.lo.denominator)), arb(fmpq(i.hi.numerator, i.hi.denominator))
    return lo.union(hi)


def meets(ball: acb, box) -> bool:
    return ball.real.overlaps(as_arb(box.re)) and ball.imag.overlaps(as_arb(box.im))


@pytest.mark.parametrize("q", range(1, 17))
def test_zeta_ball_meets_the_krawczyk_box(q):
    old = ctx.prec
    try:
        ctx.prec = 200
        assert meets(zeta_ball(q, 200), zeta_box(q))
        for p in range(1, q + 1):
            if gcd(p, q) == 1:
                assert meets(lambda_ball(p, q, 200), lambda_box(p, q))
    finally:
        ctx.prec = old


def test_zeta_ball_is_tight_and_primitive():
    old = ctx.prec
    try:
        ctx.prec = 300
        z = zeta_ball(7, 300)
        assert z.rad() < arb(2) ** -250
        assert (z**7 - 1).contains(acb(0))
        assert not (z - 1).contains(acb(0))
        assert z.imag > 0
    finally:
        ctx.prec = old


@pytest.mark.parametrize("q", [3, 5, 8, 13, 31, 64])
def test_unity_seeds_are_close_to_the_certified_roots(q):
    boxes = unity_boxes(q)
    for k in range(q):
        seed, box = _unity_seed(k, q), boxes[k]
        assert abs(seed.real - float(box.re.lo)) < 1e-12
        assert abs(seed.imag - float(box.im.lo)) < 1e-12
