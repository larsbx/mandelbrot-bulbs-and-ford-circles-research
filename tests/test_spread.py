"""Laws of kernel/bulbford/spread.py: the bulb-Ford quantities in exact rational trigonometry.

No angle is evaluated anywhere below: cos(pi p/q) is the p-th largest root of the
integer polynomial U_{q-1}, isolated over Q by a Sturm sequence, and every bracket is
a pair of Fractions.
"""
from __future__ import annotations

from fractions import Fraction
from functools import reduce
from math import gcd

import pytest

from bulbford.antipode import ONE, lambda_box, quadrance
from bulbford.spread import (
    chebyshev_u,
    compose,
    ford_prediction_squared_bracket,
    multiplier_quadrance_bracket,
    poly_mul,
    roots_in,
    spread_polynomial,
    turn_cosine_bracket,
)

X = (0, 1)


def contains(bracket, x) -> bool:
    lo, hi = bracket
    return lo <= x <= hi


def strictly_contains_root(bracket, a: Fraction, b: Fraction, n: int) -> bool:
    """lo < a + b sqrt(n) < hi, decided by exact squaring (b > 0)."""
    lo, hi = bracket
    below = lambda r: r <= a or (r - a) ** 2 < b * b * n   # r < a + b sqrt(n)
    above = lambda r: r > a and (r - a) ** 2 > b * b * n   # r > a + b sqrt(n)
    return below(lo) and above(hi)


# --- integer polynomials -------------------------------------------------------------------------


def test_chebyshev_u_small():
    assert chebyshev_u(0) == (1,)
    assert chebyshev_u(1) == (0, 2)
    assert chebyshev_u(2) == (-1, 0, 4)
    assert chebyshev_u(3) == (0, -4, 0, 8)


def test_spread_polynomials_small_and_composition():
    assert spread_polynomial(0) == (0,)
    assert spread_polynomial(1) == (0, 1)
    assert spread_polynomial(2) == (0, 4, -4)                 # 4s(1-s)
    assert spread_polynomial(3) == (0, 9, -24, 16)            # s(3-4s)^2
    for n in range(5):
        for m in range(5):
            assert compose(spread_polynomial(n), spread_polynomial(m)) == spread_polynomial(n * m)


@pytest.mark.parametrize("n", range(1, 13))
def test_spread_is_s_times_u_squared(n):
    # sin(n t) = sin t U_{n-1}(cos t), so S_n(1 - c^2) = (1 - c^2) U_{n-1}(c)^2 identically in c.
    one_minus_c2 = (1, 0, -1)
    lhs = compose(spread_polynomial(n), one_minus_c2)
    u = chebyshev_u(n - 1)
    assert lhs == poly_mul(one_minus_c2, poly_mul(u, u))


def from_roots(lead: int, roots) -> tuple[int, ...]:
    """lead * prod (d x - n) for roots n/d, as integer coefficients."""
    out = (lead,)
    for r in roots:
        r = Fraction(r)
        out = poly_mul(out, (-r.numerator, r.denominator))
    return out


@pytest.mark.parametrize(
    "lead, roots",
    [
        (1, [-3, 1, 2]),
        (-1, [-3, 1, 2]),
        (-5, [Fraction(-7, 3), Fraction(1, 2), 4, 9]),
        (3, [Fraction(-1, 9), Fraction(1, 8), Fraction(1, 7), 2, 5]),
        (-2, [0, Fraction(1, 1000), Fraction(2, 1000), Fraction(3, 1000)]),
    ],
)
def test_sturm_counts_any_squarefree_polynomial(lead, roots):
    # Leading coefficients of either sign, so the sign of each pseudo-remainder matters.
    p = from_roots(lead, roots)
    grid = sorted({Fraction(r) for r in roots} | {Fraction(-10), Fraction(10)})
    probes = sorted(set(grid) | {(a + b) / 2 for a, b in zip(grid, grid[1:])})
    for lo in probes:
        for hi in probes:
            if lo < hi:
                assert roots_in(p, lo, hi) == sum(1 for r in roots if lo < r <= hi)
    # A quadratic factor with no real roots changes nothing.
    assert roots_in(poly_mul(p, (1, 0, 1)), Fraction(-10), Fraction(10)) == len(roots)


# --- brackets ------------------------------------------------------------------------------------


def test_brackets_are_exact_and_narrow():
    for bits in (8, 64, 200):
        lo, hi = multiplier_quadrance_bracket(1, 7, bits)
        assert isinstance(lo, Fraction) and isinstance(hi, Fraction)
        assert lo <= hi and hi - lo <= Fraction(1, 2**bits) * 8


def test_rational_cases():
    # Qd(1 - lambda0) at the rational spreads: q = 2, 3, 4, 6.
    assert contains(multiplier_quadrance_bracket(1, 2), 4)
    assert contains(multiplier_quadrance_bracket(1, 3), 3)
    assert contains(multiplier_quadrance_bracket(1, 4), 2)
    assert contains(multiplier_quadrance_bracket(1, 6), 1)
    assert contains(multiplier_quadrance_bracket(2, 6), 3)
    assert contains(turn_cosine_bracket(1, 3), Fraction(1, 2))
    assert contains(turn_cosine_bracket(1, 2), 0)


def test_quadratic_irrational_cases():
    # Qd(1 - zeta_12) = 2 - sqrt 3; Qd(1 - zeta_5) = (5 - sqrt 5)/2; Qd(1 - zeta_8) = 2 - sqrt 2.
    neg = lambda br: (-br[1], -br[0])
    assert strictly_contains_root(neg(multiplier_quadrance_bracket(1, 12)), Fraction(-2), Fraction(1), 3)
    assert strictly_contains_root(neg(multiplier_quadrance_bracket(1, 5)), Fraction(-5, 2), Fraction(1, 2), 5)
    assert strictly_contains_root(neg(multiplier_quadrance_bracket(1, 8)), Fraction(-2), Fraction(1), 2)
    assert strictly_contains_root(multiplier_quadrance_bracket(3, 8), Fraction(2), Fraction(1), 2)


@pytest.mark.parametrize("q", [5, 7, 12, 17])
def test_symmetry_and_order(q):
    brackets = [multiplier_quadrance_bracket(p, q, 80) for p in range(1, q)]
    for p in range(1, q):
        assert multiplier_quadrance_bracket(p, q, 80) == multiplier_quadrance_bracket(q - p, q, 80)
        assert multiplier_quadrance_bracket(p + q, q, 80) == brackets[p - 1]
    # Strictly increasing for 1 <= p <= q/2: the brackets are disjoint and ordered.
    for p in range(1, q // 2):
        assert brackets[p - 1][1] < brackets[p][0]


@pytest.mark.parametrize("q", [3, 5, 8, 9, 12])
def test_global_identities(q):
    # prod_{p=1}^{q-1} (1 - zeta^p) = q, so the quadrances multiply to q^2; and
    # sum_{p=0}^{q-1} Qd(1 - zeta^p) = sum (2 - zeta^p - zeta^-p) = 2q.
    brackets = [multiplier_quadrance_bracket(p, q, 120) for p in range(1, q)]
    lo = reduce(lambda a, b: a * b[0], brackets, Fraction(1))
    hi = reduce(lambda a, b: a * b[1], brackets, Fraction(1))
    assert lo <= q * q <= hi
    assert sum(b[0] for b in brackets) <= 2 * q <= sum(b[1] for b in brackets)


@pytest.mark.parametrize("q", range(2, 13))
def test_agrees_with_the_krawczyk_root_of_unity_boxes(q):
    # Two independent angle-free certificates of the same algebraic number must meet:
    # Sturm on U_{q-1} over Q, and Krawczyk boxes on X^q - 1 over Q(i).
    for p in range(1, q):
        if gcd(p, q) != 1:
            continue
        lo, hi = multiplier_quadrance_bracket(p, q, 90)
        box = quadrance(ONE - lambda_box(p, q))
        assert lo <= box.hi and box.lo <= hi


def test_ford_prediction_matches_p1():
    # P1: the period-2 disk has diameter 1/2 = 2 sin(pi/2)/4, so pred^2 = 1/4.
    assert contains(ford_prediction_squared_bracket(1, 2), Fraction(1, 4))
    # pred^2 = Qd(1 - lambda0)/q^4 bracket for bracket.
    lo, hi = multiplier_quadrance_bracket(2, 7, 70)
    assert ford_prediction_squared_bracket(2, 7, 70) == (lo / 7**4, hi / 7**4)


def test_refusals():
    for bad in [(0, 5), (5, 5), (1, 0), (1, -3)]:
        with pytest.raises(ValueError):
            multiplier_quadrance_bracket(*bad)


def test_moderate_denominator():
    lo, hi = multiplier_quadrance_bracket(1, 64, 64)
    assert 0 < lo <= hi < Fraction(1, 100)


@pytest.mark.parametrize("q", [3, 5, 8, 12])
def test_cosine_has_period_2q_and_is_even(q):
    neg = lambda br: (-br[1], -br[0])
    for p in range(-2 * q, 3 * q):
        c = turn_cosine_bracket(p, q, 60)
        assert turn_cosine_bracket(p + 2 * q, q, 60) == c
        assert turn_cosine_bracket(-p, q, 60) == c
        assert turn_cosine_bracket(p + q, q, 60) == neg(c)
    assert turn_cosine_bracket(0, q) == (1, 1)
    assert turn_cosine_bracket(q, q) == (-1, -1)
