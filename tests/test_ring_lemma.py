"""P8, the second-ring lemma, one test per step of its proof, in integers for 5 ≤ q ≤ 60."""

from math import gcd

import pytest

from bulbford.wake import heights, mechanical

FRACTIONS = [(p, q) for q in range(5, 61) for p in range(1, q) if gcd(p, q) == 1]


def shift(x: int, s: int, q: int) -> int:
    return (x * pow(2, s, 2**q - 1)) % (2**q - 1)


def normal(x: int, q: int) -> frozenset[int]:
    """The height set of x, translated so that its minimum is 0."""
    g = heights(x, q)
    return frozenset(v - min(g) for v in g)


def ring2(p: int, q: int) -> tuple[int, int]:
    M, w = 2**q - 1, p - 1
    return (mechanical(p, q, (w - 2) % q) + 1) % M, (mechanical(p, q, (w + 3) % q) - 1) % M


def gaps(x: int, q: int) -> list[int]:
    ones = [k for k, b in enumerate(format(x, f"0{q}b")) if b == "1"]
    return sorted((ones[(i + 1) % len(ones)] - ones[i]) % q for i in range(len(ones)))


@pytest.mark.parametrize("p,q", FRACTIONS)
def test_step_i_the_height_set_is_a_rotation_invariant(p, q):
    x, _ = ring2(p, q)
    g = heights(x, q)
    for t in range(q):
        assert sorted(heights(shift(x, t, q), q)) == sorted(v - g[t] for v in g)


@pytest.mark.parametrize("p,q", FRACTIONS)
def test_step_ii_mechanical_words_have_an_interval_of_heights(p, q):
    for r in range(q):
        assert set(heights(mechanical(p, q, r), q)) == set(range(r - q + 1, r + 1))


@pytest.mark.parametrize("p,q", [(p, q) for p, q in FRACTIONS if 3 <= p <= q - 3])
def test_step_iii_the_two_ring_two_height_sets(p, q):
    x, y = ring2(p, q)
    assert set(heights(x, q)) == {p - q - 2, p - q - 1, p} | set(range(p - q + 1, p - 2))
    assert set(heights(y, q)) == {p - q, p + 1, p + 2} | set(range(p - q + 3, p))
    assert normal(x, q) != normal(y, q)


@pytest.mark.parametrize("q", range(5, 61))
def test_step_iv_p_equal_one_differs_in_weight(q):
    x, y = ring2(1, q)
    assert (x, y) == (2 ** (q - 2) + 1, 7)


@pytest.mark.parametrize("q", [q for q in range(5, 61, 2)])
def test_step_iv_p_equal_two_differs_in_gaps(q):
    x, y = ring2(2, q)
    assert gaps(x, q) == sorted([1, (q - 1) // 2, (q - 1) // 2])
    assert gaps(y, q) == sorted([1, (q - 3) // 2, (q + 1) // 2])


@pytest.mark.parametrize("p,q", FRACTIONS)
def test_step_iv_conjugation_swaps_the_ring_two_angles(p, q):
    M = 2**q - 1
    x, y = ring2(p, q)
    xc, yc = ring2(q - p, q)
    assert (xc, yc) == ((M - y) % M, (M - x) % M)


@pytest.mark.parametrize("p,q", FRACTIONS)
def test_p8_distinct_orbits_of_period_q(p, q):
    x, y = ring2(p, q)
    orbit = {shift(x, d, q) for d in range(q)}
    assert y not in orbit and len(orbit) == q and len({shift(y, d, q) for d in range(q)}) == q
