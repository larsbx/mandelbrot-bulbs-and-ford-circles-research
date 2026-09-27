"""P7, the flank-orbit lemma, one test per step of its proof, in integers for q ≤ 60."""

from math import gcd

import pytest

from bulbford.cycles import flank_angles
from bulbford.wake import mechanical, rotation_cycle, wake

FRACTIONS = [(p, q) for q in range(3, 61) for p in range(1, q) if gcd(p, q) == 1]


def bits(x: int, q: int) -> list[int]:
    return [int(b) for b in format(x, f"0{q}b")]


def shift(x: int, s: int, q: int) -> int:
    """σ^s: the cyclic left shift of the q-bit word, i.e. doubling s times modulo 2^q − 1."""
    return (x * pow(2, s, 2**q - 1)) % (2**q - 1)


@pytest.mark.parametrize("p,q", FRACTIONS)
def test_step_i_doubling_moves_the_intercept(p, q):
    pbar = pow(p, -1, q)
    for r in range(q):
        assert shift(mechanical(p, q, r), 1, q) == mechanical(p, q, (r + p) % q)
        assert shift(mechanical(p, q, r), pbar, q) == mechanical(p, q, (r + 1) % q)


@pytest.mark.parametrize("p,q", FRACTIONS)
def test_step_ii_neighbours_differ_by_one_adjacent_pair(p, q):
    pbar = pow(p, -1, q)
    for r in range(q - 1):
        k1 = (q - p - 1 - r) * pbar % q
        a, b = bits(mechanical(p, q, r), q), bits(mechanical(p, q, r + 1), q)
        assert [k for k in range(q) if a[k] != b[k]] == sorted({k1, (k1 + 1) % q})
        assert (a[k1], b[k1]) == (0, 1)


@pytest.mark.parametrize("p,q", FRACTIONS)
def test_step_iii_order_and_the_characteristic_arc(p, q):
    A = [mechanical(p, q, r) for r in range(q)]
    assert all(x < y for x, y in zip(A, A[1:]))
    assert [r for r in range(q - 1) if A[r + 1] - A[r] == 1] == [p - 1]


@pytest.mark.parametrize("p,q", [(p, q) for p, q in FRACTIONS if q <= 16])
def test_the_words_are_the_rotation_cycle_of_wake(p, q):
    M = 2**q - 1
    A = [mechanical(p, q, r) for r in range(q)]
    assert A == sorted(int(x * M) for x in rotation_cycle(p, q))
    assert int(wake(p, q)[0] * M) == A[p - 1]
    assert flank_angles(p, q) == ((A[(p - 2) % q] + 1) % M, (A[(p + 1) % q] - 1) % M)


@pytest.mark.parametrize("p,q", FRACTIONS)
def test_p7_the_second_flank_angle_is_pbar_doublings_after_the_first(p, q):
    M = 2**q - 1
    x = (mechanical(p, q, (p - 2) % q) + 1) % M
    y = (mechanical(p, q, (p + 1) % q) - 1) % M
    assert y == shift(x, pow(p, -1, q), q)
    assert all(shift(x, d, q) != x for d in range(1, q) if q % d == 0)  # exact period q


@pytest.mark.parametrize("q", range(3, 61))
def test_step_v_p_equal_one_in_closed_form(q):
    M = 2**q - 1
    assert [mechanical(1, q, r) for r in range(q)] == [2**r for r in range(q)]
    assert (mechanical(1, q, q - 1) + 1, mechanical(1, q, 2) - 1) == (2 ** (q - 1) + 1, 3)
    assert shift(2 ** (q - 1) + 1, 1, q) == 3 % M
