"""Rotation cycles of doubling, characteristic arcs, and their Farey (Ford) order."""

import json
from fractions import Fraction as F
from pathlib import Path

import pytest

from bulbford.cf import coprime_numerators
from bulbford.wake import acts_as_rotation, double, farey, rotation_cycle, wake

PAIRS = [(p, q) for q in range(2, 17) for p in coprime_numerators(q)]


@pytest.mark.parametrize(
    "p,q,expected",
    [(1, 2, (F(1, 3), F(2, 3))), (1, 3, (F(1, 7), F(2, 7))), (2, 3, (F(5, 7), F(6, 7))),
     (1, 4, (F(1, 15), F(2, 15))), (2, 5, (F(9, 31), F(10, 31)))],
)
def test_known_wakes(p, q, expected):
    assert wake(p, q) == expected


@pytest.mark.parametrize("p,q", PAIRS)
def test_rotation_cycle_and_characteristic_arc(p, q):
    cycle = rotation_cycle(p, q)
    lo, hi = wake(p, q)
    assert len(set(cycle)) == q and all((2**q - 1) % x.denominator == 0 for x in cycle)
    assert acts_as_rotation(p, q, cycle)
    assert hi - lo == F(1, 2**q - 1)
    # the wrap-around arc through 0 is never the characteristic one
    assert 1 + cycle[0] - cycle[-1] > hi - lo


@pytest.mark.parametrize("p,q", PAIRS)
def test_conjugate_wake_is_the_mirror(p, q):
    lo, hi = wake(p, q)
    assert wake(q - p, q) == (1 - hi, 1 - lo)


@pytest.mark.parametrize("n", [5, 8, 12])
def test_wakes_are_disjoint_and_ordered_like_farey_fractions(n):
    wakes = [wake(x.numerator, x.denominator) for x in farey(n)]
    assert all(a[1] < b[0] for a, b in zip(wakes, wakes[1:]))


@pytest.mark.parametrize("n", [4, 7, 10])
def test_mediant_wake_sits_in_the_gap_between_farey_neighbours(n):
    for a, b in zip(farey(n), farey(n)[1:]):
        assert a.numerator * b.denominator - a.denominator * b.numerator == -1
        m = F(a.numerator + b.numerator, a.denominator + b.denominator)
        assert wake(a.numerator, a.denominator)[1] < wake(m.numerator, m.denominator)[0]
        assert wake(m.numerator, m.denominator)[1] < wake(b.numerator, b.denominator)[0]


def test_doubling_agrees_with_shared_r1_vectors():
    path = Path(__file__).resolve().parents[1] / "data" / "rational_dynamics_r1_vectors.json"
    for row in json.loads(path.read_text())["vectors"]:
        assert double(F(row["p"], row["q"])) == F(row["double_num"], row["double_den"])


def test_out_of_domain_rotation_number_is_refused():
    for p, q in [(0, 3), (3, 3), (2, 4), (-1, 5)]:
        with pytest.raises(ValueError):
            rotation_cycle(p, q)
