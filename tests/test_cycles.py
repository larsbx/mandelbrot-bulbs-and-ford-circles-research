"""P6 of the finite register: ι_{p/q} = −Σ_{z ≠ z₀} 1/(1 − ρ_z) over the fixed points of f^q."""

from fractions import Fraction as F

import numpy as np
import pytest

from bulbford.cf import coprime_numerators
from bulbford.cycles import (
    _orbit, cycle_terms, cycle_through, fixed_points, flank_angles, index_by_cycles, rotation_number,
)
from bulbford.index import index_complex

PAIRS = [(p, q) for q in range(2, 13) for p in coprime_numerators(q)]


def test_q2_is_the_beta_point_alone():
    terms = cycle_terms(1, 2)
    assert len(terms) == 1 and terms[0].period == 1 and terms[0].angles == (0,)
    assert abs(terms[0].multiplier - 3.0) < 1e-12  # f'(β) = 2β, β = 3/2
    assert abs(index_by_cycles(1, 2) - 0.125) < 1e-12


@pytest.mark.parametrize("p,q", PAIRS)
def test_the_cycle_sum_is_the_parabolic_index(p, q):
    assert abs(index_by_cycles(p, q) - index_complex(p, q)) < 1e-12


@pytest.mark.parametrize("p,q", [(1, 7), (3, 8), (4, 9), (7, 12)])
def test_every_other_fixed_point_once_on_its_own_ray(p, q):
    fp = fixed_points(p, q)
    M = 2**q - 1
    others = np.setdiff1d(np.arange(M), fp.alpha_angles)
    z = fp.z[others]
    assert len(others) == 2**q - q - 1
    # Distinct points, each a fixed point of f^q, and f carries the ray of j to the ray of 2j.
    assert np.min(np.abs(z[:, None] - z[None, :]) + np.eye(len(z))) > 1e-6
    assert np.max(np.abs(fp.f(fp.z[others]) - fp.z[(2 * others) % M])) < 1e-10


def test_conjugate_numerators_give_conjugate_sums():
    for q in (5, 7, 9):
        for p in coprime_numerators(q):
            assert abs(index_by_cycles(q - p, q) - np.conj(index_by_cycles(p, q))) < 1e-10


def test_rotation_cycles_of_the_same_denominator_are_all_present():
    q, p = 7, 2
    found = {t.rotation for t in cycle_terms(p, q) if t.rotation is not None and t.period == q}
    assert found == {F(k, q) for k in range(1, q) if k != p}


def test_rotation_number_of_an_angle_orbit():
    # Angles 1/7, 2/7, 4/7: doubling shifts the sorted list by one place.
    assert rotation_number((1, 2, 4), 7) == F(1, 3)
    # 3/15 -> 6/15 -> 12/15 -> 9/15: sorted 3, 6, 9, 12; doubling sends positions 0,1,3,2 -> 1,3,2,0.
    assert rotation_number((3, 6, 12, 9), 15) is None


@pytest.mark.parametrize("q", range(3, 21))
def test_the_two_flank_angles_lie_on_one_doubling_orbit(q):
    """Exact, in integers: for every p, the orbit of α_{w−1} + 1 passes through α_{w+2} − 1."""
    M = 2**q - 1
    for p in coprime_numerators(q):
        a, b = flank_angles(p, q)
        orbit = _orbit(a, M)
        assert b in orbit and len(orbit) == q


@pytest.mark.parametrize("p,q", [(p, q) for p, q in PAIRS if q >= 3])
def test_the_flank_cycle_carries_the_largest_term(p, q):
    top = max(cycle_terms(p, q), key=lambda t: abs(t.contribution))
    assert set(flank_angles(p, q)) <= set(top.angles)


@pytest.mark.parametrize("p,q", [(2, 7), (4, 9), (3, 11), (5, 12)])
def test_a_single_orbit_gives_the_same_term_as_the_full_sum(p, q):
    for t in cycle_terms(p, q):
        assert abs(cycle_through(p, q, t.angles[-1]).contribution - t.contribution) < 1e-12


@pytest.mark.parametrize("q", [40, 64])
def test_single_orbits_verify_beyond_the_full_sum(q):
    assert cycle_through(1, q, 2 ** (q - 2) + 1).period == q  # ring-2 X₂ for p = 1


def test_the_two_crossings_of_v26():
    from scripts.crossings import record

    assert record(1, 7)["ratio"] < 1 < record(1, 8)["ratio"]
    assert record(5, 11)["ratio"] < 1 < record(6, 13)["ratio"]


@pytest.mark.parametrize("p,q", [(1, 5), (2, 7), (3, 8), (1, 13), (6, 13), (1, 40), (5, 64)])
def test_p11_the_beta_term_in_closed_form(p, q):
    """β = 1 − λ₀/2 has multiplier 2 − λ₀, so its term is −1/(1 − (2 − λ₀)^q)."""
    lam = np.exp(2j * np.pi * p / q)
    assert abs(cycle_through(p, q, 0).contribution - (-1 / (1 - (2 - lam) ** q))) < 1e-9
    assert abs(abs(2 - lam) ** 2 - (1 + 8 * np.sin(np.pi * p / q) ** 2)) < 1e-12


def test_beta_takes_the_lead_among_p_equal_one_intruders_at_q_22():
    from scripts.crossings import record

    beta = lambda q: abs(1 / (1 - (2 - np.exp(2j * np.pi / q)) ** q))  # noqa: E731
    assert record(1, 21)["intruders"][0] > beta(21)
    assert abs(record(1, 22)["intruders"][0] - beta(22)) < 1e-9
