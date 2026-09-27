"""V28 (C20′, C21): the p = 1 plateaus and κ₀ from the horn map of w + w²."""
import json
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from bulbford.implosion import (
    DPS, _expansion, fatou_coefficients, horn, kappa0, lavaurs_phase, phi_in, psi_out,
)

KAPPA0 = mp.mpc("0.02382588740", "-0.05230465911")


@pytest.fixture(autouse=True)
def precision():
    with mp.workdps(DPS):
        yield


def test_the_first_fatou_coefficients():
    assert fatou_coefficients()[:4] == (F(-1, 2), F(1, 3), F(-13, 36), F(113, 240))


@pytest.mark.parametrize("w", [mp.mpc("-2e-3", "5e-4"), mp.mpc("3e-3", "-1e-3")])
def test_the_expansion_solves_abel_to_its_order(w):
    outgoing = w.real > 0
    shift = _expansion(w + w * w, outgoing)[0] - _expansion(w, outgoing)[0] - 1
    assert abs(shift) < abs(w) ** 13


def test_the_coordinates_conjugate_g_to_a_translation():
    w, Z = mp.mpc("-0.3", "0.2"), mp.mpc("0.3", "1.1")
    assert abs(phi_in(w + w * w)[0] - phi_in(w)[0] - 1) < 1e-30
    a, b = psi_out(Z)[0], psi_out(Z + 1)[0]
    assert abs(b - a - a * a) < 1e-30


def test_the_horn_map_ends():
    assert abs(horn(mp.mpc(0.3, 8))[0] + 1j * mp.pi) < 1e-12
    assert abs(horn(mp.mpc(0.3, -8))[0] - 1j * mp.pi) < 1e-12


def test_the_lavaurs_phase_of_the_1_over_q_root():
    q = 1009
    eps = (mp.expjpi(mp.mpf(2) / q) - 1) / 2j
    assert abs(lavaurs_phase(q) - (q - mp.pi / eps)) < 1e-25
    assert abs(lavaurs_phase(q) - 1j * mp.pi) < 4 / q


def test_kappa0_does_not_depend_on_the_sampling_line():
    assert abs(kappa0(2.8) - kappa0(3.2)) < 1e-25
    assert abs(kappa0() - KAPPA0) < 1e-11


def test_kappa0_meets_the_exact_index_at_q_1009():
    """|κ(1/q) − κ₀ − 1/(2πiq)| = 2.9e-6 at q = 1009 (V28b)."""
    rows = json.loads((Path(__file__).resolve().parents[1] / "experiments/data/kappa_q1009.json").read_text())
    k = mp.mpc(*next(r["kappa"] for r in rows if r["p"] == 1))
    assert abs(k - kappa0() + 1j / (2 * mp.pi * 1009)) < 4e-6
    assert abs(k - kappa0()) > 1e-4


def test_the_flank_plateau_is_a_lavaurs_multiplier():
    """V27's flank intercept is |1 − L'| at the fixed point of L_{iπ}, and F'_q meets L' at rate q⁻²."""
    from bulbford.cycles import cycle_points
    from bulbford.implosion import fixed_point_of_cycle
    from multiplier_growth import named_angles

    gaps = []
    for q in (64, 128):
        z = cycle_points(1, q, named_angles(1, q)["flank"])[1]
        found = fixed_point_of_cycle(z, seeds=2)
        slopes = [s for _, s in found]
        assert len(slopes) == 2 and abs(slopes[0] - slopes[1]) < 1e-20
        assert abs(abs(1 - slopes[0]) - mp.mpf("51.1881731")) < 1e-6
        gaps.append(abs(slopes[0] - complex(np.prod(2 * z))))
    assert 3.5 < gaps[0] / gaps[1] < 4.5


def test_psi_out_needs_no_iterates_deep_in_the_petal():
    Z = mp.mpc(-2000, 0)
    a, b = psi_out(Z)[0], psi_out(Z + 1)[0]
    assert abs(b - a - a * a) < 1e-30
    assert abs(_expansion(a, outgoing=True)[0] - Z) < 1e-30
