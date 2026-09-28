"""V28–V30 (C20′, C21, C22): the plateaus, κ₀ and the 1/q constant of the two V25 families from their horn maps."""
import json
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from bulbford.implosion import (
    DPS, HALF, ONE, _expansion, fatou_coefficients, horn, kappa0, lavaurs_phase, phase_curvature, phi_in, psi_out,
)

KAPPA0 = mp.mpc("0.02382588740", "-0.05230465911")
KAPPA0_HALF = mp.mpc("0.01840526162", "0.04290372388")


@pytest.fixture(autouse=True)
def precision():
    with mp.workdps(DPS):
        yield


def test_the_first_fatou_coefficients():
    one, half = fatou_coefficients(ONE), fatou_coefficients(HALF)
    assert (one.principal, one.log, one.taylor[:4]) == ((F(-1),), F(1), (F(-1, 2), F(1, 3), F(-13, 36), F(113, 240)))
    assert (half.principal, half.log, half.taylor[:3]) == ((F(1, 4), F(1, 4)), F(11, 16), (F(-5, 16), F(75, 64), F(-149, 192)))


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
    """|κ(1/q) − κ₀ − 1/(2πiq)| = 2.9e-6 at q = 1009 (V28b); the 1/q term is the phase curvature (V30)."""
    rows = json.loads((Path(__file__).resolve().parents[1] / "experiments/data/kappa_q1009.json").read_text())
    k = mp.mpc(*next(r["kappa"] for r in rows if r["p"] == 1))
    assert abs(k - kappa0() - phase_curvature(1009, ONE) / 1009) < 4e-6
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


def test_the_half_germ_coordinates_conjugate_f_to_a_half_step():
    w, Z = mp.mpc("0.3", "0.1"), mp.mpc("0.2", "-3")
    assert abs(phi_in(HALF.f(w), HALF)[0] - phi_in(w, HALF)[0] - mp.mpf(1) / 2) < 1e-30
    a, b = psi_out(Z, HALF)[0], psi_out(Z + 1, HALF)[0]
    assert abs(b - HALF.f(HALF.f(a))) < 1e-30


def test_the_two_repelling_petals_carry_different_coordinates():
    """Ψ on R₋ is f ∘ Ψ(· − ½) on R₊, and it is not Ψ on R₊: the horn map has period 1."""
    Z = mp.mpc("0.2", "-3")
    assert abs(psi_out(Z, HALF, 1)[0] - HALF.f(psi_out(Z - mp.mpf(1) / 2, HALF)[0])) < 1e-30
    assert abs(psi_out(Z, HALF, 1)[0] - psi_out(Z, HALF)[0]) > 0.1


def test_the_half_horn_map_ends_are_the_log_branches():
    end = 11j * mp.pi / 16
    assert abs(horn(mp.mpc(0.3, 6), HALF)[0] + end) < 1e-11
    assert abs(horn(mp.mpc(0.3, -6), HALF)[0] - end) < 1e-11


def test_kappa0_half_does_not_depend_on_the_sampling_line():
    assert abs(kappa0(-2.4, HALF) - kappa0(-2.9, HALF)) < 1e-25
    assert abs(kappa0(germ=HALF) - KAPPA0_HALF) < 1e-11


def test_kappa0_half_meets_the_exact_index_at_q_1009():
    """q(κ − κ₀) = 0.356 + 0.241i at q = 1009 (V29b)."""
    rows = json.loads((Path(__file__).resolve().parents[1] / "experiments/data/kappa_q1009.json").read_text())
    k = mp.mpc(*next(r["kappa"] for r in rows if r["p"] == 504))
    assert abs(1009 * (k - kappa0(germ=HALF)) - mp.mpc("0.3558", "0.2414")) < 1e-3


def test_the_half_flank_plateau_is_a_lavaurs_multiplier():
    """|1 − L'| = 56.2344028 at σ∞ = −11πi/16, and F'_q meets L' at rate 1/q."""
    from bulbford.cycles import cycle_points
    from bulbford.implosion import fixed_point_of_cycle
    from multiplier_growth import named_angles

    gaps = []
    for q in (513, 1025):
        z = cycle_points((q - 1) // 2, q, named_angles((q - 1) // 2, q)["flank"])[1]
        slope = fixed_point_of_cycle(z, seeds=2, germ=HALF)[0][1]
        assert abs(abs(1 - slope) - mp.mpf("56.2344028")) < 1e-6
        gaps.append(abs(slope - complex(np.prod(2 * z))))
    assert 1.6 < gaps[0] / gaps[1] < 2.4


def test_the_phase_curvature_in_closed_form():
    """1/(2πi) for p = 1 (to O(q⁻⁴)) and i/π for p = (q − 1)/2 at every q."""
    assert abs(phase_curvature(1009, ONE) - 1 / (2j * mp.pi)) < 1e-11
    assert abs(phase_curvature(101, HALF) - 1j / mp.pi) < 1e-30


def test_the_phase_is_chosen_by_germ_value():
    from dataclasses import replace

    assert abs(phase_curvature(1009, replace(ONE)) - phase_curvature(1009, ONE)) < 1e-30
    with pytest.raises(ValueError):
        phase_curvature(1009, replace(ONE, name="other", end=-1))


def test_the_half_constant_is_the_curvature_plus_a_remainder():
    """V30: C = i/π + R, R = 0.354856 − 0.069314i, fitted on the FFT κ for q = 1025 … 8193."""
    from bulbford.taylor import kappa_fft

    data = json.loads((Path(__file__).resolve().parents[1] / "experiments/data/kappa_constant.json").read_text())
    C = complex(*data["fits"]["3 terms"]["C"])
    assert abs(C - complex(phase_curvature(4097, HALF)) - complex(0.354856, -0.069314)) < 2e-6
    assert data["fits"]["3 terms"]["max_residual"] < 2e-6
    row = next(r for r in data["rows"] if r["q"] == 1025)
    fresh = 1025 * (complex(kappa_fft(512, 1025).coeffs[2]) - complex(kappa0(germ=HALF)))
    assert abs(fresh - complex(*row["q_times_gap"])) < 1e-8
