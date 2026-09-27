"""Certified ζ_q, antipodes (ρ = −1) and G_ant enclosures: replay, cross-checks, fail-closed behaviour."""

import cmath
import json
import runpy
from fractions import Fraction as F
from math import pi
from pathlib import Path

import pytest

from bulbford.antipode import (
    check_antipode,
    divide_positive,
    lambda_box,
    power,
    quadrance,
    replay,
    root_box,
    sqrt_bounds,
    unity_boxes,
    zeta_box,
)
from bulbford.certify import ONE, Box, I, Verdict, box_from_numerators
from bulbford.dynamics import MAIN2, bulb

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "antipode_certificates.json"
RECORDS = json.loads(DATA.read_text())["certificates"]


def contains(box: Box, z: complex, slack: float = 0.0) -> bool:
    return box.re.lo - slack <= z.real <= box.re.hi + slack and box.im.lo - slack <= z.imag <= box.im.hi + slack


# --- ζ_q without trigonometry ------------------------------------------------------------------


@pytest.mark.parametrize("q", range(1, 17))
def test_zeta_box_is_a_primitive_root_with_smallest_positive_argument(q):
    z = zeta_box(q)
    assert contains(z, cmath.exp(2j * pi / q), 1e-15)                      # cross-check only
    assert (power(z, q) - ONE).re.contains_zero() and (power(z, q) - ONE).im.contains_zero()
    assert all((power(z, d) - ONE).excludes_zero() for d in range(1, q))   # exact order q


def test_unity_boxes_are_disjoint_and_complete():
    boxes = unity_boxes(7)
    assert len(boxes) == 7
    assert all(b.re.hi < o.re.lo or o.re.hi < b.re.lo or b.im.hi < o.im.lo or o.im.hi < b.im.lo
               for i, b in enumerate(boxes) for o in boxes[i + 1:])


def test_root_box_encloses_the_cardioid_root():
    for p, q in [(1, 3), (2, 5), (5, 12)]:
        lam = cmath.exp(2j * pi * p / q)
        assert contains(root_box(lambda_box(p, q)), lam / 2 - lam * lam / 4, 1e-15)


# --- stored antipode certificates -----------------------------------------------------------------


def test_every_stored_antipode_is_accepted_on_replay_and_bounds_are_reproduced():
    for row in RECORDS:
        cert = replay(row)
        assert cert.verdict is Verdict.ACCEPTED, (row["p"], row["q"])
        g = cert.g()
        assert [str(g.lo), str(g.hi)] == row["G_ant_bounds"]
        assert g.hi - g.lo < F(1, 10**15)


def test_bounds_agree_with_the_floating_point_instrument():
    for row in RECORDS[::7]:
        lo, hi = (float(F(x)) for x in row["G_ant_bounds"])
        assert lo - 1e-12 <= bulb(MAIN2, row["p"], row["q"]).G_ant <= hi + 1e-12


def test_period_two_bulb_is_exact():
    """P1: the period-2 disk has G_ant = 1 exactly; the enclosure must contain 1."""
    g = replay(next(r for r in RECORDS if r["q"] == 2)).g()
    assert g.lo <= 1 <= g.hi


def test_conjugate_bulbs_have_equal_bounds_up_to_width():
    rows = {(r["p"], r["q"]): r for r in RECORDS}
    for (p, q), r in rows.items():
        a, b = (I(*(F(x) for x in rows[k]["G_ant_bounds"])) for k in [(p, q), (q - p, q)])
        assert not (a.hi < b.lo or b.hi < a.lo)


def test_stored_records_are_current():
    script = runpy.run_path(str(ROOT / "scripts" / "certify_bulbs.py"))
    assert DATA.read_text() == script["OUTPUTS"][DATA]()


# --- fail closed --------------------------------------------------------------------------------


def _boxes(p, q):
    row = next(r for r in RECORDS if (r["p"], r["q"]) == (p, q))
    return box_from_numerators(row["z_box_numerators"], row["prec"]), box_from_numerators(row["c_box_numerators"], row["prec"])


def test_wrong_period_is_not_accepted():
    """The period-3 antipode cycle also solves f^6(z) = z; its multiplier there is +1, so K fails and the
    (0, 3) exclusion fails: two independent refusals."""
    z, c = _boxes(1, 3)
    cert = check_antipode(1, 6, z, c)
    assert cert.verdict is Verdict.INCONCLUSIVE
    assert not cert.krawczyk_inclusion
    assert (0, 3) in cert.failed_exclusions


def test_displaced_box_fails_krawczyk():
    z, c = _boxes(2, 5)
    shift = Box.point(F(1, 2**20))
    assert not check_antipode(2, 5, z, c + shift).krawczyk_inclusion


def test_oversized_box_is_inconclusive():
    z, c = _boxes(1, 4)
    grow = lambda b: Box(I(b.re.lo - F(1, 8), b.re.hi + F(1, 8)), I(b.im.lo - F(1, 8), b.im.hi + F(1, 8)))
    assert check_antipode(1, 4, grow(z), grow(c)).verdict is Verdict.INCONCLUSIVE


def test_interval_helpers_refuse_out_of_domain():
    with pytest.raises(ValueError):
        divide_positive(I.point(1), I(F(-1), F(1)))
    with pytest.raises(ValueError):
        sqrt_bounds(I(F(-1), F(1)))
    assert quadrance(Box(I(F(-1), F(2)), I.point(0))) == I(F(0), F(4))
    s = sqrt_bounds(I(F(2), F(2)), 40)
    assert s.lo**2 <= 2 <= s.hi**2 and s.hi - s.lo <= F(1, 2**39)


@pytest.mark.parametrize("q", [127, 251])
def test_long_powers_do_not_wrap(q):
    """Square-and-multiply keeps ζ_q^q tight; n successive box products would inflate a 2^-100 box past use."""
    z = zeta_box(q)
    w = power(z, q) - ONE
    assert w.re.contains_zero() and w.im.contains_zero()
    assert w.re.hi - w.re.lo < F(1, 2**80)
