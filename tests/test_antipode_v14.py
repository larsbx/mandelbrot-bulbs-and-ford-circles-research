"""Certified antipodes at the V14 denominators q = 59, 127, 251 (every unit p ≤ q/2).

CI replays all of q = 59 and a sample at q = 127, 251 that includes the records needing the smaller
box of the precision ladder; `python experiments/scripts/certify_bulbs.py --check` replays everything.
"""

import json
from fractions import Fraction as F
from pathlib import Path

import pytest

from bulbford.antipode import replay
from bulbford.certify import Verdict

ROOT = Path(__file__).resolve().parents[1]
RECORDS = json.loads((ROOT / "experiments/data/antipode_certificates_v14.json").read_text())["certificates"]
QS = tuple(sorted({r["q"] for r in RECORDS}))  # the V14 denominators certified so far


def taylor(q: int) -> dict[int, dict]:
    return {r["p"]: r for r in json.loads((ROOT / "experiments/data" / f"taylor_q{q}.json").read_text())}


def bracket(row) -> tuple[F, F]:
    lo, hi = row["G_ant_bounds"]
    return F(lo), F(hi)


def _sample():
    picked = []
    for q in QS[1:]:
        rows = [r for r in RECORDS if r["q"] == q]
        picked += rows[:1] + [r for r in rows if r["prec"] > 160][:1]
    return picked


def test_coverage_matches_the_v14_sweeps():
    for q in QS:
        assert sorted(r["p"] for r in RECORDS if r["q"] == q) == sorted(taylor(q))
    assert all(r["verdict"] == "ACCEPTED" for r in RECORDS)


@pytest.mark.parametrize("row", [r for r in RECORDS if r["q"] == 59] + _sample(), ids=lambda r: f'{r["p"]}/{r["q"]}')
def test_replay_from_boxes_alone(row):
    cert = replay(row)
    assert cert.verdict is Verdict.ACCEPTED
    g = cert.g()
    assert [str(g.lo), str(g.hi)] == row["G_ant_bounds"]


def test_floating_point_instrument_lies_within_1e_11_of_every_bracket():
    for q in QS:
        tay = taylor(q)
        for row in (r for r in RECORDS if r["q"] == q):
            lo, hi = bracket(row)
            g = F(tay[row["p"]]["G_ant"])
            assert lo - F(1, 10**11) <= g <= hi + F(1, 10**11), (row["p"], q)
            assert hi - lo < F(1, 10**12)


def test_v14_column_from_certified_g_ant():
    """Δ = max_p |(|u_a|/2) − G_ant| with G_ant certified (|u_a|/2 is still the FFT instrument).
    q²Δ rounds to 1.99, 2.06, 3.47; the previously tabulated 1.98, 2.10 came from rounding Δ first."""
    expected = {59: 1.99, 127: 2.06, 251: 3.47}
    for q in QS:
        tay = taylor(q)
        delta = max(
            max(abs(F(tay[r["p"]]["G_ant_series"]) - e) for e in bracket(r)) for r in RECORDS if r["q"] == q
        )
        assert round(float(delta) * q * q, 2) == expected[q]


def test_an_exploding_enclosure_is_inconclusive_not_slow():
    """A box far too wide for a long orbit trips the 2^40 guard: INCONCLUSIVE, with the BLOWN marker."""
    import time

    from bulbford.antipode import BLOWN, check_antipode
    from bulbford.certify import I, box_from_numerators

    row = next(r for r in RECORDS if r["q"] == max(QS))
    grow = lambda b: type(b)(I(b.re.lo - F(1, 2**12), b.re.hi + F(1, 2**12)), I(b.im.lo - F(1, 2**12), b.im.hi + F(1, 2**12)))
    z, c = (grow(box_from_numerators(row[k], row["prec"])) for k in ("z_box_numerators", "c_box_numerators"))
    t0 = time.perf_counter()
    cert = check_antipode(row["p"], row["q"], z, c, row["prec"])
    assert cert.verdict is Verdict.INCONCLUSIVE and cert.failed_exclusions == BLOWN
    assert time.perf_counter() - t0 < 30
