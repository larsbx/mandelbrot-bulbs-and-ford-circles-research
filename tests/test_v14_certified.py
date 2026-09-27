"""V14 with both proxies certified: |u_a|/2 on the P10 antipode boxes (Arb) against G_ant (P10)."""

import json
import runpy
from fractions import Fraction as F
from pathlib import Path

import pytest

from bulbford.antipode import u_half_abs
from bulbford.certify import Box, I

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = runpy.run_path(str(ROOT / "scripts" / "v14_certified.py"))
DOC = json.loads((ROOT / "data" / "v14_certified.json").read_text())


def test_table_is_current_and_recomputed_from_the_boxes():
    assert (ROOT / "data" / "v14_certified.json").read_text() == SCRIPT["render"]()


def test_fft_series_value_lies_within_1e_11_of_every_certified_u_bracket():
    """The V14 |u_a|/2 was a root of the FFT-truncated Taylor series; it is the same root to 1e-11."""
    for q in (59, 127, 251):
        fft = {r["p"]: F(r["G_ant_series"]) for r in json.loads((ROOT / "data" / f"taylor_q{q}.json").read_text())}
        for r in (r for r in DOC["rows"] if r["q"] == q):
            lo, hi = (F(x) for x in r["u_a_half"])
            assert lo - F(1, 10**11) <= fft[r["p"]] <= hi + F(1, 10**11), (r["p"], q)
            assert hi - lo < F(1, 10**12)


def test_certified_v14_column():
    """q²·max_p Δ is certified in intervals that round to the corrected register values 1.99, 2.06, 3.47,
    and Δ is strictly positive everywhere: the two proxies differ on every bulb."""
    expected = {59: (1.99, 8), 127: (2.06, 6), 251: (3.47, 1)}
    for s in DOC["summary"]:
        lo, hi = (float(F(x)) for x in s["q2_max_delta"])
        assert round(lo, 2) == round(hi, 2) == expected[s["q"]][0]
        assert s["argmax_p"] == expected[s["q"]][1]
    assert all(F(r["delta"][0]) > 0 for r in DOC["rows"])


def test_branch_that_cannot_be_separated_is_refused():
    """At c = 1/4 (λ = 1) the two branches of λ = 1 ∓ √(1 − 4c) meet; the enclosure refuses."""
    c = Box(I(F(1, 4) - F(1, 2**30), F(1, 4) + F(1, 2**30)), I(-F(1, 2**30), F(1, 2**30)))
    with pytest.raises(ValueError):
        u_half_abs(1, 3, c, 160)


def test_quoted_intervals_contain_the_certified_enclosures():
    """Every q²Δ interval quoted in the registers is rounded outward: it contains the stored enclosure."""
    import re

    stored = {s["q"]: tuple(F(x) for x in s["q2_max_delta"]) for s in DOC["summary"]}
    quoted = []
    for name in ("RESEARCH_bulb-ford-correction.md", "RESEARCH_bulb-ford-correction_finite.md"):
        text = (ROOT / name).read_text()
        quoted += [(F(a), F(b)) for a, b in re.findall(r"\[(\d\.\d{9,}), (\d\.\d{9,})\]", text)]
    assert len(quoted) >= 6
    for lo, hi in quoted:
        match = [q for q, (slo, shi) in stored.items() if lo <= slo and shi <= hi]
        assert match, (lo, hi)
    assert SCRIPT["outward"](F(-1, 3), 3, up=False) == "-0.334" and SCRIPT["outward"](F(-1, 3), 3, up=True) == "-0.333"
