"""P12: certified continuation from the certified centre (P3) to the certified antipode (P10)."""

import json
import runpy
from fractions import Fraction as F
from pathlib import Path

import pytest

from bulbford.continuation import continue_centre_to_antipode

ROOT = Path(__file__).resolve().parents[1]
DOC = json.loads((ROOT / "experiments/data/continuation_certificates.json").read_text())
SCRIPT = runpy.run_path(str(ROOT / "experiments/scripts/certify_continuation.py"))
BOXES = {(b[0], b[1]): b for b in SCRIPT["boxes"]()}


def test_every_bulb_up_to_q16_is_accepted():
    assert len(DOC["certificates"]) == 79
    assert all(r["accepted"] for r in DOC["certificates"])


@pytest.mark.parametrize("row", [r for r in DOC["certificates"] if r["q"] <= 7], ids=lambda r: f'{r["p"]}/{r["q"]}')
def test_records_recompute_from_the_certificate_boxes(row):
    """`experiments/scripts/certify_continuation.py --check` recomputes all 79; CI recomputes q ≤ 7."""
    assert SCRIPT["record"](BOXES[(row["p"], row["q"])]) == row


def test_mismatched_antipode_is_not_accepted():
    """From the 1/3 centre to the 2/3 antipode the path leaves B_{1/3}: no continuation closes."""
    p, q, centre, _, _ = BOXES[(1, 3)]
    *_, ant_c, ant_z = BOXES[(2, 3)]
    assert not continue_centre_to_antipode(p, q, centre, ant_c, ant_z, min_ds=F(1, 2**14)).accepted

