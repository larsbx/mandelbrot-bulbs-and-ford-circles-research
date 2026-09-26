"""Joint box certificates: calculus replay, bulb-centre records, and fail-closed behaviour."""

import json
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from pathlib import Path

import pytest

from bulbford.certify import (
    Box,
    I,
    JointCertificate,
    Verdict,
    box_from_record,
    certify_center,
    certify_type,
    failed_exclusions,
    forbidden_pairs,
    intended_pairs,
    replay,
)

DATA = Path(__file__).resolve().parents[1] / "data" / "center_certificates.json"
RECORDS = json.loads(DATA.read_text())["certificates"]


# --- replay of finite-mandelbrot-research reference semantics ----------------------------------
# reference/python/interval/interval_exclusion_reference.py pins these two boxes and counts.


def test_reference_c_minus_2_box():
    beta = Box(I(F(-33, 16), F(-31, 16)), I(F(-1, 16), F(1, 16)))
    assert len(forbidden_pairs(2, 1, 3)) == 5
    assert failed_exclusions(beta, 2, 1, 3) == ()


def test_reference_m41_box():
    h, den = F(1, 2**25), 2**49
    re, im = F(-56912193317957, den), F(538341446717435, den)
    beta = Box(I(re - h, re + h), I(im - h, im + h))
    assert len(forbidden_pairs(4, 1, 6)) == 18
    assert failed_exclusions(beta, 4, 1, 6) == ()


def test_c_minus_2_joint_certificate_type_2_1():
    beta = Box(I(F(-2) - F(1, 64), F(-2) + F(1, 64)), I(F(-1, 64), F(1, 64)))
    assert certify_type(beta, 2, 1, 3).verdict is Verdict.ACCEPTED


def test_intended_pairs_at_minimal_and_larger_horizon():
    assert intended_pairs(0, 3, 3) == {(0, 3)}
    assert intended_pairs(0, 3, 6) == {(0, 3), (1, 4), (2, 5), (3, 6), (0, 6)}
    assert forbidden_pairs(0, 3, 3) == set(combinations(range(4), 2)) - {(0, 3)}


# --- stored bulb-centre certificates --------------------------------------------------------------


def test_every_stored_centre_is_accepted_on_replay():
    for row in RECORDS:
        assert row["verdict"] == "ACCEPTED"
        assert replay(row).verdict is Verdict.ACCEPTED, (row["p"], row["q"])


def test_one_certified_centre_per_satellite_and_boxes_disjoint():
    by_q: dict[int, list[Box]] = {}
    for row in RECORDS:
        by_q.setdefault(row["q"], []).append(box_from_record(row))
    for q, boxes in by_q.items():
        assert len(boxes) == sum(1 for p in range(1, q) if gcd(p, q) == 1)
        for a, b in combinations(boxes, 2):
            assert a.re.hi < b.re.lo or b.re.hi < a.re.lo or a.im.hi < b.im.lo or b.im.hi < a.im.lo


def test_stored_records_are_current():
    import runpy

    script = runpy.run_path(str(Path(__file__).resolve().parents[1] / "scripts" / "certify_centers.py"))
    assert DATA.read_text() == script["render"]()


# --- fail closed --------------------------------------------------------------------------------


def _row(p, q):
    return next(r for r in RECORDS if (r["p"], r["q"]) == (p, q))


def test_wrong_period_claim_is_not_accepted():
    """The period-3 centre is a root of Q_6 too; type (0, 6) must fail the (0, 3) exclusion."""
    beta = box_from_record(_row(1, 3))
    cert = certify_type(beta, 0, 6, 6)
    assert cert.krawczyk_inclusion
    assert (0, 3) in cert.failed_exclusions
    assert cert.verdict is Verdict.INCONCLUSIVE


def test_oversized_box_is_inconclusive_not_rejected_claim():
    cert = certify_center(1, 3, complex(-0.1225611668766536, 0.7448617666197442), radius_bits=3)
    assert cert.verdict is Verdict.INCONCLUSIVE


def test_box_without_the_root_fails_krawczyk():
    beta = Box.around(F(1, 2), F(1, 2), F(1, 2**40))
    assert not certify_type(beta, 0, 5).krawczyk_inclusion


@pytest.mark.parametrize("inclusion,failures", [(True, ((1, 2),)), (False, ()), (False, ((0, 1),))])
def test_joint_verdict_needs_both_halves_on_the_same_box(inclusion, failures):
    beta = Box.point(0)
    assert JointCertificate(0, 1, 1, beta, inclusion, failures, 64).verdict is Verdict.INCONCLUSIVE


def test_malformed_type_is_refused():
    with pytest.raises(ValueError):
        certify_type(Box.point(0), 0, 0)
    with pytest.raises(ValueError):
        certify_type(Box.point(0), 1, 2, horizon=2)
    with pytest.raises(ValueError):
        I(F(1), F(0))
