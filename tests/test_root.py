"""P14: the certified centre's component has the satellite root c_root(p/q) on its boundary."""

import json
from fractions import Fraction as F
from pathlib import Path

import pytest
from flint import acb, arb, ctx

from bulbford.antipode import lambda_box, root_box
from bulbford.continuation import _acb, _arb_exact
from bulbford.index import lambda_ball
from bulbford.root import _circle, _count, _s_and_ds, _units, certify_root_side

ROOT = Path(__file__).resolve().parents[1]
DOC = json.loads((ROOT / "experiments/data/root_certificates.json").read_text())

import certify_root  # noqa: E402  (experiments/scripts is on the pytest path)

BY_KEY = {(r["p"], r["q"]): r for r in DOC["certificates"]}
CENTRES = {(p, q): centre for p, q, centre in certify_root.jobs()}


@pytest.fixture(autouse=True)
def work_precision():
    old = ctx.prec
    try:
        ctx.prec = 128
        yield
    finally:
        ctx.prec = old


def _ball(row):
    m = acb(_arb_exact(F(row["disk_centre"][0])), _arb_exact(F(row["disk_centre"][1])))
    c1 = acb(_arb_exact(F(row["c1"][0])), _arb_exact(F(row["c1"][1])))
    return m, c1, _arb_exact(F(row["r"])), _arb_exact(F(row["r_inner"]))


def test_every_stored_bulb_is_accepted():
    assert all(r["accepted"] for r in DOC["certificates"])
    assert {(r["p"], r["q"]) for r in DOC["certificates"]} == {(1, 2), (1, 3), (2, 3), (1, 4), (3, 4)}


@pytest.mark.parametrize("key", sorted(BY_KEY))
def test_every_stored_record_replays_from_exact_endpoints(key):
    row = BY_KEY[key]
    assert certify_root.replay(row, CENTRES[key]) == row


def test_record_for_one_half_recomputes():
    """Also check the untrusted survey regenerates the same exact input for 1/2."""
    job = next(j for j in certify_root.jobs(2) if j[:2] == (1, 2))
    assert certify_root.record(job) == BY_KEY[(1, 2)]


def test_argument_principle_counts_q_plus_one_at_the_root():
    for p, q in [(1, 3), (1, 4)]:
        m, _, r, r_in = _ball(BY_KEY[(p, q)])
        c_root = _acb(root_box(lambda_box(p, q)))
        for radius in (r, r_in):
            n = _count(q, c_root, m, radius)
            assert abs(n - (q + 1)).upper() < 1e-10


@pytest.mark.parametrize("pole, expected", [(acb(0), 1), (acb("0.25", "0.25"), 1), (acb(2), 0)])
def test_algebraic_contour_mean_counts_a_known_pole(pole, expected):
    result = _circle(arb(1), lambda u: u / (u - pole))
    assert result.contains(acb(expected))
    assert result.rad() < arb(2) ** -35


def test_contour_with_a_pole_on_it_is_refused():
    assert not _circle(arb(1), lambda u: u / (u - 1)).is_finite()


def test_contour_error_bound_covers_a_nonzero_sampling_tail():
    # At this precision the finite root mean's bias exceeds its roundoff.
    # The exact residue is 1, so it must be covered by the added Laurent error.
    ctx.prec = 256
    pole = acb(arb(55) / 64)
    result = _circle(arb(1), lambda u: u / (u - pole))
    assert result.contains(acb(1))
    assert result.rad() < arb(2) ** -35


def test_contour_crossing_the_log_cut_is_refused():
    assert not _circle(arb(1), lambda u: u.log(analytic=True)).is_finite()


def test_log_moment_and_its_parameter_derivative_at_period_doubling():
    lam = lambda_ball(1, 2, 128)
    s, ds = _s_and_ds(2, acb("-0.75"), acb("-0.5"), arb("0.3"), lam)
    assert s.contains(acb(0)) and s.rad() < arb(2) ** -35
    assert ds.contains(acb(3)) and ds.rad() < arb(2) ** -35


def test_contour_nodes_use_the_requested_precision_and_restore_context():
    ctx.prec = 53
    _units.cache_clear()
    units = _units(128, 400)
    assert ctx.prec == 53
    ctx.prec = 400
    assert all((u ** 128).contains(acb(1)) and u.rad() < arb(2) ** -350 for u in units)


def test_disk_catching_other_periodic_points_is_refused():
    p, q = 1, 3
    m, c1, _, r_in = _ball(BY_KEY[(p, q)])
    res = certify_root_side(p, q, _acb(root_box(lambda_box(p, q))), c1, m, arb(3), r_in, min_dt=F(1, 2**6))
    assert not res.accepted


def test_log_not_holomorphic_on_the_inner_disk_is_refused():
    """2r′ ≥ 1 lets D̄′ reach z = 0, where log(2z/λ₀) is singular: refused before any integration."""
    p, q = 1, 2
    m, c1, r, _ = _ball(BY_KEY[(p, q)])
    res = certify_root_side(p, q, _acb(root_box(lambda_box(p, q))), c1, m, r, arb("0.6"))
    assert not res.accepted and "holomorphic" in res.reason


def test_inner_disk_whose_image_leaves_the_outer_disk_is_refused():
    p, q = 1, 3
    m, c1, _, r_in = _ball(BY_KEY[(p, q)])
    res = certify_root_side(p, q, _acb(root_box(lambda_box(p, q))), c1, m, r_in, r_in, min_dt=F(1, 2**6))
    assert not res.accepted and "f(D′)" in res.reason
