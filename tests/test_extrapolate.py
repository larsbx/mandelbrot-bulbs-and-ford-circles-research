import mpmath as mp
import pytest
from bulbford.extrapolate import power_fit, linear_limit


def test_power_fit_exact_polynomial():
    with mp.workdps(90):                                    # inputs exact to working precision
        c = [mp.mpc(1, -2), mp.mpf(3), mp.mpc(0, 5)]
        rows = [(q, c[0] + c[1] / q + c[2] / q ** 2) for q in (10, 20, 40, 80)]
        assert all(abs(x - y) < 1e-60 for x, y in zip(power_fit(rows, 2), c))


def _differences(step, nu, d, n):
    for _ in range(n):
        yield d
        d = step(nu, d)


@pytest.mark.parametrize("step", [lambda nu, d: nu * mp.conj(d), lambda nu, d: nu * d], ids=["antilinear", "linear"])
def test_linear_limit_exact_for_real_linear_contractions(step):
    """Antilinear (Δ ↦ νΔ̄: the conjugated hierarchy) and holomorphic (Δ ↦ νΔ) contractions are both real-linear
    on ℂ ≅ ℝ², so the limit is recovered exactly from three differences."""
    with mp.workdps(60):
        nu, x0, d0 = mp.mpc(-0.0155, 0.0389), mp.mpc(0.3, -0.7), mp.mpc(0.02, 0.01)
        ds = list(_differences(step, nu, d0, 300))
        xs = [x0 + mp.fsum(ds[:k]) for k in range(6)]
        got, ev = linear_limit(xs)
        assert abs(got - (x0 + mp.fsum(ds))) < 1e-45
        assert [abs(e) for e in ev] == pytest.approx([abs(nu)] * 2, rel=1e-30)
