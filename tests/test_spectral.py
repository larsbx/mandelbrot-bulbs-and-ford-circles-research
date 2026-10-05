"""Ramanujan sums and the truncated Brjuno sum used by the spectral analyses (V31, V33)."""
from fractions import Fraction
from math import log, sqrt

import numpy as np
import pytest
from sympy import mobius, totient

from jump_law import centered_r2
from spectral import brjuno, joint_brjuno_fit, ramanujan


@pytest.mark.parametrize("q", range(1, 31))
def test_ramanujan_sum_values(q):
    m = np.array([1, q, 2 * q])
    assert np.allclose(ramanujan(q, m), [int(mobius(q)), int(totient(q)), int(totient(q))])


def test_brjuno_golden_mean():
    """B(1/φ) = log φ · Σ φ^{-n} = log φ / (1 − 1/φ); Fibonacci ratios F_n/F_{n+1} approximate 1/φ."""
    fib = [1, 1]
    while len(fib) < 60:
        fib.append(fib[-1] + fib[-2])
    phi = (1 + sqrt(5)) / 2
    assert brjuno(Fraction(fib[-2], fib[-1]), T=10**9) == pytest.approx(log(phi) / (1 - 1 / phi), rel=1e-9)


def test_brjuno_truncation_and_period():
    x = Fraction(355, 1133)
    assert brjuno(x + 3, 50) == brjuno(x, 50)
    assert brjuno(x, 1) == pytest.approx(log(1 / float(x)))


def test_joint_fit_recovers_two_families():
    """y = Σ A c_{q'}(m)/m + Σ B c_{q'}(m)/m² is recovered exactly by the joint (log, V) fit."""
    from spectral import ramanujan_fit
    m = np.arange(1, 201)
    A, B = np.array([0.3, -0.2, 0.1, 0.05]), np.array([-0.1, 0.4, 0.0, 0.2])
    y = sum(A[b - 1] * ramanujan(b, m) / m + B[b - 1] * ramanujan(b, m) / m**2 for b in range(1, 5))
    (fa, fb), r2 = ramanujan_fit(y, 4, 200, 4, g=(lambda k: 1.0 / k, lambda k: 1.0 / k**2), background=lambda k: 1.0 / k**4)
    assert np.allclose(fa, A) and np.allclose(fb, B) and r2 == pytest.approx(1.0)


def test_joint_brjuno_fit_uses_shared_kappa_smooth_series(tmp_path, monkeypatch):
    """The shared κ Fourier component does not leak a stable Brjuno amplitude."""
    import spectral
    qs, amplitude = [101, 103, 107], -0.004
    def table(q):
        p = np.arange(1, q)
        x = np.array([pow(int(a), -1, q) for a in p]) / q
        b = np.array([brjuno(Fraction(min(a, q-a), q), 45) for a in np.rint(x*q).astype(int)])
        y = 0.03 + 0.01*np.cos(2*np.pi*x) - 0.002*np.cos(4*np.pi*x) + amplitude*b
        return p, x, y.astype(complex)
    monkeypatch.setattr(spectral, "kappa_table", table)
    fit = joint_brjuno_fit(qs, 2, 45, bounded=0)
    assert all(v == pytest.approx(amplitude, abs=2e-10) for v in fit["amplitude"].values())


def test_centered_r2_counts_residual_mean_bias():
    y = np.array([1.0, 2.0, 3.0])
    residual = np.array([1.0, 1.0, 1.0])
    assert centered_r2(y, residual) == pytest.approx(-0.5)
