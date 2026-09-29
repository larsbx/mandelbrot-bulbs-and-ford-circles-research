"""Spectral tools for κ(p̄/q): Fourier coefficients in x̃ = p̄/q and Ramanujan-sum fits (V31, next move 2).

A singularity at every reduced p'/q' with amplitude A_{q'} (depending on q' only) and profile with
Fourier decay g(m) contributes Σ_{q'} A_{q'} c_{q'}(m) g(m) to the m-th coefficient, c_{q'} the Ramanujan
sum. Profiles: a jump in a sine series (g = 1/(πm)), a log singularity (g = 1/m), a V-cusp (g = 1/m²),
|δ| log|δ| (g = log m / m²), and the Hölder cusp |δ|^{s−1} (g = m^{-s}).
"""
from __future__ import annotations

import json
from fractions import Fraction
from math import gcd, log

import numpy as np

from paths import DATA


def kappa_table(q: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(p, p̄/q, κ) over all units p mod q, from kappa_q<q>.json (p ≤ q/2) and κ(q − p) = conj κ(p)."""
    half = {r["p"]: complex(*r["kappa"]) for r in json.loads((DATA / f"kappa_q{q}.json").read_text())}
    full = {**half, **{q - p: k.conjugate() for p, k in half.items()}}
    p = np.array(sorted(full))
    return p, np.array([pow(int(a), -1, q) for a in p]) / q, np.array([full[a] for a in p])


def fourier(x: np.ndarray, y: np.ndarray, m: np.ndarray, trig) -> np.ndarray:
    return np.array([2 * np.mean(y * trig(2 * np.pi * k * x)) for k in m])


def coefficients(q: int, mmax: int, bounded: int = 4) -> dict:
    """Sine coefficients S_m of Im κ and cosine coefficients C_m of Re κ, bulbs with p ≤ bounded excluded."""
    p, x, k = kappa_table(q)
    keep = np.minimum(p, q - p) > bounded
    m = np.arange(1, mmax + 1)
    return {"m": m, "S": fourier(x[keep], k[keep].imag, m, np.sin),
            "C": fourier(x[keep], k[keep].real - k[keep].real.mean(), m, np.cos)}


def ramanujan(qq: int, m: np.ndarray) -> np.ndarray:
    return sum(np.cos(2 * np.pi * a * m / qq) for a in range(qq) if gcd(a, qq) == 1)


def ramanujan_fit(y: np.ndarray, lo: int, hi: int, qmax: int, g=lambda m: np.ones_like(m, float),
                  background=lambda m: 1.0 / m**2) -> tuple[np.ndarray, float]:
    """Least squares y(m) = g(m) Σ_{q' ≤ qmax} A_{q'} c_{q'}(m) + b·background(m) on m ∈ [lo, hi]; (A, R²)."""
    m = np.arange(lo, hi + 1)
    X = np.column_stack([g(m) * ramanujan(b, m) for b in range(1, qmax + 1)] + [background(m)])
    c, *_ = np.linalg.lstsq(X, y[lo - 1:hi], rcond=None)
    r = y[lo - 1:hi] - X @ c
    return c[:qmax], float(1 - r.var() / y[lo - 1:hi].var())


def brjuno(x: Fraction, T: int) -> float:
    """Yoccoz's Brjuno sum B(x) = Σ_{n≥0} β_{n−1} log(1/x_n) (Gauss map x_{n+1} = {1/x_n}, β_n = x_0 ⋯ x_n),
    truncated once the convergent denominator exceeds T; a rational x also stops at its last partial quotient."""
    x = x - x.numerator // x.denominator
    b, beta, q_prev, q_cur = 0.0, 1.0, 0, 1
    while x:
        b += beta * log(x.denominator / x.numerator)
        a = x.denominator // x.numerator
        q_prev, q_cur = q_cur, a * q_cur + q_prev
        if q_cur > T:
            break
        beta *= float(x)
        x = 1 / x - a
    return b
