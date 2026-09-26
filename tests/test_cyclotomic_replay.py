"""Replay finite-math-kernels cyclotomic-germ/v1 vectors against two independent local instruments.

The kernel's `[w^q] 1/P` for g(w) = ζ_q^p w + w² is, in this repository's reading, the holomorphic index
ι_{p/q} (bulbford/index.py docstring).  The vectors are checked against

  * scripts/exact_index.py — separately written ℚ(ζ_q) series code (sympy for the field inverse);
  * bulbford.index.index    — Arb ball arithmetic, with a rigorous radius;
  * data/parabolic_index_exact_vectors.json — the pinned Q(i) fixtures.
"""

import cmath
import importlib.util
import json
from fractions import Fraction as F
from math import gcd, pi
from pathlib import Path

import pytest

from bulbford.index import index

ROOT = Path(__file__).resolve().parents[1]
PAYLOAD = json.loads((ROOT / "data" / "cyclotomic_germ_v1_vectors.json").read_text())
VECTORS = PAYLOAD["vectors"]
SOURCE_COMMIT = "16bf936f308303a0012ceb0bc388400db6595d46"  # larsbx/finite-math-kernels, fixtures/cyclotomic_germ_v1.json

_spec = importlib.util.spec_from_file_location("exact_index", ROOT / "scripts" / "exact_index.py")
exact_index = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(exact_index)


def coords(row) -> tuple[F, ...]:
    return tuple(F(c) for c in row["reciprocal_coefficient"])


def evaluate(v: tuple[F, ...], q: int) -> complex:
    """Cross-check only: numerical evaluation at ζ_q = e^{2πi/q} (never a certificate)."""
    return sum(complex(c) * cmath.exp(2j * pi * k / q) for k, c in enumerate(v))


def test_payload_shape():
    assert PAYLOAD["schema"] == "cyclotomic-germ/v1"
    assert PAYLOAD["authority"] == "none"
    assert len(SOURCE_COMMIT) == 40
    expected = {(1, 1)} | {(p, q) for q in range(2, 9) for p in range(1, q) if gcd(p, q) == 1}
    assert {(r["p"], r["q"]) for r in VECTORS} == expected


@pytest.mark.parametrize("row", [r for r in VECTORS if r["q"] > 1], ids=lambda r: f'{r["p"]}/{r["q"]}')
def test_vectors_match_independent_exact_series(row):
    assert tuple(exact_index.exact_index(row["p"], row["q"])) == coords(row)


@pytest.mark.parametrize("row", VECTORS, ids=lambda r: f'{r["p"]}/{r["q"]}')
def test_vectors_lie_in_the_arb_ball(row):
    ball, value = index(row["p"], row["q"]), evaluate(coords(row), row["q"])
    assert abs(complex(ball) - value) <= float(ball.rad()) + 1e-13


def test_vectors_agree_with_pinned_gaussian_fixtures():
    fixtures = json.loads((ROOT / "data" / "parabolic_index_exact_vectors.json").read_text())["vectors"]
    kernel = {(r["p"], r["q"]): coords(r) for r in VECTORS}
    for f in fixtures:
        re, im = F(f["re_num"], f["re_den"]), F(f["im_num"], f["im_den"])
        # basis of Q(ζ_q): q ∈ {1, 2} has φ = 1; q = 4 has (1, i)
        assert kernel[(f["p"], f["q"])] == ((re,) if f["q"] < 4 else (re, im))


def test_conjugate_rotation_gives_conjugate_value():
    """ι_{(q−p)/q} = conj ι_{p/q}: the automorphism ζ ↦ ζ^{−1} (kernel C2), read numerically here."""
    kernel = {(r["p"], r["q"]): coords(r) for r in VECTORS}
    for (p, q), v in kernel.items():
        if q > 2:
            assert abs(evaluate(kernel[(q - p, q)], q) - evaluate(v, q).conjugate()) < 1e-13
