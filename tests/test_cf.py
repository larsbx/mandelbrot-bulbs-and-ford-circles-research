from fractions import Fraction
from bulbford.cf import cf, from_cf, xstar, modinv, convergent_denominators, coprime_numerators


def test_cf_roundtrip():
    for q in range(2, 60):
        for p in coprime_numerators(q):
            a = cf(p, q)
            assert from_cf(a) == (p, q)
            assert len(a) == 1 or a[-1] >= 2


def test_xstar_is_previous_convergent_denominator():
    """x* = q_{n-1}/q for the canonical CF (a_n ≥ 2)."""
    for q in range(2, 60):
        for p in coprime_numerators(q):
            qs = convergent_denominators(cf(p, q))
            assert qs[-1] == q
            assert xstar(p, q) == Fraction(qs[-2], q)


def test_modinv():
    assert modinv(3, 7) == 5 and (3 * 5) % 7 == 1
