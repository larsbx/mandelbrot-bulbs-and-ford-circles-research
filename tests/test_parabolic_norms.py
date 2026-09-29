"""The parabolic coefficient a_q at the solvable multipliers, and the norm-divisibility lemma (B3)."""
from math import gcd

import pytest
from sympy import cyclotomic_poly, factorint, n_order

from bulbford.norms import a_at_2, a_at_4, a_coefficient, norm_a

# Exact resultants Res(Φ_q, a_q) from sympy (bridges spike), for the replay of norm_a.
NORMS = {2: 2, 3: 21, 4: 136, 5: 845680, 6: 10128, 7: 206087845452517, 8: 459109339648,
         9: 238356750916093565376, 10: 9764277910343680}


@pytest.mark.parametrize("q", range(1, 25))
def test_power_map_closed_form(q):
    assert a_coefficient(2, q) == a_at_2(q)


@pytest.mark.parametrize("q", range(1, 12))
def test_chebyshev_closed_form(q):
    assert a_coefficient(4, q) == a_at_4(q)


@pytest.mark.parametrize("q", sorted(NORMS))
def test_norm_matches_resultant(q):
    assert norm_a(q) == NORMS[q]


def lemma_rows(q):
    """(ℓ, zeros, v_ℓ(N(a_q))) for primes ℓ > 2q + 2 dividing Φ_q(2)Φ_q(4)Φ_q(−2)."""
    n = norm_a(q)
    primes = {l for b in (2, 4, -2) for l in factorint(abs(int(cyclotomic_poly(q, b)))) if l > 2 * q + 2}
    for l in sorted(primes):
        zeros = {x % l for x in (2, 4, -2) if n_order(x % l, l) == q}
        v = next(v for v in range(64) if n % l ** (v + 1))
        yield l, zeros, v


@pytest.mark.parametrize("q", range(3, 15))
def test_norm_divisibility(q):
    """Proven for λ ∈ {2, 4}; λ = −2 is validated (the equality is what q ≤ 20 shows)."""
    for l, zeros, v in lemma_rows(q):
        assert v == len(zeros), (q, l, zeros, v)
        for lam in zeros:
            assert a_coefficient(lam, q) % l == 0


def test_mod_l_reduction_is_the_solvable_map():
    """For ℓ ≡ 1 (mod q), reduction ζ_q ↦ 2 (mod ℓ) is a ring map, so a_q(ζ) ≡ a_q(2) = C(2^q, q+1)."""
    for q, l in ((5, 31), (7, 127), (13, 8191)):
        assert n_order(2, l) == q and (l - 1) % q == 0 and gcd(a_at_2(q), l) == l
