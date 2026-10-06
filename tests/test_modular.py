"""V57: exact normal-form invariants in ℚ(ζ_q) by multimodular arithmetic (bulbford/modular.py)."""
from fractions import Fraction
from math import gcd, lcm
import mpmath as mp
import pytest
from bulbford.modular import exact_invariants, embed
from bulbford.normal_form import normal_form

# data/exact_denominators.txt (exact ℚ(ζ_q) series arithmetic, scripts/exact_index.py), power basis 1, ζ, …, λ₀ = ζ
A_EXACT = {3: (5, 1), 4: (-10, -6), 5: (36, 22, 22, 10), 6: (64, -116), 7: (220, 163, 258, 247, 276, 124),
           8: (384, 520, 108, -708)}
DENOM_IOTA_A2 = {3: 3, 4: 2, 5: 5, 6: 3, 7: 7, 8: 1}


@pytest.mark.parametrize("q", sorted(A_EXACT))
def test_reproduces_exact_small_q(q):
    """a = [z^{q+1}] f^q and denom(ι·a²) = denom(X/q) agree with the exact series computation for q ≤ 8."""
    e = exact_invariants(q)
    assert e.a == A_EXACT[q]
    assert lcm(*(Fraction(x, q).denominator for x in e.X)) == DENOM_IOTA_A2[q]


@pytest.mark.parametrize("q", [7, 12, 30])
def test_embeddings_are_the_galois_orbit(q):
    """σ_j(X) = q·ι_{j/q}·a(j/q)² for every j coprime to q, against the floating normal form at 60 digits."""
    e = exact_invariants(q)
    for j in (j for j in range(1, q) if gcd(j, q) == 1):
        nf = normal_form(j, q, dps=60)
        with mp.workdps(60):
            assert abs(embed(e.a, q, j) - nf.a) <= mp.mpf(10) ** -40 * abs(nf.a)
            ref = q * nf.iota * nf.a ** 2
            assert abs(embed(e.X, q, j) - ref) <= mp.mpf(10) ** -40 * abs(ref)


def _rad_odd(q):
    """Product of the odd primes dividing q."""
    out, l = 1, 3
    while q % 2 == 0:
        q //= 2
    while q > 1:
        if q % l == 0:
            out *= l
            while q % l == 0:
                q //= l
        l += 2
    return out


@pytest.mark.parametrize("q", range(2, 41))
def test_integrality_and_denominator_law(q):
    """C22 in sharp form (V57, validated for q ≤ 200, 256, 384): q·ι·a² ∈ ℤ[ζ_q], and the least d with d·ι·a² ∈ ℤ[ζ_q]
    is the product of the odd primes dividing q, doubled for q = 2, 4."""
    e = exact_invariants(q)
    assert all(x.denominator == 1 for x in e.X)
    assert lcm(*(Fraction(x, q).denominator for x in e.X)) == _rad_odd(q) * (2 if q in (2, 4) else 1)
