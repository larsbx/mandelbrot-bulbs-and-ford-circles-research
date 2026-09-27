"""P9, the arc-width lemma: the arc j places from the characteristic arc has width 2^{(j·p̄) mod q}/M."""

from math import gcd

import pytest

from bulbford.wake import mechanical

FRACTIONS = [(p, q) for q in range(3, 61) for p in range(1, q) if gcd(p, q) == 1]


def exponents(p: int, q: int) -> dict[int, int]:
    """j ↦ log₂ of the width (in units 1/M) of the arc from α_{w+j} to α_{w+j+1}, indices mod q."""
    M, w = 2**q - 1, p - 1
    A = [mechanical(p, q, r) for r in range(q)]
    widths = {j: (A[(w + j + 1) % q] - A[(w + j) % q]) % M for j in range(q)}
    assert all(x & (x - 1) == 0 for x in widths.values())  # every width is a power of two
    return {j: x.bit_length() - 1 for j, x in widths.items()}


@pytest.mark.parametrize("p,q", FRACTIONS)
def test_p9_arc_widths_are_placed_by_pbar(p, q):
    pbar = pow(p, -1, q)
    assert exponents(p, q) == {j: j * pbar % q for j in range(q)}


@pytest.mark.parametrize("p,q", FRACTIONS)
def test_the_widest_arcs_sit_at_minus_p_and_minus_two_p(p, q):
    e = exponents(p, q)
    assert e[(-p) % q] == q - 1 and e[(-2 * p) % q] == q - 2



@pytest.mark.parametrize("p,q", FRACTIONS)
def test_a_flanking_arc_is_among_the_two_widest_iff_xstar_at_most_two_over_q(p, q):
    e, pbar = exponents(p, q), pow(p, -1, q)
    assert (e[1], e[q - 1]) == (pbar, q - pbar)
    assert (max(e[1], e[q - 1]) >= q - 2) == (min(pbar, q - pbar) <= 2)
