"""Arithmetic of the internal angle p/q: modular inverse, x*, continued fractions.

Thin adapters over the vendored `rational_dynamics_py` (larsbx/finite-math-kernels,
vendor/python, pinned in vendored.toml); this module keeps the (p, q) signatures
the research code uses. Contracts that differ from the vendored functions:

- `coprime_numerators(1) == ()` here (the units 1 ≤ p < q); `units(1) == (0,)`.
- `modinv` and `xstar` refuse gcd(p, q) > 1 instead of reducing p/q first.
- `convergent_denominators` reads any valid expansion, canonical or not.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd

from rational_dynamics_py import (
    address,
    continued_fraction,
    from_continued_fraction,
    mod_inverse,
    signed_mod_inverse,
    units,
)
from rational_dynamics_py.farey import require_int


def _unit(p: int, q: int):
    """The address of p mod q, refusing a p that is not a unit mod q (q ≥ 2) and any non-int."""
    require_int(p, q)
    if q < 2 or gcd(p, q) != 1:
        raise ValueError(f"{p} is not a unit modulo {q}")
    return address(p % q, q)


def modinv(p: int, q: int) -> int:
    """p̄ ∈ (0, q) with p·p̄ ≡ 1 (mod q)."""
    return mod_inverse(_unit(p, q))


def xstar(p: int, q: int) -> Fraction:
    """x*(p/q) = ‖p⁻¹ mod q‖ / q ∈ (0, 1/2]."""
    return Fraction(abs(signed_mod_inverse(_unit(p, q))), q)


def cf(p: int, q: int) -> tuple[int, ...]:
    """Canonical continued fraction [a0; a1, …, an] of p/q (last quotient ≥ 2 unless q = 1)."""
    return continued_fraction(address(p, q))


def from_cf(a: tuple[int, ...]) -> tuple[int, int]:
    """(p, q) with p/q = [a0; a1, …, an], in lowest terms."""
    x = from_continued_fraction(a)
    return x.numerator, x.denominator


def convergent_denominators(a: tuple[int, ...]) -> tuple[int, ...]:
    """(q0, q1, …, qn), the denominators of the convergents of [a0; a1, …, an].

    Each convergent p_k/q_k of the recurrence is already in lowest terms
    (p_k q_{k−1} − p_{k−1} q_k = ±1), so q_k is the denominator of the prefix value.
    """
    return tuple(from_continued_fraction(a[: k + 1]).denominator for k in range(len(a)))


def coprime_numerators(q: int) -> tuple[int, ...]:
    """The units 1 ≤ p < q of ℤ/qℤ, increasing; empty for q = 1 (no bulb has denominator 1)."""
    require_int(q)
    return units(q) if q > 1 else ()
