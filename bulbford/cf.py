"""Arithmetic of the internal angle p/q: modular inverse, x*, continued fractions."""
from __future__ import annotations
from fractions import Fraction
from math import gcd


def modinv(p: int, q: int) -> int:
    """p̄ ∈ (0, q) with p·p̄ ≡ 1 (mod q)."""
    return pow(p, -1, q)


def xstar(p: int, q: int) -> Fraction:
    """x*(p/q) = ‖p⁻¹ mod q‖ / q ∈ (0, 1/2]."""
    pb = modinv(p, q)
    return Fraction(min(pb, q - pb), q)


def cf(p: int, q: int) -> tuple[int, ...]:
    """Canonical continued fraction [a0; a1, …, an] of p/q (last quotient ≥ 2 unless q = 1)."""
    return (p // q,) + (cf(q, p % q) if p % q else ())


def from_cf(a: tuple[int, ...]) -> tuple[int, int]:
    """(p, q) with p/q = [a0; a1, …, an]."""
    if len(a) == 1:
        return a[0], 1
    p, q = from_cf(a[1:])          # tail = a1 + 1/(a2 + …) = p/q
    return a[0] * p + q, p


def convergent_denominators(a: tuple[int, ...]) -> tuple[int, ...]:
    """q_{-1}=0? no: (q0, q1, …, qn) of the convergents of [a0; a1, …, an]."""
    qs: tuple[int, ...] = (1,)
    prev = 0
    for ak in a[1:]:
        qs, prev = qs + (ak * qs[-1] + prev,), qs[-1]
    return qs


def coprime_numerators(q: int) -> tuple[int, ...]:
    return tuple(p for p in range(1, q) if gcd(p, q) == 1)
