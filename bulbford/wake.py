"""Exact rotation cycles of angle doubling and the combinatorial p/q-wake.

For 0 < p/q < 1 in lowest terms, the doubling map θ ↦ 2θ mod 1 has a cycle
x₀ < … < x_{q−1} on which it acts as the rotation x_i ↦ x_{i+p mod q}.  Its
angles are rationals with denominator 2^q − 1, built from the mechanical word
b_k = [k·p mod q ≥ q − p].  The shortest arc between consecutive cycle points,
(θ₋, θ₊), has length 1/(2^q − 1); it is the characteristic arc.

Everything here is finite arithmetic in ℚ/ℤ, the same doubling that
`larsbx/finite-math-kernels` R1 owns.  Two readings are imports, not
consequences:

  [Gol92]  the rotation cycle with rotation number p/q is unique;
  [DH/Mil00] the parameter rays at θ₋, θ₊ land at the root of B_{p/q}, so the
             characteristic arc is the angular width of the p/q-wake.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd

TAGS = ("Gol92-rotation-cycle-uniqueness", "DH-Mil00-rational-ray-landing")


def double(theta: Fraction) -> Fraction:
    return 2 * theta % 1


def _reduced(p: int, q: int) -> None:
    if not (0 < p < q and gcd(p, q) == 1):
        raise ValueError("need 0 < p < q with gcd(p, q) = 1")


def rotation_cycle(p: int, q: int) -> tuple[Fraction, ...]:
    """The doubling cycle x₀ < … < x_{q−1} of rotation number p/q."""
    _reduced(p, q)
    bits = "".join("1" if k * p % q >= q - p else "0" for k in range(q))
    x0 = Fraction(int(bits, 2), 2**q - 1)
    orbit = [x0]
    for _ in range(q - 1):
        orbit.append(double(orbit[-1]))
    return tuple(sorted(orbit))


def acts_as_rotation(p: int, q: int, cycle: tuple[Fraction, ...]) -> bool:
    """Finite check: doubling maps x_i to x_{i+p mod q} for every i."""
    return all(double(x) == cycle[(i + p) % q] for i, x in enumerate(cycle))


def wake(p: int, q: int) -> tuple[Fraction, Fraction]:
    """(θ₋, θ₊): the consecutive cycle points bounding the shortest arc."""
    cycle = rotation_cycle(p, q)
    return min(zip(cycle, cycle[1:]), key=lambda pair: pair[1] - pair[0])


def farey(n: int) -> tuple[Fraction, ...]:
    """Interior Farey fractions of order n, increasing."""
    return tuple(sorted({Fraction(p, q) for q in range(2, n + 1) for p in range(1, q)}))
