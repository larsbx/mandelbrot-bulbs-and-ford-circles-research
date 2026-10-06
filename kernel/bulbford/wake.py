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

The generic finite arithmetic — `rotation_cycle`, `mechanical`, `wake`, `farey` —
is the vendored `rational_dynamics_py` (larsbx/finite-math-kernels, vendor/python,
pinned in vendored.toml); the names here are its adapters, kept because
larsbx/math-vizops loads this file on its own from a sibling checkout and reads
`wake`, `mechanical`, `rotation_cycle` and `double`. This repository keeps the
TAGS, `heights` (P8) and the finite check `acts_as_rotation`.
"""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

if __package__ != "bulbford":  # loaded as a bare file (math-vizops): use this checkout's bulbford
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import bulbford  # noqa: E402,F401  (puts the pinned vendor/python first and refuses any other copy)
import rational_dynamics_py as _rd  # noqa: E402

TAGS = ("Gol92-rotation-cycle-uniqueness", "DH-Mil00-rational-ray-landing")


def double(theta: Fraction) -> Fraction:
    """2θ mod 1 on any Fraction (kept local: the vendored `double_mod_one` takes a
    non-negative reduced `Address`, and math-vizops calls this on Fractions)."""
    return 2 * theta % 1


def rotation_cycle(p: int, q: int) -> tuple[Fraction, ...]:
    """The doubling cycle x₀ < … < x_{q−1} of rotation number p/q (refuses p/q not reduced in (0, 1))."""
    return _rd.rotation_cycle(p, q)


def mechanical(p: int, q: int, r: int) -> int:
    """The conjugate c(r) of the rotation cycle as a numerator over 2^q − 1 (P7).

    Bit k (most significant first) is [(r + k·p) mod q ≥ q − p]. c(0) is the word
    `rotation_cycle` starts from, doubling sends c(r) to c(r + p), and the sorted
    cycle is c(0) < c(1) < … < c(q − 1). Refuses p/q not reduced in (0, 1).
    """
    return _rd.mechanical_word(p, q, r)


def heights(word: int, q: int) -> list[int]:
    """g_k = q·H_k − k·m for the cyclic q-bit word (most significant bit first), k < q (P8).

    H_k is the number of 1s among the first k bits and m the weight. Rotating the
    word by t translates the set {g_k} by −g_t, so that set up to translation is a
    doubling invariant.
    """
    bits = [int(b) for b in format(word, f"0{q}b")]
    m, out, h = sum(bits), [], 0
    for k, bit in enumerate(bits):
        out.append(q * h - k * m)
        h += bit
    return out


def acts_as_rotation(p: int, q: int, cycle: tuple[Fraction, ...]) -> bool:
    """Finite check: doubling maps x_i to x_{i+p mod q} for every i."""
    return all(double(x) == cycle[(i + p) % q] for i, x in enumerate(cycle))


def wake(p: int, q: int) -> tuple[Fraction, Fraction]:
    """(θ₋, θ₊): the consecutive cycle points bounding the shortest arc, (x_{p−1}, x_p)."""
    return _rd.wake(p, q)


def farey(n: int) -> tuple[Fraction, ...]:
    """Interior Farey fractions of order n, increasing (the interior of F_n; n < 1 is refused)."""
    return _rd.farey_sequence(n, interior=True)
