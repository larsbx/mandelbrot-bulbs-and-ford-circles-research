"""P6: the parabolic index as a finite sum over the other cycles.

For the quadratic map at the root of B_{p/q}, f(z) = z² + c with
c = λ₀/2 − λ₀²/4 and λ₀ = e^{2πip/q}, the parabolic point z₀ = λ₀/2 is a
fixed point of F = f^q of multiplicity q + 1. Every other fixed point of F is
simple and repelling, and the holomorphic index formula on ℙ¹ (∞ is fixed with
multiplier 0, index 1) gives the exact identity

    ι_{p/q} = −Σ_{F(z) = z, z ≠ z₀} 1/(1 − F'(z)).

Ray bookkeeping. The 2^q − 1 external angles j/(2^q − 1) have period dividing
q under doubling. The q angles of the p/q rotation cycle land at z₀; the
remaining 2^q − 1 − q land one each on the 2^q − q − 1 other fixed points of F
[DH/Mil00, Gol92]. So pulling back every ray at once and following angle j to
its landing point z_j enumerates the terms of the sum, each exactly once, and
the doubling orbit of j is the cycle of z_j.

Numerics: rays are pulled back in double precision (8 sub-steps per halving of
the potential, nearest square-root branch), then polished by Newton on
F(z) = z. The output is VALIDATED evidence: the tests check it against the Arb
ball of `bulbford.index` and the ray relation f(z_j) = z_{2j}.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

import numpy as np

from bulbford.wake import rotation_cycle, wake

TAGS = ("DH-Mil00-rational-ray-landing", "Gol92-rotation-cycle-uniqueness", "Mil-12-holomorphic-index-formula")

SUBSTEPS = 8
OCTAVES = 56
LOG_RADIUS = np.log(1e4)


@dataclass(frozen=True)
class FixedPoints:
    p: int
    q: int
    c: complex
    z: np.ndarray  # z[j] = landing point of the ray of angle j/(2^q − 1)
    alpha_angles: np.ndarray  # the q angles landing at z₀

    def f(self, z: np.ndarray) -> np.ndarray:
        return z * z + self.c


@dataclass(frozen=True)
class CycleTerm:
    angles: tuple[int, ...]  # numerators over 2^q − 1, in orbit order
    period: int
    multiplier: complex  # (f^period)' along the cycle
    rotation: Fraction | None  # rotation number if the angles form a rotation set
    contribution: complex  # −period / (1 − multiplier^{q/period})


def parameter(p: int, q: int) -> complex:
    lam = np.exp(2j * np.pi * p / q)
    return lam / 2 - lam * lam / 4


def _pull_back_rays(c: complex, q: int, angles: list[int] | None = None, octaves: int = OCTAVES) -> np.ndarray:
    """Landing points of the rays of angles j/(2^q − 1), for a doubling-closed set of j.

    `angles` defaults to all j; a single doubling orbit is also closed, which makes
    one cycle cheap for large q. Angles are exact integers; only the starting
    points of the rays are rounded to floats.
    """
    M = 2**q - 1
    if angles is None:
        j = np.arange(M)
        doubled, frac = (2 * j) % M, j / M
    else:
        position = {a: i for i, a in enumerate(angles)}
        doubled = np.array([position[2 * a % M] for a in angles])
        frac = np.array([a / M for a in angles])
    theta = 2j * np.pi * frac
    ring = []
    for m in range(SUBSTEPS):
        w = np.exp(LOG_RADIUS * 2.0 ** (-m / SUBSTEPS) + theta)
        ring.append(w - c / (2 * w))  # inverse Böttcher map to O(w⁻³)
    for m in range(SUBSTEPS, SUBSTEPS * octaves):
        previous = ring[(m - 1) % SUBSTEPS]
        root = np.sqrt(ring[m % SUBSTEPS][doubled] - c)
        ring[m % SUBSTEPS] = np.where(np.abs(root - previous) <= np.abs(root + previous), root, -root)
    return ring[(SUBSTEPS * octaves - 1) % SUBSTEPS]


def _iterate(z: np.ndarray, c: complex, n: int) -> tuple[np.ndarray, np.ndarray]:
    w, dw = z.copy(), np.ones_like(z)
    for _ in range(n):
        dw = 2 * w * dw
        w = w * w + c
    return w, dw


def _polish(z: np.ndarray, c: complex, q: int, sweeps: int = 40) -> np.ndarray:
    for _ in range(sweeps):
        w, dw = _iterate(z, c, q)
        step = (w - z) / (dw - 1)
        z = z - step
        if np.max(np.abs(step)) < 1e-15:
            break
    return z


@lru_cache(maxsize=None)
def fixed_points(p: int, q: int) -> FixedPoints:
    c = parameter(p, q)
    M = 2**q - 1
    alpha = np.array(sorted(int(x * M) for x in rotation_cycle(p, q)))
    z = _pull_back_rays(c, q)
    others = np.setdiff1d(np.arange(M), alpha)
    z[others] = _polish(z[others], c, q)
    z[alpha] = np.exp(2j * np.pi * p / q) / 2
    return FixedPoints(p, q, c, z, alpha)


def flank_angles(p: int, q: int) -> tuple[int, int]:
    """α_{w−1} + 1 and α_{w+2} − 1 (numerators over 2^q − 1), exact.

    α_0 < … < α_{q−1} are the angles of the p/q rotation cycle and (α_w, α_{w+1})
    is the characteristic arc. These are the two angles just inside the outer
    ends of the three arcs around it.
    """
    M = 2**q - 1
    alpha = sorted(int(x * M) for x in rotation_cycle(p, q))
    w = alpha.index(int(wake(p, q)[0] * M))
    return (alpha[(w - 1) % q] + 1) % M, (alpha[(w + 2) % q] - 1) % M


def rotation_number(angles: tuple[int, ...], M: int) -> Fraction | None:
    """p'/k if doubling shifts the sorted angles cyclically by p' places, else None."""
    ordered = sorted(angles)
    k = len(ordered)
    position = {a: i for i, a in enumerate(ordered)}
    shifts = {(position[(2 * a) % M] - i) % k for i, a in enumerate(ordered)}
    return Fraction(shifts.pop(), k) if len(shifts) == 1 else None


def _orbit(j: int, M: int) -> tuple[int, ...]:
    orbit = [j]
    while (nxt := (2 * orbit[-1]) % M) != j:
        orbit.append(nxt)
    return tuple(orbit)


def _term(orbit: tuple[int, ...], z: np.ndarray, q: int, M: int) -> CycleTerm:
    k = len(orbit)
    rho = complex(np.prod(2 * z))
    return CycleTerm(orbit, k, rho, rotation_number(orbit, M), -k / (1 - rho ** (q // k)))


@lru_cache(maxsize=None)
def cycle_terms(p: int, q: int) -> tuple[CycleTerm, ...]:
    fp = fixed_points(p, q)
    M = 2**q - 1
    alpha = set(fp.alpha_angles.tolist())
    seen = np.zeros(M, dtype=bool)
    seen[fp.alpha_angles] = True
    terms = []
    for j in range(M):
        if seen[j]:
            continue
        orbit = _orbit(j, M)
        seen[list(orbit)] = True
        assert not alpha.intersection(orbit)
        terms.append(_term(orbit, fp.z[list(orbit)], q, M))
    return tuple(terms)


def cycle_through(p: int, q: int, angle: int) -> CycleTerm:
    """The term of the one cycle whose rays include angle/(2^q − 1), from its q rays alone.

    The angle must not land at z₀ (not in the p/q rotation cycle). A ray of period q
    needs several periods of pull-back to land, so the depth grows with q. The result
    is checked (the ray relation f(z_j) = z_{2j} and |F(z) − z|); a check that fails
    raises instead of returning an unverified term.
    """
    M, c = 2**q - 1, parameter(p, q)
    orbit = _orbit(angle % M, M)
    z = _polish(_pull_back_rays(c, q, list(orbit), octaves=max(OCTAVES, 6 * q)), c, q)
    position = {a: i for i, a in enumerate(orbit)}
    ray = max(abs(z[i] ** 2 + c - z[position[2 * a % M]]) for i, a in enumerate(orbit))
    residual = float(np.max(np.abs(_iterate(z, c, q)[0] - z)))
    if not (ray < 1e-9 and residual < 1e-9):
        raise ValueError(f"cycle through {angle}/{M} did not verify: ray {ray:.1e}, residual {residual:.1e}")
    return _term(orbit, z, q, M)


def index_by_cycles(p: int, q: int) -> complex:
    return complex(sum(t.contribution for t in cycle_terms(p, q)))
