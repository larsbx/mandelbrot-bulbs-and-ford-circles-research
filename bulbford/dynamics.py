"""Cycles, multipliers and bulb-size proxies for unicritical families z ↦ z^d + c.

A `Family` is one hyperbolic component H parametrised by its multiplier λ ∈ 𝔻:
c = c_of(λ). Its satellite at internal angle p/q has period q·period and root
c_of(λ₀), λ₀ = e^{2πip/q}. All state is carried in immutable `Cycle` records.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import pi
from typing import Callable, Iterable, NamedTuple
import cmath
import numpy as np

try:                                   # compiled kernel when numba is present; identical code runs as Python otherwise
    from numba import njit
except ImportError:                    # pragma: no cover
    def njit(*args, **kwargs):
        return args[0] if args and callable(args[0]) else (lambda f: f)


@dataclass(frozen=True)
class Family:
    name: str
    d: int
    period: int
    c_of: Callable[[complex], complex]
    dc_of: Callable[[complex], complex]     # c'(λ)


MAIN2 = Family("main2", 2, 1, lambda l: l / 2 - l * l / 4, lambda l: (1 - l) / 2)
DISK2 = Family("disk2", 2, 2, lambda l: l / 4 - 1, lambda l: 0.25 + 0j)
MAIN3 = Family("main3", 3, 1,
               lambda l: (lambda z: z - z ** 3)(cmath.sqrt(l / 3)),
               lambda l: (1 - l) / (6 * cmath.sqrt(l / 3)))
FAMILIES = {f.name: f for f in (MAIN2, DISK2, MAIN3)}


class Orbit(NamedTuple):
    """f_c^n(z) and its partials: z_z = ∂/∂z, z_c = ∂/∂c, z_zz, z_zc."""
    z: complex
    z_z: complex
    z_c: complex
    z_zz: complex
    z_zc: complex


@njit(cache=True)
def _orbit_kernel(d, c, z, n):
    """f^n(z), f(z) = z^d + c, with ∂z, ∂c, ∂zz, ∂zc — the hot loop of every cycle computation (compiled by numba)."""
    z_z, z_c, z_zz, z_zc = 1 + 0j, 0j, 0j, 0j
    for _ in range(n):
        zd2 = 1 + 0j                                  # z^{d−2}
        for _ in range(d - 2):
            zd2 *= z
        w, w1 = d * zd2 * z, d * (d - 1) * zd2       # f'(z), f''(z)
        z, z_z, z_c, z_zz, z_zc = (zd2 * z * z + c, w * z_z, w * z_c + 1,
                                   w1 * z_z * z_z + w * z_zz, w1 * z_z * z_c + w * z_zc)
    return z, z_z, z_c, z_zz, z_zc


def orbit(fam: Family, c: complex, z: complex, n: int) -> Orbit:
    return Orbit(*_orbit_kernel(fam.d, complex(c), complex(z), int(n)))


def newton(step: Callable[[complex], complex], x0: complex, tol: float = 1e-15, maxit: int = 80) -> complex:
    x = x0
    for _ in range(maxit):
        s = step(x)
        x = x - s
        if abs(s) < tol * max(1.0, abs(x)):
            break
    return x


class Cycle(NamedTuple):
    """One point z of a period-n cycle of z ↦ z^d + c, with its multiplier ρ and dρ/dc."""
    c: complex
    z: complex
    n: int
    rho: complex
    drho_dc: complex


def cycle_at(fam: Family, c: complex, z0: complex, n: int) -> Cycle:
    z = newton(lambda z: (lambda o: (o.z - z) / (o.z_z - 1))(orbit(fam, c, z, n)), z0)
    o = orbit(fam, c, z, n)
    dz_dc = -o.z_c / (o.z_z - 1)                    # implicit differentiation of f^n(z) = z
    return Cycle(c, z, n, o.z_z, o.z_zc + o.z_zz * dz_dc)


def center(fam: Family, c0: complex, n: int) -> complex:
    """Newton on f_c^n(0) = 0 from c0."""
    return newton(lambda c: (lambda o: o.z / o.z_c)(orbit(fam, c, 0j, n)), c0)


def has_exact_period(fam: Family, cyc: Cycle, tol: float = 1e-9) -> bool:
    """z is not a fixed point of f^m for a proper divisor m of n (guards Newton against period collapse)."""
    return all(abs(orbit(fam, cyc.c, cyc.z, m).z - cyc.z) > tol * max(1.0, abs(cyc.z))
               for m in range(1, cyc.n) if cyc.n % m == 0)


def step(fam: Family, cyc: Cycle, c: complex, depth: int = 0) -> Cycle:
    """Continue the cycle to parameter c, bisecting the step whenever the period collapses."""
    nxt = cycle_at(fam, c, cyc.z, cyc.n)
    if has_exact_period(fam, nxt) or depth >= 20:
        return nxt
    return step(fam, step(fam, cyc, (cyc.c + c) / 2, depth + 1), c, depth + 1)


def track(fam: Family, cyc: Cycle, path: Iterable[complex]) -> Cycle:
    """Continue the cycle along a parameter path (last point returned)."""
    for c in path:
        cyc = step(fam, cyc, c)
    return cyc


def track_all(fam: Family, cyc: Cycle, path: Iterable[complex]) -> tuple[Cycle, ...]:
    out: tuple[Cycle, ...] = ()
    for c in path:
        cyc = step(fam, cyc, c)
        out = out + (cyc,)
    return out


def solve_rho(fam: Family, cyc: Cycle, target: complex, tol: float = 1e-14) -> Cycle:
    """Newton in c on ρ(c) = target, cycle tracked at each step."""
    for _ in range(60):
        s = (cyc.rho - target) / cyc.drho_dc
        cyc = cycle_at(fam, cyc.c - s, cyc.z, cyc.n)
        if abs(s) < tol * abs(cyc.c):
            break
    return cyc


class Bulb(NamedTuple):
    fam: Family
    p: int
    q: int
    root: complex
    cen: Cycle          # ρ = 0 (critical point on the cycle: ρ vanishes to order d − 1 here)
    ants: tuple         # the d − 1 boundary points with ρ = −1 (ρ is a (d−1)-fold branched cover of the disc)

    @property
    def ant(self) -> Cycle:
        return self.ants[0]

    @property
    def lam0(self) -> complex:
        return cmath.exp(2j * pi * self.p / self.q)

    @property
    def pred(self) -> float:
        """first-order diameter 2|c'(λ₀)| q⁻²."""
        return 2 * abs(self.fam.dc_of(self.lam0)) / self.q ** 2

    @property
    def G_ants(self) -> tuple:
        return tuple(abs(a.c - self.root) / self.pred for a in self.ants)

    @property
    def G_ant(self) -> float:
        return self.G_ants[0]

    @property
    def G_cen(self) -> float:
        return 2 * abs(self.cen.c - self.root) / self.pred


def bulb(fam: Family, p: int, q: int, steps: int = 40) -> Bulb:
    """Centre (ρ = 0) and all d − 1 antipodes (ρ = −1) of the satellite B_{p/q}.  Near the centre
    ρ ≈ A(c − c_cen)^{d−1}; each branch is seeded at ρ = −ε and continued in the target ρ* ∈ [−ε, −1], so the
    solution stays in the bulb (the root→centre ray need not pass an antipode in degree ≥ 3)."""
    lam0 = cmath.exp(2j * pi * p / q)
    n = q * fam.period
    root = fam.c_of(lam0)
    c_cen = center(fam, fam.c_of(lam0 * (1 + 1 / q ** 2)), n)
    cen = cycle_at(fam, c_cen, 0j, n)
    k, eps = fam.d - 1, 1e-3
    delta = 1e-3 * abs(c_cen - root)
    A = cycle_at(fam, c_cen + delta, cen.z, n).rho / delta ** k
    ants = ()
    for j in range(k):
        c0 = c_cen + ((-eps / A) ** (1 / k)) * cmath.exp(2j * pi * j / k)
        cyc = solve_rho(fam, cycle_at(fam, c0, cen.z, n), -eps)
        for target in np.linspace(-eps, -1, steps + 1)[1:]:
            cyc = solve_rho(fam, cyc, target)
        ants = ants + (cyc,)
    return Bulb(fam, p, q, root, cen, ants)


def rho_on_path(fam: Family, p: int, q: int, us: np.ndarray) -> np.ndarray:
    """ρ at c = c_of(λ₀ e^{u/q²}) along the u-path `us` (tracked continuously from the centre)."""
    lam0 = cmath.exp(2j * pi * p / q)
    n = q * fam.period
    cen = cycle_at(fam, center(fam, fam.c_of(lam0 * (1 + 1 / q ** 2)), n), 0j, n)
    c_of_u = lambda u: fam.c_of(lam0 * cmath.exp(u / q ** 2))
    lead_in = (cen.c + (c_of_u(us[0]) - cen.c) * t for t in np.linspace(0, 1, 9)[1:])
    start = track(fam, cen, lead_in)
    return np.array([cy.rho for cy in track_all(fam, start, map(c_of_u, us))])
