"""Parabolic implosion at c = 1/4: the horn map of g(w) = w + w² and the 1/q roots.

In w = z − 1/2 the map z² + c is w ↦ w + w² + δ, δ = c − 1/4. At the root of
B_{1/q}, δ = −(1 − λ₀)²/4 = ε² with ε = (λ₀ − 1)/(2i), and the two fixed points
are w = ±iε (α = λ₀/2 above, β = 1 − λ₀/2 below).

Fatou coordinates of g: Φ(g(w)) = Φ(w) + 1 with the asymptotic expansion
Φ(w) = −1/w + log(∓w) + Σ_k c_k w^k (`fatou_coefficients`), log(−w) in the
attracting petal (w < 0) and log(w) in the repelling petal (w > 0). With that
normalisation, matching Φ_in and Φ_out to the perturbed coordinate
arctan(w/ε)/ε + ½ log(w² + ε²) gives Lavaurs' phase for n iterates as
σ = n − π/ε [Lav89, Shi00]: g_ε^n ≈ L_σ = Ψ_out ∘ (· + σ) ∘ Φ_in. For n = q at
the root of B_{1/q}, σ_q = q − π cot(π/q) + iπ (`lavaurs_phase`).

The horn map E = Φ_in ∘ Ψ_out satisfies E(Z + 1) = E(Z) + 1, and
E(Z) − Z = −iπ + Σ_{k ≥ 1} a_k e^{2πikZ} on the upper end (+iπ on the lower).
A fixed point w of L_σ (`lavaurs_fixed_point`, found in the w-plane) gives
Z = Φ_in(w) + σ with E(Z) − Z = −σ, and L_σ'(w) = E'(Z). The q → ∞ phase of the
1/q roots is σ = iπ, the constant of the upper end: the parabolic cycle sits at
that end. `kappa0` is a₂/(2πi a₁²): the second-order coefficient
of the multiplier of the fixed point leaving the upper end, E' = 1 − u + κ₀u² + …,
in the bulb coordinate u = 2πi(σ − iπ). The ratio is invariant under both
translations of the Fatou coordinates.

Arithmetic is mpmath at `DPS` digits; the orbits run until |w| < `RADIUS`, where
the expansion to order `ORDER` has truncation error below RADIUS^(ORDER+1).
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import comb

import mpmath as mp

TAGS = ("Lav89-Lavaurs-phase", "Shi00-parabolic-implosion", "EV-horn-map")

DPS = 40
ORDER = 12
RADIUS = mp.mpf("1e-3")
ESCAPE = 10
STEPS = 200_000


@lru_cache(maxsize=None)
def fatou_coefficients(order: int = ORDER) -> tuple[Fraction, ...]:
    """c_1 … c_order in Φ(w) = −1/w + log(∓w) + Σ c_k w^k, exact.

    Φ(w + w²) − Φ(w) − 1 = 1/(1 + w) + log(1 + w) − 1 + Σ c_k ((w + w²)^k − w^k),
    and (w + w²)^k − w^k = k w^{k+1} + …, so the w^{k+1} coefficient fixes c_k.
    """
    residual = [Fraction(0)] + [Fraction((-1) ** j) - Fraction((-1) ** j, j) for j in range(1, order + 2)]
    coefficients = []
    for k in range(1, order + 1):
        ck = -residual[k + 1] / k
        for i in range(1, order + 2 - k):
            residual[k + i] += ck * comb(k, i)
        coefficients.append(ck)
    return tuple(coefficients)


def _expansion(w, outgoing: bool):
    cs = [mp.mpf(c.numerator) / c.denominator for c in fatou_coefficients()]
    log = mp.log(w) if outgoing else mp.log(-w)
    value = -1 / w + log + sum(c * w ** (k + 1) for k, c in enumerate(cs))
    slope = 1 / w**2 + 1 / w + sum((k + 1) * c * w**k for k, c in enumerate(cs))
    return value, slope


def phi_in(w):
    """(Φ_in(w), Φ_in'(w)) for w in the basin of 0; ValueError if the orbit escapes."""
    with mp.workdps(DPS):
        w, dw, n = mp.mpc(w), mp.mpc(1), 0
        while not (abs(w) < RADIUS and abs(w.imag) < -w.real):
            dw, w, n = dw * (1 + 2 * w), w + w * w, n + 1
            if abs(w) > ESCAPE or n > STEPS:
                raise ValueError("the orbit escapes or lingers: not verified in the parabolic basin")
        value, slope = _expansion(w, outgoing=False)
        return value - n, slope * dw


def psi_out(Z):
    """(Ψ_out(Z), Ψ_out'(Z)): g^m of the repelling-petal inverse at Z − m."""
    with mp.workdps(DPS):
        Z = mp.mpc(Z)
        m = max(0, int(mp.ceil(Z.real + 1 / RADIUS + 2 * abs(Z.imag))))
        target, w = Z - m, -1 / (Z - m)
        for _ in range(100):
            value, slope = _expansion(w, outgoing=True)
            step = (value - target) / slope
            w -= step
            if abs(step) < mp.mpf(10) ** (5 - DPS):
                break
        dw = 1 / _expansion(w, outgoing=True)[1]
        for _ in range(m):
            dw, w = dw * (1 + 2 * w), w + w * w
        return w, dw


def horn(Z):
    """(E(Z) − Z, E'(Z)) for the horn map E = Φ_in ∘ Ψ_out."""
    with mp.workdps(DPS):
        w, dw = psi_out(Z)
        value, slope = phi_in(w)
        return value - Z, slope * dw


def lavaurs_phase(q: int):
    """σ_q = q − π/ε at the root of B_{1/q}: q − π cot(π/q) + iπ."""
    with mp.workdps(DPS):
        return q - mp.pi * mp.cot(mp.pi / q) + 1j * mp.pi


def horn_coefficients(kmax: int, height: float = 3.0, nodes: int = 48) -> list:
    """a_1 … a_kmax of E(Z) − Z + iπ = Σ a_k e^{2πikZ}, sampled on Im Z = height."""
    with mp.workdps(DPS):
        y = mp.mpf(height)
        h = [horn(mp.mpc(mp.mpf(j) / nodes, y))[0] + 1j * mp.pi for j in range(nodes)]
        return [sum(h[j] * mp.expjpi(-2 * k * mp.mpf(j) / nodes) for j in range(nodes)) / nodes
                * mp.exp(2 * mp.pi * k * y) for k in range(1, kmax + 1)]


def kappa0(height: float = 3.0):
    with mp.workdps(DPS):
        a1, a2 = horn_coefficients(2, height)
        return a2 / (2j * mp.pi * a1**2)


def lavaurs_map(w, sigma=None):
    """(L_σ(w), L_σ'(w)) for L_σ = Ψ_out ∘ (· + σ) ∘ Φ_in; default σ = iπ."""
    with mp.workdps(DPS):
        sigma = 1j * mp.pi if sigma is None else mp.mpc(sigma)
        value, slope = phi_in(w)
        image, dimage = psi_out(value + sigma)
        return image, dimage * slope


def lavaurs_fixed_point(w, sigma=None, tol: float = 1e-25, cap: float = 0.2):
    """(w*, L_σ'(w*)) with L_σ(w*) = w*, by Newton from w with steps capped at `cap`.

    L_σ'(w*) = E'(Z*) at Z* = Φ_in(w*) + σ, the matching fixed point of the horn map.
    """
    with mp.workdps(DPS):
        w = mp.mpc(w)
        for _ in range(100):
            image, slope = lavaurs_map(w, sigma)
            step = (image - w) / (slope - 1)
            w -= step * min(1, cap / abs(step))
            if abs(step) < tol:
                return w, lavaurs_map(w, sigma)[1]
        raise ValueError(f"Newton did not converge from {w}")


def fixed_point_of_cycle(z, sigma=None, seeds: int = 4):
    """The fixed points of L_σ seeded by the `seeds` points of a p = 1 cycle farthest from 1/2.

    Returns the list of (w*, L_σ'(w*)), one per seed that converges; for a cycle that
    L_σ models, they coincide.
    """
    with mp.workdps(DPS):
        w = sorted((mp.mpc(complex(x)) - mp.mpf(1) / 2 for x in z), key=abs, reverse=True)
        found = []
        for x in w[:seeds]:
            try:
                found.append(lavaurs_fixed_point(x, sigma))
            except ValueError:
                pass
        return found
