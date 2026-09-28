"""Parabolic implosion for the two V25 families: the horn maps at c = 1/4 and c = −3/4.

A `Germ` is the parabolic point that the roots of a family approach, in the local
coordinate w = z − z₀, with f(w) = w·(w + 2z₀) the map z² + c₀ there:

- `ONE` (p = 1): c₀ = 1/4, z₀ = 1/2, f(w) = w + w², one attracting petal (w < 0)
  and one repelling petal (w > 0).
- `HALF` (p = (q − 1)/2): c₀ = −3/4, z₀ = −1/2, f(w) = w² − w with multiplier −1.
  f² = w − 2w³ + w⁴ has two attracting petals (w ≷ 0) and two repelling ones
  (Im w ≷ 0), and f swaps each pair.

Fatou coordinates solve Φ(f(w)) = Φ(w) + `step` (1, or ½ for the half-step of
`HALF`). They are normalised by the expansion Φ = Σ_{e<0} a_e w^e + c·log(±w^m) +
Σ_k d_k w^k (`fatou_coefficients`, exact), and used only at |w| < `radius`, where
the truncation error is below (2·radius)^(order+1). Φ_in is taken on the
incoming petal (`incoming` direction), and Ψ_out = Φ_out⁻¹ is taken on the
outgoing one (`outgoing`). For `HALF` the other petal of each pair is reached by
one half-step, and Φ on it is a different function. So `phi_in` also returns
the parity of the half-steps it took, and `psi_out(Z, parity=1)` is
f ∘ Ψ_out(Z − ½), the outgoing coordinate of the other petal.

The Lavaurs map is L_σ(w) = Ψ_out,π(Φ_in(w) + σ), where π is the parity of w.
The horn map E = Φ_in ∘ Ψ_out satisfies E(Z + 1) = E(Z) + 1, and at the end
where the parabolic cycle sits (`end` = +1 upper, −1 lower),
E(Z) − Z = −σ∞ + Σ_{k≥1} a_k e^{2πi·end·kZ}, with σ∞ = end·c·iπ: the log
branches of the two petals differ by iπ there. σ∞ is the phase of the roots as q
grows: iπ for `ONE` (σ_q = q − π/ε, `lavaurs_phase`), and −11πi/16 for `HALF`.
Along λ = λ₀e^{u/q²} the phase moves by end·u/2πi, so the fixed point leaving
the end has multiplier 1 − u + κ₀u² + … with κ₀ = end·a₂/(2πi a₁²) (`kappa0`).
The ratio is invariant under translating either Fatou coordinate.

Arithmetic is mpmath at `DPS` digits.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

import mpmath as mp

TAGS = ("Lav89-Lavaurs-phase", "Shi00-parabolic-implosion", "EV-horn-map")

DPS = 40
ESCAPE = 10
STEPS = 400_000


@dataclass(frozen=True)
class Germ:
    name: str
    z0: Fraction
    step: Fraction
    petals: int  # attracting petals of f^petals; also the pole order of Φ
    log_power: int  # c·log(sign·w^m)
    log_sign_in: int
    log_sign_out: int
    incoming: complex  # direction of the incoming petal used by Φ_in
    outgoing: complex  # direction of the outgoing petal used by Ψ_out
    end: int  # +1: the parabolic cycle sits at the upper end of the horn map, −1: lower
    height: float  # a line Im Z = height inside that end, for Fourier sampling
    radius: str
    order: int

    def f(self, w):
        return w * (w + 2 * self.z0)


ONE = Germ("p = 1", Fraction(1, 2), Fraction(1), 1, 1, -1, 1, -1, 1, 1, 3.0, "1e-3", 12)
HALF = Germ("p = (q - 1)/2", Fraction(-1, 2), Fraction(1, 2), 2, 2, 1, -1, 1, 1j, -1, -2.6, "0.012", 26)


def germ_of(p: int, q: int) -> Germ:
    if p == 1:
        return ONE
    if 2 * p + 1 == q:
        return HALF
    raise ValueError(f"no germ for {p}/{q}")


@dataclass(frozen=True)
class Expansion:
    principal: tuple[Fraction, ...]  # a_{−P}, …, a_{−1}
    log: Fraction
    taylor: tuple[Fraction, ...]  # d_1, d_2, …


def _mul(a, b, n):
    return [sum(a[i] * b[k - i] for i in range(k + 1) if i < len(a) and k - i < len(b)) for k in range(n)]


def _power(u, e, n):
    """Coefficients 0 … n−1 of u^e for a series u with u[0] = ±1."""
    base = u
    if e < 0:
        inverse = [1 / u[0]]
        for k in range(1, n):
            inverse.append(-sum(u[i] * inverse[k - i] for i in range(1, min(k, len(u) - 1) + 1)) / u[0])
        base, e = inverse, -e
    out = [Fraction(1)] + [Fraction(0)] * (n - 1)
    for _ in range(e):
        out = _mul(out, base, n)
    return out


def _log(u, n):
    """Coefficients 0 … n−1 of log(u/u[0])."""
    x = [Fraction(0)] + [c / u[0] for c in u[1:]] + [Fraction(0)] * n
    out, term = [Fraction(0)] * n, [Fraction(1)] + [Fraction(0)] * (n - 1)
    for j in range(1, n):
        term = _mul(term, x, n)
        out = [o + Fraction((-1) ** (j + 1), j) * t for o, t in zip(out, term)]
    return out


def _solve(A, b):
    n = len(b)
    M = [row[:] + [bi] for row, bi in zip(A, b)]
    for col in range(n):
        pivot = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[pivot] = M[pivot], M[col]
        for r in range(n):
            if r != col and M[r][col] != 0:
                factor = M[r][col] / M[col][col]
                M[r] = [x - factor * y for x, y in zip(M[r], M[col])]
    return [M[i][n] / M[i][i] for i in range(n)]


@lru_cache(maxsize=None)
def fatou_coefficients(germ: Germ = ONE) -> Expansion:
    """The expansion of Φ, exact: Φ(f(w)) − Φ(w) = step, order by order.

    With u = f(w)/w, each basis term changes by w^e(u^e − 1), m·log(u/u(0)), or
    w^k(u^k − 1). The equations at orders 1 − P … K + 1 are square in the
    unknowns a_{−P} … a_{−1}, c, d_1 … d_K.
    """
    P, K = germ.petals, germ.order
    u = [2 * germ.z0, Fraction(1)]
    lo, hi = 1 - P, K + 1
    n = hi - lo + 1 + P
    columns = []
    for e in list(range(-P, 0)) + [None] + list(range(1, K + 1)):
        if e is None:
            series, shift = [germ.log_power * c for c in _log(u, n)], 0
        else:
            series, shift = [c - (i == 0) for i, c in enumerate(_power(u, e, n))], e
        columns.append([series[o - shift] if 0 <= o - shift < n else Fraction(0) for o in range(lo, hi + 1)])
    A = [[col[r] for col in columns] for r in range(hi - lo + 1)]
    b = [germ.step if o == 0 else Fraction(0) for o in range(lo, hi + 1)]
    x = _solve(A, b)
    return Expansion(tuple(x[:P]), x[P], tuple(x[P + 1:]))


def num(c: Fraction):
    return mp.mpf(c.numerator) / c.denominator


def _expansion(w, outgoing: bool, germ: Germ = ONE):
    e = fatou_coefficients(germ)
    P, m = germ.petals, germ.log_power
    sign = germ.log_sign_out if outgoing else germ.log_sign_in
    value = (sum(num(a) * w ** (j - P) for j, a in enumerate(e.principal)) + num(e.log) * mp.log(sign * w**m)
             + sum(num(d) * w ** (k + 1) for k, d in enumerate(e.taylor)))
    slope = (sum((j - P) * num(a) * w ** (j - P - 1) for j, a in enumerate(e.principal)) + num(e.log) * m / w
             + sum((k + 1) * num(d) * w**k for k, d in enumerate(e.taylor)))
    return value, slope


def _inside(w, direction, radius) -> bool:
    v = w / direction
    return abs(w) < radius and abs(v.imag) < v.real


def phi_in(w, germ: Germ = ONE):
    """(Φ_in(w), Φ_in'(w), parity) for w in the parabolic basin; ValueError if the orbit escapes.

    The parity is the number of steps to the incoming petal, mod `petals`.
    """
    with mp.workdps(DPS):
        radius = mp.mpf(germ.radius)
        w, dw, n = mp.mpc(w), mp.mpc(1), 0
        while not _inside(w, germ.incoming, radius):
            dw, w, n = dw * (2 * w + 2 * germ.z0), germ.f(w), n + 1
            if abs(w) > ESCAPE or n > STEPS:
                raise ValueError("the orbit escapes or lingers: not verified in the parabolic basin")
        value, slope = _expansion(w, False, germ)
        return value - n * num(germ.step), slope * dw, n % germ.petals


def psi_out(Z, germ: Germ = ONE, parity: int = 0):
    """(Ψ_out(Z), Ψ_out'(Z)): f^m of the outgoing-petal inverse at Z − m·step, m a multiple of `petals`.

    `parity` j gives f^j ∘ Ψ_out(Z − j·step), the coordinate of the j-th image petal.
    """
    with mp.workdps(DPS):
        radius, step, P = mp.mpf(germ.radius), num(germ.step), germ.petals
        lead = num(fatou_coefficients(germ).principal[0])
        Z = mp.mpc(Z) - parity * step
        rounds = int(mp.ceil((Z.real + abs(lead) / radius**P + 2 * abs(Z.imag)) / (P * step)))
        m = P * max(0, rounds)
        target = Z - m * step
        root = mp.root(lead / target, P)
        w = max((root * mp.expjpi(2 * mp.mpf(j) / P) for j in range(P)), key=lambda r: (r / germ.outgoing).real)
        for _ in range(100):
            value, slope = _expansion(w, True, germ)
            delta = (value - target) / slope
            w -= delta
            if abs(delta) < mp.mpf(10) ** (5 - DPS):
                break
        dw = 1 / _expansion(w, True, germ)[1]
        for _ in range(m + parity):
            dw, w = dw * (2 * w + 2 * germ.z0), germ.f(w)
        return w, dw


def horn(Z, germ: Germ = ONE):
    """(E(Z) − Z, E'(Z)) for the horn map E = Φ_in ∘ Ψ_out."""
    with mp.workdps(DPS):
        w, dw = psi_out(Z, germ)
        value, slope, _ = phi_in(w, germ)
        return value - Z, slope * dw


def phase(germ: Germ = ONE):
    """σ∞ = end·c·iπ: iπ for `ONE`, −11πi/16 for `HALF`."""
    with mp.workdps(DPS):
        return germ.end * num(fatou_coefficients(germ).log) * 1j * mp.pi


def lavaurs_phase(q: int):
    """σ_q = q − π/ε at the root of B_{1/q} (`ONE`): q − π cot(π/q) + iπ."""
    with mp.workdps(DPS):
        return q - mp.pi * mp.cot(mp.pi / q) + 1j * mp.pi


def phase_along(u, q: int, germ: Germ = ONE):
    """σ_q(u) along λ = λ₀e^{u/q²}, up to a constant in u.

    `ONE`: q − π/ε = q − 2πi/(λ − 1). `HALF`: q/2 − T with the transit time
    T = −πi/log λ² of f², log λ² = −2πi/q + 2u/q² (the residue of α).
    """
    with mp.workdps(max(DPS, mp.mp.dps)):
        if germ == ONE:
            return q - 2j * mp.pi / (mp.expjpi(mp.mpf(2) / q) * mp.exp(u / mp.mpf(q) ** 2) - 1)
        if germ == HALF:
            return mp.mpf(q) / 2 + 1j * mp.pi / (-2j * mp.pi / q + 2 * u / mp.mpf(q) ** 2)
        raise ValueError(f"no phase formula for germ {germ.name}")


def phase_curvature(q: int, germ: Germ = ONE):
    """−end·2πi·q·[u²]σ_q(u): the part of q(κ − κ₀) that comes from the phase alone.

    The fixed point leaving the end has multiplier 1 − end·2πi·s + …, s = σ − σ∞,
    so the u² term of σ adds −end·2πi·[u²]σ to κ. It is i/π for `HALF` at every q
    and 1/(2πi) + O(q⁻⁴) for `ONE`.
    """
    with mp.workdps(2 * DPS):
        h = mp.mpf(10) ** (-DPS // 2)
        second = (phase_along(h, q, germ) - 2 * phase_along(0, q, germ) + phase_along(-h, q, germ)) / (2 * h * h)
        return -germ.end * 2j * mp.pi * q * second


def horn_coefficients(kmax: int, height: float | None = None, nodes: int = 48, germ: Germ = ONE) -> list:
    """a_1 … a_kmax of E(Z) − Z + σ∞ = Σ a_k e^{2πi·end·kZ}, sampled on Im Z = height."""
    with mp.workdps(DPS):
        y = mp.mpf(germ.height if height is None else height)
        s = phase(germ)
        h = [horn(mp.mpc(mp.mpf(j) / nodes, y), germ)[0] + s for j in range(nodes)]
        return [sum(h[j] * mp.expjpi(-2 * germ.end * k * mp.mpf(j) / nodes) for j in range(nodes)) / nodes
                * mp.exp(2 * mp.pi * germ.end * k * y) for k in range(1, kmax + 1)]


def kappa0(height: float | None = None, germ: Germ = ONE):
    """end·a₂/(2πi a₁²)."""
    with mp.workdps(DPS):
        a1, a2 = horn_coefficients(2, height, germ=germ)
        return germ.end * a2 / (2j * mp.pi * a1**2)


def lavaurs_map(w, sigma=None, germ: Germ = ONE):
    """(L_σ(w), L_σ'(w)) for L_σ = Ψ_out,parity ∘ (· + σ) ∘ Φ_in; default σ = σ∞."""
    with mp.workdps(DPS):
        sigma = phase(germ) if sigma is None else mp.mpc(sigma)
        value, slope, parity = phi_in(w, germ)
        image, dimage = psi_out(value + sigma, germ, parity)
        return image, dimage * slope


def lavaurs_fixed_point(w, sigma=None, germ: Germ = ONE, tol: float = 1e-25, cap: float = 0.2):
    """(w*, L_σ'(w*)) with L_σ(w*) = w*, by Newton from w with steps capped at `cap`.

    L_σ'(w*) = E'(Z*) at Z* = Φ_in(w*) + σ, the matching fixed point of the horn map.
    """
    with mp.workdps(DPS):
        w = mp.mpc(w)
        for _ in range(100):
            image, slope = lavaurs_map(w, sigma, germ)
            delta = (image - w) / (slope - 1)
            w -= delta * min(1, cap / abs(delta))
            if abs(delta) < tol:
                return w, lavaurs_map(w, sigma, germ)[1]
        raise ValueError(f"Newton did not converge from {w}")


def fixed_point_of_cycle(z, sigma=None, seeds: int = 4, germ: Germ = ONE, cap: float = 0.2):
    """The fixed points of L_σ seeded by the `seeds` points of a cycle farthest from z₀.

    Returns (w*, L_σ'(w*)) for each seed that converges, ordered by the seed's
    residual |L_σ(w) − w|; for a cycle that L_σ models they coincide, and the first
    is the one to use when a poor seed reaches a neighbouring fixed point.
    """
    with mp.workdps(DPS):
        w = sorted((mp.mpc(complex(x)) - num(germ.z0) for x in z), key=abs, reverse=True)
        found = []
        for x in w[:seeds]:
            try:
                residual = abs(lavaurs_map(x, sigma, germ)[0] - x)
                found.append((residual, lavaurs_fixed_point(x, sigma, germ, cap=cap)))
            except ValueError:
                pass
        return [point for _, point in sorted(found, key=lambda f: f[0])]
