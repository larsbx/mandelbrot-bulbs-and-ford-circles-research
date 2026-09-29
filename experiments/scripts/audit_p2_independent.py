"""Independent re-derivation of P2 for the contributions audit.

Deliberately imports nothing from `bulbford`: the satellite multiplier ρ(ε) is
computed by deflated Newton on f^q(z) = z for f(z) = λz + z², λ = λ₀e^ε, and
the index ι by a contour integral of 1/(z − f^q(z)) on |z| = 10⁻³ at λ₀. A
cubic through three ε values gives [ε]ρ and [ε²]ρ, compared with −q² and
q³(ι − ½) (P2). Also prints the exact V13 values of ι_{1/3}, ι_{1/4}.

    python3 experiments/scripts/audit_p2_independent.py
"""
from __future__ import annotations

import mpmath as mp

mp.mp.dps = 60
FRACTIONS = ((1, 3), (1, 4), (2, 5))
EPS = tuple(mp.mpf(s) * mp.mpf("1e-6") for s in (1, -1, 2))
SEEDS = tuple(mp.mpf("0.02") * mp.expj(2 * mp.pi * k / 24) for k in range(24))


def orbit_jet(lam, q: int, z):
    """(f^q(z), (f^q)'(z)) for f(z) = λz + z²."""
    w, d = z, mp.mpc(1)
    for _ in range(q):
        d, w = d * (lam + 2 * w), lam * w + w * w
    return w, d


def satellite_multiplier(lam, q: int, z):
    """Multiplier of the cycle reached by Newton on (f^q(z) − z)/z, which deflates z = 0."""
    for _ in range(200):
        w, d = orbit_jet(lam, q, z)
        z -= (w - z) * z / ((d - 1) * z - (w - z))
    return orbit_jet(lam, q, z)[1]


def nearest_to_one(lam, q: int):
    return min((satellite_multiplier(lam, q, s) for s in SEEDS), key=lambda m: abs(m - 1))


def index(lam0, q: int):
    r = mp.mpf("1e-3")
    integrand = lambda t: 1j * r * mp.expj(t) / (r * mp.expj(t) - orbit_jet(lam0, q, r * mp.expj(t))[0])
    return mp.quad(integrand, [0, 2 * mp.pi]) / (2j * mp.pi)


def check(p: int, q: int) -> str:
    lam0 = mp.expj(2 * mp.pi * p / q)
    rho = [nearest_to_one(lam0 * mp.exp(e), q) for e in EPS]
    a, b, _ = mp.lu_solve(mp.matrix([[e, e**2, e**3] for e in EPS]), mp.matrix([r - 1 for r in rho]))
    iota = index(lam0, q)
    return (f"{p}/{q}: [ε]ρ = {mp.nstr(a, 10)} (P2: {-q * q});  [ε²]ρ = {mp.nstr(b, 8)}, "
            f"q³(ι − ½) = {mp.nstr(q**3 * (iota - mp.mpf(1) / 2), 8)};  ι = {mp.nstr(iota, 12)}")


if __name__ == "__main__":
    print(*map(lambda f: check(*f), FRACTIONS), sep="\n")
    print("V13 exact: ι_{1/3} =", mp.nstr((92 - 16 * mp.expj(2 * mp.pi / 3)) / 441, 12),
          " ι_{1/4} =", mp.nstr(mp.mpc(1447, -365) / 4624, 12))
