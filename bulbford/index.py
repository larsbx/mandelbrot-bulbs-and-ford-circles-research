"""Holomorphic index ι_{p/q} of the parabolic fixed point 0 of f(z) = λ₀ z + z², λ₀ = e^{2πip/q}.

ι = Res_{z=0} dz / (z − f^q(z)).  Since f^q(z) = z + a z^{q+1} + O(z^{q+2}) with a ≠ 0,
z − f^q(z) = z^{q+1} P(z), P(0) = −a, and ι = [z^q] (1/P).  Ball arithmetic (flint/arb);
the returned value carries a rigorous error radius.
"""
from __future__ import annotations
from functools import lru_cache
from flint import acb, acb_series, arb, ctx


def _series_fq(p: int, q: int, nterms: int) -> acb_series:
    with_prec = lambda coeffs: acb_series(coeffs, prec=nterms)
    lam0 = (acb(0, 2) * arb.pi() * p / q).exp()
    z = with_prec([0, 1])
    s = z
    for _ in range(q):
        s = lam0 * s + s * s
    return s


def _residue_series(p: int, q: int) -> acb:
    n = 2 * q + 2
    s = _series_fq(p, q, n)
    coeffs = (acb_series([0, 1], prec=n) - s).coeffs()
    coeffs = coeffs + [acb(0)] * (n - len(coeffs))
    P = acb_series(coeffs[q + 1:], prec=q + 1)        # (z − f^q)/z^{q+1}
    inv = P.inv().coeffs()
    return inv[q] if len(inv) > q else acb(0)


def index(p: int, q: int, dps: int | None = None) -> acb:
    """ι_{p/q} as an acb ball; precision raised automatically until the ball radius < 1e-12."""
    old, old_cap = ctx.dps, ctx.cap
    try:
        dps = dps or max(30, int(1.5 * q))
        ctx.cap = 2 * q + 2
        while True:
            ctx.dps = dps
            r = _residue_series(p, q)
            if r.rad() < 1e-12 or dps > 60 * q + 200:
                return r
            dps = 2 * dps
    finally:
        ctx.dps, ctx.cap = old, old_cap


@lru_cache(maxsize=None)
def index_complex(p: int, q: int) -> complex:
    return complex(index(p, q))


def kappa(p: int, q: int) -> complex:
    """κ(p/q) = (ι − ½)/q: the ε²-coefficient of ρ divided by q⁴ (theorem P2)."""
    return (index_complex(p, q) - 0.5) / q


def index_at(fam, p: int, q: int) -> complex:
    """Holomorphic index of the parabolic cycle point of the satellite root c_fam(λ₀), λ₀ = e^{2πip/q},
    for the unicritical family z ↦ z^d + c, at a point z₀ of the parabolic cycle of period per(H).
    ι = Res_{w=0} dw / (w − (F(z₀+w) − z₀)), F = f^{q·per}.  Ball arithmetic; z₀ refined by Newton in acb."""
    import cmath
    from flint import acb_poly
    n = q * fam.period
    old, oc = ctx.dps, ctx.cap
    try:
        dps = max(30, int(1.5 * n))
        while True:
            ctx.dps, ctx.cap = dps, 2 * n + 2
            lam0 = (acb(0, 2) * arb.pi() * p / q).exp()
            c = _c_of(fam, lam0)
            # cycle point: Newton on f^{per}(z) = z from the double-precision estimate
            z0 = acb(_cycle_point_estimate(fam, p, q))
            for _ in range(200):
                w, A = z0, acb(1)
                for _ in range(fam.period):
                    A = fam.d * w ** (fam.d - 1) * A; w = w ** fam.d + c
                z0 = z0 - (w - z0) / (A - 1)
            s = acb_series([z0, 1], prec=2 * n + 2)
            for _ in range(n):
                s = s ** fam.d + c
            coeffs = (acb_series([z0, 1], prec=2 * n + 2) - s).coeffs()
            coeffs = coeffs + [acb(0)] * (2 * n + 2 - len(coeffs))
            # parabolic with q petals at the cycle point: coefficients 1..q vanish; P = (w − F)/w^{q+1}
            P = acb_series(coeffs[q + 1:], prec=q + 1)
            inv = P.inv().coeffs()
            r = inv[q] if len(inv) > q else acb(0)
            if r.rad() < 1e-10 or dps > 60 * n + 200:
                return complex(r)
            dps = 2 * dps
    finally:
        ctx.dps, ctx.cap = old, oc


def _c_of(fam, lam0):
    """c_fam(λ₀) in ball arithmetic for the built-in families."""
    if fam.name == "main2": return lam0 / 2 - lam0 * lam0 / 4
    if fam.name == "disk2": return lam0 / 4 - 1
    if fam.name == "main3":
        z = (lam0 / 3).sqrt(); return z - z ** 3
    raise ValueError(fam.name)


def _cycle_point_estimate(fam, p: int, q: int) -> complex:
    import cmath
    lam0 = cmath.exp(2j * cmath.pi * p / q)
    if fam.name == "main2": return 0j
    if fam.name == "main3": return cmath.sqrt(lam0 / 3)
    if fam.name == "disk2":
        c = lam0 / 4 - 1; return (-1 + cmath.sqrt(-3 - 4 * c)) / 2
    raise ValueError(fam.name)
