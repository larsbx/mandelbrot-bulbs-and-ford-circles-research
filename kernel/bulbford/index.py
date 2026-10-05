"""Holomorphic index ι_{p/q} of the parabolic fixed point 0 of f(z) = λ₀ z + z², λ₀ = ζ_q^p.

ι = Res_{z=0} dz / (z − f^q(z)).  Since f^q(z) = z + a z^{q+1} + O(z^{q+2}) with a ≠ 0,
z − f^q(z) = z^{q+1} P(z), P(0) = −a, and ι = [z^q] (1/P).  Ball arithmetic (flint/arb);
the returned value carries a rigorous error radius.

ζ_q is built from its polynomial, not from an angle: Arb isolates the roots of Φ_q, and
ζ_q is the one in the open upper half plane with the largest real part, decided by
certain ball comparisons, the same order check as `antipode.zeta_box`.
"""
from __future__ import annotations
from functools import lru_cache
from flint import acb, acb_series, ctx, fmpz_poly


@lru_cache(maxsize=None)
def zeta_ball(q: int, prec: int) -> acb:
    """ζ_q as an Arb ball at `prec` bits, selected from the certified roots of Φ_q by order alone."""
    if q <= 2:
        return acb(1 if q == 1 else -1)
    old = ctx.prec
    try:
        ctx.prec = prec
        roots = [r for r, _ in fmpz_poly.cyclotomic(q).complex_roots()]
        upper = [r for r in roots if r.imag > 0]
        best = [r for r in upper if all(r is o or r.real > o.real for o in upper)]
    finally:
        ctx.prec = old
    if len(upper) != len(roots) // 2 or len(best) != 1:
        raise ValueError(f"ζ_{q} not separated by the order check at {prec} bits")
    return best[0]


def lambda_ball(p: int, q: int, prec: int) -> acb:
    """λ₀ = ζ_q^p as an Arb ball."""
    return zeta_ball(q, prec) ** (p % q)


def _series_fq(p: int, q: int, nterms: int) -> acb_series:
    with_prec = lambda coeffs: acb_series(coeffs, prec=nterms)
    lam0 = lambda_ball(p, q, ctx.prec)
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
