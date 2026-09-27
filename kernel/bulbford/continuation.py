"""Certified continuation of the attracting cycle from a certified centre to a certified antipode.

For the segment c(s) = c₀ + s(c₁ − c₀), s ∈ [0, 1], from the midpoint c₀ of the P3 centre box to the
midpoint c₁ of the P10 antipode box, the path is covered by parameter boxes C_k with cycle boxes Z_k:

  1. parametric Krawczyk: K(Z_k; C_k) ⊂ int Z_k, so for every c ∈ C_k the map f_c^q has exactly one
     fixed point z_k(c) in Z_k, and z_k is analytic on C_k (the inclusion makes (f^q)′ − 1 invertible);
  2. |ρ| < 1 on Z_k × C_k, ρ = (f^q)′(z): the cycle attracts on all of C_k (interior pieces);
  3. junctions: at the shared exact point c(s_k), K(Z_{k+1}; {c(s_k)}) ⊂ Z_k or K(Z_k; {c(s_k)}) ⊂ Z_{k+1},
     so the two fixed points coincide there (each box holds only one), and
     then on the whole convex overlap, since both solve the same nonsingular equation;
  4. ends: C_0 ⊇ the P3 centre box; the last box C_L ⊇ the P10 antipode box and Z_L ⊇ its cycle box,
     and on Z_L × C_L the derivative of |ρ|² along v = c₁ − c₀ is certified positive:
     Re(conj(ρ)·(dρ/dc)·v) > 0, dρ/dc = ∂_zz f^q·dz/dc + ∂_zc f^q, dz/dc = −∂_c f^q/(∂_z f^q − 1).

Then every point of the half-open path from the true centre to the true antipode c_ant carries the same
attracting cycle: on the last piece |ρ| rises strictly to 1 along −v from c_ant, reaching back into the
previous piece. So c_ant lies on the boundary of the hyperbolic component of the certified centre, where
ρ = −1. This is finite; no theorem about which bulb that component is enters here.

Arithmetic: Arb complex balls (python-flint `acb`), rigorous and outward-rounded like the rational boxes
of bulbford.certify but in C. Boxes coming from the rational certificates enter exactly.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from flint import acb, arb, ctx

from .certify import Box

WORK_PREC = 192
MARGIN = Fraction(1, 2**60)  # ≫ the 2^-64 certificate boxes, so C_0 and C_L contain them


def _arb_exact(x: Fraction) -> arb:
    return arb(x.numerator) / arb(x.denominator)


def _ball(lo: Fraction, hi: Fraction) -> arb:
    return arb(_arb_exact((lo + hi) / 2), _arb_exact((hi - lo) / 2))


def _acb(b: Box) -> acb:
    return acb(_ball(b.re.lo, b.re.hi), _ball(b.im.lo, b.im.hi))


def _point(x: tuple[Fraction, Fraction]) -> acb:
    return acb(_arb_exact(x[0]), _arb_exact(x[1]))


def _widen(x: acb, re_rad: arb, im_rad: arb) -> acb:
    return acb(arb(x.real.mid(), x.real.rad() + re_rad), arb(x.imag.mid(), x.imag.rad() + im_rad))


def _inside(k: acb, z: acb) -> bool:
    """k ⊂ int z, decided on exact endpoints."""
    return all(bool(kp.lower() > zp.lower()) and bool(kp.upper() < zp.upper()) for kp, zp in ((k.real, z.real), (k.imag, z.imag)))


def _within(k: acb, z: acb) -> bool:
    return all(bool(kp.lower() >= zp.lower()) and bool(kp.upper() <= zp.upper()) for kp, zp in ((k.real, z.real), (k.imag, z.imag)))


def jet(z: acb, c: acb, q: int) -> tuple[acb, acb, acb, acb, acb]:
    """(f^q, ∂_z, ∂_c, ∂_zz, ∂_zc) of f_c(z) = z² + c at the balls (z, c)."""
    w, a, b, zz, zc = z, acb(1), acb(0), acb(0), acb(0)
    for _ in range(q):
        w, a, b, zz, zc = w * w + c, 2 * w * a, 2 * w * b + 1, 2 * (a * a + w * zz), 2 * (a * b + w * zc)
    return w, a, b, zz, zc


def param_krawczyk(q: int, z: acb, c: acb) -> acb:
    """K(Z; C) = m − A(f_C^q(m) − m) + (1 − A((f^q)′(Z; C) − 1))(Z − m), A ≈ 1/((f^q)′(m; mid C) − 1)."""
    m = acb(z.mid())
    w, *_ = jet(m, c, q)
    _, am, *_ = jet(m, acb(c.mid()), q)
    _, a_z, *_ = jet(z, c, q)
    inv = acb(1 / (am - 1).mid()).mid()
    return m - inv * (w - m) + (1 - inv * (a_z - 1)) * (z - m)


def _float_cycle_point(q: int, c: complex, z: complex) -> complex:
    """Untrusted Newton for f_c^q(z) = z (only chooses where to put a box)."""
    for _ in range(60):
        w, a = z, 1 + 0j
        for _ in range(q):
            a, w = 2 * w * a, w * w + c
        step = (w - z) / (a - 1)
        z -= step
        if abs(step) < 1e-15 * (1 + abs(z)):
            break
    return z


@dataclass(frozen=True, slots=True)
class Continuation:
    p: int
    q: int
    pieces: int
    last_s0: Fraction | None
    accepted: bool
    reason: str


def continue_centre_to_antipode(p: int, q: int, centre: Box, ant_c: Box, ant_z: Box,
                                min_ds: Fraction = Fraction(1, 2**26)) -> Continuation:
    old = ctx.prec
    ctx.prec = WORK_PREC
    try:
        return _continue(p, q, centre, ant_c, ant_z, min_ds)
    finally:
        ctx.prec = old


def _continue(p: int, q: int, centre: Box, ant_c: Box, ant_z: Box, min_ds: Fraction) -> Continuation:
    c0, c1 = centre.mid, ant_c.mid
    v = (c1[0] - c0[0], c1[1] - c0[1])
    vmax = max(abs(v[0]), abs(v[1]))
    at = lambda s: (c0[0] + s * v[0], c0[1] + s * v[1])
    fc = lambda s: complex(float(at(s)[0]), float(at(s)[1]))
    centre_b, ant_cb, ant_zb, v_b = _acb(centre), _acb(ant_c), _acb(ant_z), _point(v)

    def cbox(s0: Fraction, s1: Fraction) -> acb:
        a, b = at(s0), at(s1)
        m = max(MARGIN, abs(s1 - s0) * vmax / 20)
        return acb(_ball(min(a[0], b[0]) - m, max(a[0], b[0]) + m), _ball(min(a[1], b[1]) - m, max(a[1], b[1]) + m))

    def inflate(z: acb, c: acb, rounds: int = 8) -> acb | None:
        """ε-inflation: Z ← K(Z) widened by 1/4 (+ MARGIN) until K(Z) ⊂ int Z."""
        for _ in range(rounds):
            k = param_krawczyk(q, z, c)
            if not (k.real.is_finite() and k.imag.is_finite()):
                return None
            if _inside(k, z):
                return z
            z = _widen(k.union(z), k.real.rad() / 4 + _arb_exact(MARGIN), k.imag.rad() / 4 + _arb_exact(MARGIN))
        return None

    def attracting(z: acb, c: acb) -> bool:
        _, a, *_ = jet(z, c, q)
        return bool(abs(a).upper() < 1)  # abs, not ** 2: arb ** 2 is NaN on a ball containing 0

    def joins(z_prev: acb, z_next: acb, s: Fraction) -> bool:
        """At c(s) the fixed point of one box lies in the other, so by uniqueness there they coincide."""
        cs = _point(at(s))
        return _within(param_krawczyk(q, z_next, cs), z_prev) or _within(param_krawczyk(q, z_prev, cs), z_next)

    def monotone(z: acb, c: acb) -> bool:
        w, a, b, zz, zc = jet(z, c, q)
        drho = zz * (-b / (a - 1)) + zc
        g = (a.conjugate() * drho * v_b).real
        return g.is_finite() and bool(g.lower() > 0)

    def finish(s: Fraction, zf: complex, z_prev: acb) -> bool:
        cl = cbox(s, Fraction(1))
        if not cl.contains(ant_cb):
            return False
        zc = _float_cycle_point(q, fc((s + 1) / 2), zf)
        zl = inflate(_widen(acb(zc).union(ant_zb), _arb_exact(MARGIN), _arb_exact(MARGIN)), cl)
        return zl is not None and zl.contains(ant_zb) and monotone(zl, cl) and joins(z_prev, zl, s)

    pieces, s, ds, zf, z_prev = 0, Fraction(0), Fraction(1, 64), 0j, None
    while True:
        s1 = min(s + ds, Fraction(1))
        c_k = cbox(s, s1)
        zc = _float_cycle_point(q, fc((s + s1) / 2), zf)
        z_k = None
        if s1 < 1 and (pieces or c_k.contains(centre_b)):
            z_k = inflate(_widen(acb(zc), _arb_exact(MARGIN), _arb_exact(MARGIN)), c_k)
        if z_k is not None and attracting(z_k, c_k) and (z_prev is None or joins(z_prev, z_k, s)):
            pieces, s, zf, z_prev, ds = pieces + 1, s1, zc, z_k, min(ds * 2, Fraction(1, 8))
            continue
        # stuck (or at the end): try to close the path at the antipode from here
        if z_prev is not None and finish(s, zf, z_prev):
            return Continuation(p, q, pieces + 1, s, True, "accepted")
        ds /= 2
        if ds < min_ds:
            return Continuation(p, q, pieces, None, False, f"step below {min_ds} at s = {float(s):.6f}")
