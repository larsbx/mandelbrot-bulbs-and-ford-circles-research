"""Certified antipodes, roots of unity, and enclosures of G_ant.

Three finite certificates, all in the rational dyadic boxes of `certify.py`:

1. **ζ_q.**  Krawczyk boxes for all q roots of X^q − 1, pairwise disjoint, hence
   one per root.  ζ_q := e^{2πi/q} is selected by an order check, not an angle:
   it is the root in the open upper half plane with the largest real part.
   λ₀ = ζ_q^p is then an interval power, c_root = λ₀/2 − λ₀²/4 an interval
   polynomial.

2. **Antipode (type ρ = −1).**  For the system in (z, c)

       F₁ = f_c^q(z) − z,   F₂ = (f_c^q)'(z) + 1,

   a two-variable Krawczyk inclusion K(X) ⊂ int X gives a unique solution in
   X = Z × C; the same-box exclusions 0 ∉ f_c^j(Z) − f_c^i(Z) for the forbidden
   pairs of type (0, q) give z exact period q.  So some c ∈ C has a q-cycle of
   multiplier exactly −1.

3. **G.**  With c'(λ) = (1 − λ)/2 on the main cardioid,
   G_ant = q²|c_ant − c_root| / |1 − λ₀|.  The enclosure is computed for the
   quadrance G² = q⁴ Qd(c_ant − c_root) / Qd(1 − λ₀) (rank-2, no square root);
   G is bracketed from it by rational square-root bounds.

Imports needed to read (2)–(3) as bulb geometry: [DH] (the multiplier map of
B_{p/q} is a conformal isomorphism onto 𝔻, so ∂B_{p/q} has exactly one point
with ρ = −1) and `SatelliteLabel` (the certified cycle is the one of B_{p/q};
VALIDATED by continuation only).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache, reduce
from math import isqrt

from root_isolation_py import disjoint

from .certify import (
    ONE,
    PREC,
    ZERO,
    Box,
    I,
    TAGS,
    Verdict,
    _box,
    _ciq,
    _complex_dyadic,
    _down,
    _preconditioner,
    box_from_numerators,
    box_numerators,
    forbidden_pairs,
    krawczyk,
    krawczyk_boxes,
    newton_refine,
)
from .spread import turn_cosine_bracket

TWO = Box.point(2)


# --- boxes: powers, quadrance, division, square-root bounds ------------------------------------


def power(b: Box, n: int, prec: int = PREC) -> Box:
    """b^n by square-and-multiply: O(log n) box products, so the rectangular wrapping of a rotating
    factor compounds O(log n) times rather than n times (n = 251 overflows a 2^-100 box otherwise)."""
    if n == 0:
        return ONE
    half = power(b, n // 2, prec)
    sq = (half * half).rounded(prec)
    return (sq * b).rounded(prec) if n % 2 else sq


def quadrance(b: Box) -> I:
    """Qd(x, y) = x² + y² on a box, with the square of an interval taken tightly."""
    sq = lambda i: I(Fraction(0) if i.contains_zero() else min(i.lo**2, i.hi**2), max(i.lo**2, i.hi**2))
    return sq(b.re) + sq(b.im)


def divide_positive(a: I, b: I) -> I:
    if b.lo <= 0:
        raise ValueError("divisor interval must be positive")
    v = (a.lo / b.lo, a.lo / b.hi, a.hi / b.lo, a.hi / b.hi)
    return I(min(v), max(v))


def sqrt_bounds(i: I, prec: int = PREC) -> I:
    """[√lo, √hi] enlarged to dyadics of `prec` bits by integer square roots."""
    if i.lo < 0:
        raise ValueError("square root of a negative interval")
    scale = 4**prec
    lo = isqrt(int(i.lo * scale))
    hi_sq = -(-i.hi.numerator * scale // i.hi.denominator)
    hi = isqrt(hi_sq) + (isqrt(hi_sq) ** 2 != hi_sq)
    return I(Fraction(lo, 2**prec), Fraction(hi, 2**prec))


# --- 1: roots of unity -------------------------------------------------------------------------


def _unity(q: int):
    def r(x: Box, prec: int = PREC) -> tuple[Box, Box]:
        return power(x, q, prec) - ONE, (Box.point(q) * power(x, q - 1, prec)).rounded(prec)

    return r


def _disjoint(a: Box, b: Box) -> bool:
    return disjoint((_ciq(a),), (_ciq(b),))


def _unity_seed(k: int, q: int) -> complex:
    """An untrusted float seed for the k-th root of X^q − 1, with no angle: cos(π·2k/q) from the
    exact Sturm bracket of `spread.py`, and the sine as ±√(1 − cos²), positive for 0 < 2k mod 2q < q."""
    lo, hi = turn_cosine_bracket(2 * k, q, 60)
    c = float((lo + hi) / 2)
    s = math.sqrt(max(0.0, 1.0 - c * c))
    return complex(c, s if 0 < (2 * k) % (2 * q) < q else -s)


def unity_boxes(q: int, radius_bits: int = 100, prec: int = PREC) -> tuple[Box, ...]:
    """Krawczyk-certified, pairwise disjoint boxes for the q roots of X^q − 1 (seeds untrusted)."""
    boxes = []
    for k in range(q):
        re, im = newton_refine(_unity(q), *_complex_dyadic(_unity_seed(k, q), prec), 4, prec)
        beta = Box.around(re, im, Fraction(1, 2**radius_bits))
        if not krawczyk(_unity(q), beta, prec).strictly_inside(beta):
            raise ValueError(f"root-of-unity box {k}/{q} not certified")
        boxes.append(beta)
    if not all(_disjoint(a, b) for i, a in enumerate(boxes) for b in boxes[i + 1 :]):
        raise ValueError("root-of-unity boxes overlap")
    return tuple(boxes)


@lru_cache(maxsize=None)
def zeta_box(q: int, radius_bits: int = 100, prec: int = PREC) -> Box:
    """ζ_q: the upper-half-plane root with the largest real part, chosen by exact box comparisons."""
    if q <= 2:
        return Box.point(1 if q == 1 else -1)
    upper = [b for b in unity_boxes(q, radius_bits, prec) if b.im.lo > 0]
    best = [b for b in upper if all(b is o or o.re.hi < b.re.lo for o in upper)]
    if len(best) != 1:
        raise ValueError("ζ_q not separated by the order check")
    return best[0]


def lambda_box(p: int, q: int, prec: int = PREC) -> Box:
    return power(zeta_box(q, prec=prec), p, prec)


def root_box(lam: Box, prec: int = PREC) -> Box:
    """c_root = λ/2 − λ²/4."""
    return (Box.point(Fraction(1, 2)) * lam - Box.point(Fraction(1, 4)) * lam * lam).rounded(prec)


# --- 2: antipode --------------------------------------------------------------------------------


BLOWUP = Fraction(2**40)


class Blowup(ArithmeticError):
    """An orbit enclosure left |x| ≤ 2^40: the box is too wide for this orbit, so no inclusion can follow."""


def _bounded(b: Box) -> Box:
    if max(abs(b.re.lo), abs(b.re.hi), abs(b.im.lo), abs(b.im.hi)) > BLOWUP:
        raise Blowup
    return b


def expansion_bits(q: int, z: complex, c: complex) -> int:
    """⌈log₂ max_k ∏_{i<k} (|Re 2z_i| + |Im 2z_i|)⌉ along the (untrusted, floating-point) orbit.

    Rectangular arithmetic widens a box multiplied by a by |Re a| + |Im a| (not |a|): the rotation is
    wrapped at every step, so this, not |(f^k)'|, is how far an enclosure grows along the orbit."""
    log2, best = 0.0, 0.0
    for _ in range(q):
        log2 += math.log2(max(abs((2 * z).real) + abs((2 * z).imag), 1e-300))
        best, z = max(best, log2), z * z + c
    return math.ceil(best)


def jet(z: Box, c: Box, n: int, prec: int = PREC) -> tuple[Box, Box, Box, Box, Box]:
    """(f^n, ∂_z f^n, ∂_c f^n, ∂_zz f^n, ∂_zc f^n) at (z, c); raises Blowup if the orbit enclosure explodes."""

    def step(s, _):
        w, a, b, zz, zc = s
        r = lambda x: x.rounded(prec)
        return (_bounded(r(w * w + c)), r(TWO * w * a), r(TWO * w * b + ONE), r(TWO * (a * a + w * zz)), r(TWO * (a * b + w * zc)))

    return reduce(step, range(n), (z, ONE, ZERO, ZERO, ZERO))


Vec = tuple[Box, Box]
Mat = tuple[tuple[Box, Box], tuple[Box, Box]]


def system(q: int):
    """(F, J) for F = (f^q(z) − z, (f^q)'(z) + 1) in the unknowns (z, c)."""

    def f(x: Vec, prec: int = PREC) -> tuple[Vec, Mat]:
        w, a, b, zz, zc = jet(x[0], x[1], q, prec)
        return (w - x[0], a + ONE), ((a - ONE, b), (zz, zc))

    return f


def _matvec(m: Mat, v: Vec) -> Vec:
    return tuple(m[i][0] * v[0] + m[i][1] * v[1] for i in range(2))


def _dyadic_point(b: Box, prec: int) -> Box:
    re, im = b.mid
    return Box.point(_down(re, prec), _down(im, prec))


def _inverse_point(m: Mat, prec: int) -> Mat:
    """Dyadic point matrix near the inverse of mid(m) (zero where mid(m) is singular); affects acceptance
    only, never soundness."""
    inverse = _preconditioner(prec)(tuple(tuple(map(_ciq, row)) for row in m))
    return tuple(tuple(map(_box, row)) for row in inverse)


def krawczyk2(f, x: Vec, prec: int = PREC) -> Vec:
    """K(X) = m − A F(m) + (I − A J(X))(X − m)."""
    return krawczyk_boxes(f, x, prec)


def failed_cycle_exclusions(z: Box, c: Box, q: int, prec: int = PREC) -> tuple[tuple[int, int], ...]:
    """Forbidden pairs of type (0, q) on the orbit of Z under f_C (empty ⇔ exact period q)."""
    orbit = reduce(lambda acc, _: acc + (_bounded((acc[-1] * acc[-1] + c).rounded(prec)),), range(q), (z,))
    return tuple(sorted(p for p in forbidden_pairs(0, q, q) if not (orbit[p[1]] - orbit[p[0]]).excludes_zero()))


@dataclass(frozen=True, slots=True)
class AntipodeCertificate:
    p: int
    q: int
    z: Box
    c: Box
    krawczyk_inclusion: bool
    failed_exclusions: tuple[tuple[int, int], ...]
    prec: int

    @property
    def verdict(self) -> Verdict:
        both = self.krawczyk_inclusion and not self.failed_exclusions
        return Verdict.ACCEPTED if both else Verdict.INCONCLUSIVE

    def g_squared(self) -> I:
        """q⁴ Qd(C − c_root) / Qd(1 − λ₀), enclosing G_ant² for every point of C."""
        lam = lambda_box(self.p, self.q, self.prec)
        num = quadrance(self.c - root_box(lam, self.prec)) * I.point(self.q**4)
        return divide_positive(num, quadrance(ONE - lam))

    def g(self) -> I:
        return sqrt_bounds(self.g_squared(), self.prec)

    def as_record(self) -> dict:
        g = self.g()
        return {
            "p": self.p,
            "q": self.q,
            "prec": self.prec,
            "z_box_numerators": box_numerators(self.z, self.prec),
            "c_box_numerators": box_numerators(self.c, self.prec),
            "krawczyk_inclusion": self.krawczyk_inclusion,
            "failed_exclusions": [list(x) for x in self.failed_exclusions],
            "verdict": self.verdict.value,
            "G_ant_bounds": [str(g.lo), str(g.hi)],
            "tags": sorted(TAGS),
        }


BLOWN = ((-1, -1),)  # failed_exclusions marker: the orbit enclosure exploded before the horizon


# --- |u_a|/2 from the certified antipode (Arb balls) -----------------------------------------------


def _arb(i: I):
    """An Arb ball equal to the dyadic interval i (exact at the working precision set by the caller)."""
    from flint import arb

    to = lambda x: arb(x.numerator) / arb(x.denominator)
    return arb(to((i.lo + i.hi) / 2), to((i.hi - i.lo) / 2))


def _fraction(x) -> Fraction:
    m, e = x.man_exp()
    return Fraction(int(m)) * Fraction(2) ** int(e)


def u_half_abs(p: int, q: int, c: Box, prec: int) -> I:
    """Enclosure of |u_a|/2 for every c_ant ∈ C, where λ₀e^{u_a/q²} is the main-cardioid parameter of c_ant.

    From c = λ/2 − λ²/4: 1 − 4c = (1 − λ)², and the branch near λ₀ has Re(1 − λ) > 0, i.e. λ = 1 − √(1 − 4c)
    with the principal root, which is checked on the ball (refused otherwise). Then u = q²·log(λ/λ₀), with λ₀
    the angle-free ball `index.lambda_ball` (a root of Φ_q selected by order). Arb's √ and log are rigorous,
    so this is a certified enclosure; reading u_a as the root of
    R_q(u) = −1 near u = 2 uses [DH] (one ρ = −1 point on ∂B_{p/q}), as P10 does for c_ant."""
    from flint import acb, ctx

    from .index import lambda_ball

    old = ctx.prec
    try:
        ctx.prec = prec + 64
        cb = acb(_arb(c.re), _arb(c.im))
        root = (1 - 4 * cb).sqrt()
        if not root.real > 0:
            raise ValueError("branch not separated: Re √(1 − 4c) does not exclude 0")
        lam = 1 - root
        lam0 = lambda_ball(p, q, ctx.prec)
        half = abs(q * q * (lam / lam0).log()) / 2
        return I(_fraction(half.lower()), _fraction(half.upper()))
    finally:
        ctx.prec = old


def check_antipode(p: int, q: int, z: Box, c: Box, prec: int = PREC) -> AntipodeCertificate:
    try:
        inside = all(k.strictly_inside(b) for k, b in zip(krawczyk2(system(q), (z, c), prec), (z, c)))
    except Blowup:
        inside = False
    try:
        failed = failed_cycle_exclusions(z, c, q, prec)
    except Blowup:
        failed = BLOWN
    return AntipodeCertificate(p, q, z, c, inside, failed, prec)


def certify_antipode(p: int, q: int, seed_z: complex, seed_c: complex, radius_bits: int | None = 64, prec: int | None = PREC, newton_steps: int = 4) -> AntipodeCertificate:
    """Polish (z, c) by untrusted dyadic Newton steps, then check on a box of half-width 2^-radius_bits.

    radius_bits=None sizes the box from the orbit: expansion_bits + 48 (at least 64). prec=None keeps
    expansion_bits + 96 bits below the box, because a point evaluation of the jet also loses the orbit's
    growth, so the Newton residual floor sits near 2^(expansion − prec). Both only choose the box;
    acceptance is still the inclusion checks."""
    growth = expansion_bits(q, seed_z, seed_c) if radius_bits is None or prec is None else 0
    radius_bits = max(64, growth + 48) if radius_bits is None else radius_bits
    prec = radius_bits + growth + 96 if prec is None else prec
    f = system(q)

    def step(x: Vec, _) -> Vec:
        value, jac = f(x, prec)
        d = _matvec(_inverse_point(jac, prec), value)
        return tuple(_dyadic_point(x[i] - d[i], prec) for i in range(2))

    x0 = tuple(Box.point(*_complex_dyadic(s, prec)) for s in (seed_z, seed_c))
    mz, mc = (b.mid for b in reduce(step, range(newton_steps), x0))
    r = Fraction(1, 2**radius_bits)
    return check_antipode(p, q, Box.around(*mz, r), Box.around(*mc, r), prec)


def replay(row: dict) -> AntipodeCertificate:
    """Re-check a stored record from its boxes alone; stored verdict and bounds are not trusted."""
    prec = row["prec"]
    return check_antipode(row["p"], row["q"], box_from_numerators(row["z_box_numerators"], prec), box_from_numerators(row["c_box_numerators"], prec), prec)
