"""Finite joint box certificates for critical-orbit types, applied to bulb centres.

The calculus is the one of `larsbx/finite-mandelbrot-research`,
`docs/finite-certificate-calculus.md`: for a claimed critical-orbit type
(ℓ, k) and a dyadic complex box β,

    W_box(β):  K_R(β) ⊂ int β                                (Krawczyk localisation)
               0 ∉ Q_j(β) − Q_i(β)  for every (i, j) ∈ F_{ℓ,k}(H)  (same-box exclusions)

with Q₀ = 0, Q_{n+1} = Q_n² + C and R = Q_{ℓ+k} − Q_ℓ.  A satellite centre of
period q is type (0, q): R = Q_q is the Gleason polynomial, and the exclusions
say the critical orbit has exact period q.  Krawczyk inclusion already proves
the root in β is simple, so no separate squarefree step is needed.

All numbers are rationals with dyadic denominators; every operation rounds
outward to `prec` bits, so each box encloses the exact value.  The Krawczyk
operator and its preconditioner are the vendored `root_isolation_py`
(larsbx/finite-math-kernels, docs/root-isolation-spec.md), on `closed_interval`
boxes converted at the boundary; this module supplies the maps, the fixed-point
rounding to the 2^-prec grid, and the zero preconditioner for a singular centre.  Seeds may come
from anywhere (floating-point Newton, continuation): they are untrusted, and
only the inclusion checks decide.  An unmet check is INCONCLUSIVE, never a
disproof.  Which bulb a certified centre belongs to is not a finite fact here;
it is recorded as a named import (`TAGS`).
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import ceil, floor
from typing import Callable

from closed_interval import IQ, ComplexIQ
from root_isolation_py import krawczyk as _krawczyk
from root_isolation_py import midpoint_inverse, strictly_inside

PREC = 160  # bits after the binary point kept by outward rounding


# --- rational intervals (docs/rational-interval-arithmetic-spec.md semantics) ---------------


def _down(x: Fraction, prec: int) -> Fraction:
    return Fraction(floor(x * 2**prec), 2**prec)


def _up(x: Fraction, prec: int) -> Fraction:
    return Fraction(ceil(x * 2**prec), 2**prec)


@dataclass(frozen=True, slots=True)
class I:
    lo: Fraction
    hi: Fraction

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError("empty interval")

    @staticmethod
    def point(x) -> "I":
        return I(Fraction(x), Fraction(x))

    def rounded(self, prec: int) -> "I":
        return I(_down(self.lo, prec), _up(self.hi, prec))

    def __add__(self, o: "I") -> "I":
        return I(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o: "I") -> "I":
        return I(self.lo - o.hi, self.hi - o.lo)

    def __mul__(self, o: "I") -> "I":
        v = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return I(min(v), max(v))

    def contains_zero(self) -> bool:
        return self.lo <= 0 <= self.hi


@dataclass(frozen=True, slots=True)
class Box:
    """Rank-2 record (re, im) of intervals; multiplication is (xu − yv, xv + yu)."""

    re: I
    im: I

    @staticmethod
    def point(re, im=0) -> "Box":
        return Box(I.point(re), I.point(im))

    @staticmethod
    def around(re: Fraction, im: Fraction, radius: Fraction) -> "Box":
        return Box(I(re - radius, re + radius), I(im - radius, im + radius))

    def rounded(self, prec: int) -> "Box":
        return Box(self.re.rounded(prec), self.im.rounded(prec))

    def __add__(self, o: "Box") -> "Box":
        return Box(self.re + o.re, self.im + o.im)

    def __sub__(self, o: "Box") -> "Box":
        return Box(self.re - o.re, self.im - o.im)

    def __mul__(self, o: "Box") -> "Box":
        return Box(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    def excludes_zero(self) -> bool:
        return not (self.re.contains_zero() and self.im.contains_zero())

    def strictly_inside(self, o: "Box") -> bool:
        """The Krawczyk-Moore test of the vendored root_isolation_py."""
        return strictly_inside((_ciq(self),), (_ciq(o),))

    @property
    def mid(self) -> tuple[Fraction, Fraction]:
        return (self.re.lo + self.re.hi) / 2, (self.im.lo + self.im.hi) / 2


ZERO, ONE = Box.point(0), Box.point(1)


def _ciq(b: Box) -> ComplexIQ:
    return ComplexIQ(IQ(b.re.lo, b.re.hi), IQ(b.im.lo, b.im.hi))


def _box(z: ComplexIQ) -> Box:
    """Back from `closed_interval`; a rejected box has lo > hi, which `I` refuses."""
    return Box(I(z.re.lo, z.re.hi), I(z.im.lo, z.im.hi))


# --- critical-orbit polynomials on boxes -------------------------------------------------------


def orbit_jet(c: Box, n: int, prec: int = PREC) -> tuple[tuple[Box, Box], ...]:
    """((Q_0, Q_0'), …, (Q_n, Q_n')) on c: Q_{m+1} = Q_m² + C, Q_{m+1}' = 2 Q_m Q_m' + 1."""
    two = Box.point(2)
    step = lambda acc, _: acc + (
        ((acc[-1][0] * acc[-1][0] + c).rounded(prec), (two * acc[-1][0] * acc[-1][1] + ONE).rounded(prec)),
    )
    return reduce(step, range(n), ((ZERO, ZERO),))


def return_map(ell: int, period: int) -> Callable[[Box, int], tuple[Box, Box]]:
    """R = Q_{ℓ+k} − Q_ℓ and R' on a box."""

    def r(c: Box, prec: int = PREC) -> tuple[Box, Box]:
        jet = orbit_jet(c, ell + period, prec)
        return jet[ell + period][0] - jet[ell][0], jet[ell + period][1] - jet[ell][1]

    return r


# --- collision sets --------------------------------------------------------------------------


def intended_pairs(ell: int, period: int, horizon: int) -> frozenset[tuple[int, int]]:
    return frozenset(
        (i, j) for i, j in combinations(range(horizon + 1), 2) if i >= ell and (j - i) % period == 0
    )


def forbidden_pairs(ell: int, period: int, horizon: int) -> frozenset[tuple[int, int]]:
    return frozenset(combinations(range(horizon + 1), 2)) - intended_pairs(ell, period, horizon)


def failed_exclusions(beta: Box, ell: int, period: int, horizon: int, prec: int = PREC) -> tuple[tuple[int, int], ...]:
    """Forbidden pairs whose difference box on β still contains 0 (empty ⇔ all excluded)."""
    orbit = [z for z, _ in orbit_jet(beta, horizon, prec)]
    return tuple(sorted((i, j) for i, j in forbidden_pairs(ell, period, horizon) if not (orbit[j] - orbit[i]).excludes_zero()))


# --- Krawczyk ------------------------------------------------------------------------------


def _complex_dyadic(z: complex, prec: int) -> tuple[Fraction, Fraction]:
    return _down(Fraction(z.real), prec), _down(Fraction(z.imag), prec)


def _preconditioner(prec: int) -> Callable:
    """`root_isolation_py.midpoint_inverse` rounded down to the 2^-prec grid, or the zero matrix where the
    centre is singular (then K(β) ⊇ β and the inclusion fails: INCONCLUSIVE, not an exception).  Any point
    matrix is admissible (docs/root-isolation-spec.md, section 4); its accuracy affects only acceptance."""

    def invert(j):
        inverse = midpoint_inverse(j, lambda x: _down(x, prec))
        if all(e.accepted() for row in inverse for e in row):
            return inverse
        return tuple(tuple(_ciq(ZERO) for _ in row) for row in j)

    return invert


def _approx_inverse(d: Box, prec: int) -> Box:
    """A dyadic point near 1/d (ZERO when mid d = 0)."""
    return _box(_preconditioner(prec)(((_ciq(d),),))[0][0])


def krawczyk_boxes(f, x: tuple[Box, ...], prec: int = PREC) -> tuple[Box, ...]:
    """K(X) = m − A F(m) + (I − A J(X))(X − m), m = mid X, A ≈ J(m)⁻¹, rounded outward to `prec` bits:
    the vendored `root_isolation_py.krawczyk` with f(boxes, prec) -> (values, Jacobian) on bulbford boxes."""

    def at(v):
        values, jacobian = f(tuple(map(_box, v)), prec)
        return tuple(map(_ciq, values)), tuple(tuple(map(_ciq, row)) for row in jacobian)

    image = _krawczyk(tuple(map(_ciq, x)), at, lambda v: at(v)[1], _preconditioner(prec),
                      lambda z: _ciq(_box(z).rounded(prec)))
    return tuple(map(_box, image))


def krawczyk(r, beta: Box, prec: int = PREC) -> Box:
    """K(β) = m − A R(m) + (1 − A R'(β))(β − m), m = mid β, A ≈ 1/R'(m)."""

    def f(v, p):
        value, deriv = r(v[0], p)
        return (value,), ((deriv,),)

    return krawczyk_boxes(f, (beta,), prec)[0]


def newton_refine(r, re: Fraction, im: Fraction, steps: int, prec: int = PREC) -> tuple[Fraction, Fraction]:
    """Untrusted high-precision seed polishing on dyadic points."""

    def step(z, _):
        value, deriv = r(Box.point(*z), prec)
        corr = _approx_inverse(deriv, prec) * value
        return _down(z[0] - corr.re.lo, prec), _down(z[1] - corr.im.lo, prec)

    return reduce(step, range(steps), (re, im))


# --- joint certificate ------------------------------------------------------------------------


class Verdict(str, Enum):
    ACCEPTED = "ACCEPTED"
    INCONCLUSIVE = "INCONCLUSIVE"  # a check failed on this box: no conclusion, never a disproof


#: Named analytic imports.  None is used to *check* a certificate; they are what a
#: consumer must add to read an accepted box as a statement about a bulb.
TAGS = {
    "DH-multiplier": "[DH] the multiplier map of a hyperbolic component is a conformal isomorphism onto 𝔻;"
    " its unique zero is the component's centre",
    "SatelliteLabel": "the period-q centre found by Newton from c_of(λ₀(1 + q⁻²)), and the ρ = −1 point"
    " continued from it, belong to B_{p/q} (VALIDATED numerically; not a finite fact of this module)",
}


@dataclass(frozen=True, slots=True)
class JointCertificate:
    ell: int
    period: int
    horizon: int
    box: Box
    krawczyk_inclusion: bool
    failed_exclusions: tuple[tuple[int, int], ...]
    prec: int

    @property
    def verdict(self) -> Verdict:
        both = self.krawczyk_inclusion and not self.failed_exclusions
        return Verdict.ACCEPTED if both else Verdict.INCONCLUSIVE


def certify_type(beta: Box, ell: int, period: int, horizon: int | None = None, prec: int = PREC) -> JointCertificate:
    """Check W_box(β) for type (ℓ, k): both halves are evaluated on the same β."""
    horizon = ell + period if horizon is None else horizon
    if horizon < ell + period or period < 1 or ell < 0:
        raise ValueError("need ℓ ≥ 0, k ≥ 1 and H ≥ ℓ + k")
    inclusion = krawczyk(return_map(ell, period), beta, prec).strictly_inside(beta)
    return JointCertificate(ell, period, horizon, beta, inclusion, failed_exclusions(beta, ell, period, horizon, prec), prec)


@dataclass(frozen=True, slots=True)
class CenterCertificate:
    p: int
    q: int
    certificate: JointCertificate

    @property
    def verdict(self) -> Verdict:
        return self.certificate.verdict

    def as_record(self) -> dict:
        c = self.certificate
        return {
            "p": self.p,
            "q": self.q,
            "type": [c.ell, c.period],
            "horizon": c.horizon,
            "prec": c.prec,
            "box_numerators": box_numerators(c.box, c.prec),
            "krawczyk_inclusion": c.krawczyk_inclusion,
            "failed_exclusions": [list(x) for x in c.failed_exclusions],
            "verdict": c.verdict.value,
            "tags": sorted(TAGS),
        }


def box_numerators(box: Box, prec: int) -> list[str]:
    """(re_lo, re_hi, im_lo, im_hi) as integer numerators over 2^prec; refuses a non-dyadic endpoint."""
    ends = [e * 2**prec for e in (box.re.lo, box.re.hi, box.im.lo, box.im.hi)]
    if any(e.denominator != 1 for e in ends):
        raise ValueError("endpoint is not dyadic at this precision")
    return [str(e.numerator) for e in ends]


def box_from_numerators(nums: list[str], prec: int) -> Box:
    re_lo, re_hi, im_lo, im_hi = (Fraction(int(n), 2**prec) for n in nums)
    return Box(I(re_lo, re_hi), I(im_lo, im_hi))


def box_from_record(row: dict) -> Box:
    return box_from_numerators(row["box_numerators"], row["prec"])


def certify_center(p: int, q: int, seed: complex, radius_bits: int = 64, prec: int = PREC, newton_steps: int = 4) -> CenterCertificate:
    """Joint type-(0, q) certificate on a dyadic box of half-width 2^-radius_bits around a polished seed."""
    r = return_map(0, q)
    re, im = newton_refine(r, *_complex_dyadic(seed, prec), newton_steps, prec)
    beta = Box.around(re, im, Fraction(1, 2**radius_bits))
    return CenterCertificate(p, q, certify_type(beta, 0, q, q, prec))


def replay(row: dict) -> CenterCertificate:
    """Re-check a stored record from its box alone; the stored verdict is not trusted."""
    ell, period = row["type"]
    return CenterCertificate(row["p"], row["q"], certify_type(box_from_record(row), ell, period, row["horizon"], row["prec"]))
