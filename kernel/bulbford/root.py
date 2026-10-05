"""The parabolic root side of the satellite label: the certified centre's component has c_root on its boundary.

Fix p/q, the root c_root = λ₀/2 − λ₀²/4 (λ₀ = e^{2πip/q}), and a circle ∂D (centre m near z₀ = λ₀/2, radius r)
with a smaller concentric disk D′ (radius r′). For parameter boxes K covering the segment from c_root to a
point c₁ on the chord towards the certified centre, with F_c = f_c^q − id and certified contour means:

  1. (1/2πi)∮_{∂D} F′/F = q + 1 = (1/2πi)∮_{∂D′} F′/F for every c ∈ K, and for each proper divisor d of q
     the count for f^d − id on ∂D is 1. With f(D̄′) ⊂ D (a ball check), the q zeros of F in D′ other than
     α(c) (the fixed point of f near z₀) are permuted by f and have exact period q: one q-cycle C(c).
  2. ρ_C·λ = λ₀^{q+1}·exp(S), S(c) = (1/2πi)∮_{∂D′} log(2z/λ₀)·F′/F dz (product of 2z over the q + 1 zeros;
     |2m/λ₀ − 1| + 2r′ < 1 is checked, so the log is holomorphic on all of D̄′), where
     λ = 2α = 1 − √(1 − 4c) on the branch Re √(1 − 4c) > 0 (checked). S is analytic in c across c_root, where
     C collides with α, and so is ρ_C.
  3. Along the segment c(t) = c_root + t·v: Re(conj(ρ_C)·ρ_C′·v) < 0 and Re(conj(λ)·λ′·v) > 0 on every box,
     so |ρ_C| decreases strictly from 1 at t = 0 and |λ| increases strictly from 1: for t > 0, C is an
     attracting q-cycle and α is not a multiple zero, so the half-open segment lies in one hyperbolic
     component H₁ with c_root ∈ ∂H₁.

The certified continuation from the P3 centre to c₁ (bulbford.continuation) then ends with its cycle box
inside D′, i.e. on C(c₁), so the centre lies in H₁. What remains imported is [DH]: the only hyperbolic
components whose closure contains the satellite root c_root(p/q) are the main cardioid and B_{p/q}, and a
component has one centre.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

from flint import acb, acb_poly, arb, ctx

from .continuation import _arb_exact
from .index import lambda_ball, zeta_ball

WORK_PREC = 128


@lru_cache(maxsize=None)
def _units(n: int, prec: int) -> tuple[acb, ...]:
    """All n-th roots of unity, from the certified polynomial/order interface."""
    old = ctx.prec
    try:
        ctx.prec = prec
        unit = zeta_ball(n, prec)
        return tuple(unit ** k for k in range(n))
    finally:
        ctx.prec = old


def _annulus(r: arb, a: arb, n: int) -> tuple[acb, ...]:
    """Boxes covering r/a ≤ |u| ≤ r*a, with no evaluation of angles.

    In a sector centred at an n-th root, its half-sector endpoint is ζ_{2n}.
    Its real coordinate is between Re(ζ_{2n})/a and a, and its imaginary
    coordinate is bounded by a*Im(ζ_{2n}). Rotating these boxes covers the annulus.
    """
    half = zeta_ball(2 * n, ctx.prec)
    sector = acb((half.real.lower() / a).union(a), arb(0, (a * half.imag).upper()))
    return tuple(r * unit * sector for unit in _units(n, ctx.prec))


def _circle(r: arb, density) -> acb:
    """Certified Laurent constant of H(u)=u*h(m+u), hence (1/2πi)∮h(z)dz.

    Each density here is rational, or rational times the principal logarithm.
    Finite ball evaluation on every box covering an annulus certifies that H
    has no pole or log cut there. If |H| ≤ M on r/a ≤ |u| ≤ r*a, Cauchy's
    coefficient bound gives |h_k*r^k| ≤ M*a^{-|k|}. The mean at the n-th
    roots of unity differs from h_0 by at most 2M/(a^n−1), since only the
    nonzero multiples of n survive that mean. All nodes and error bounds
    are algebraic balls. Failure to certify the annulus or error is NaN.
    """
    if not bool(r > 0):
        return acb("nan")
    for denominator in (16, 32, 64, 128):
        a = arb(denominator + 1) / denominator
        for sectors in (64, 128, 256, 512, 1024):
            bound = arb(0)
            for u in _annulus(r, a, sectors):
                value = density(u)
                if not value.is_finite():
                    break
                bound = max(bound, abs(value).upper())
            else:
                nodes = 64
                error = 2 * bound / (a ** nodes - 1)
                while not bool(error < arb(2) ** -40) and nodes < 8192:
                    nodes *= 2
                    error = 2 * bound / (a ** nodes - 1)
                if bool(error < arb(2) ** -40):
                    mean = sum((density(r * unit) for unit in _units(nodes, ctx.prec)), acb(0)) / nodes
                    return mean + acb(arb(0, error.upper()), arb(0, error.upper()))
    return acb("nan")


def _polynomials(n: int, c: acb, m: acb) -> tuple[acb_poly, ...]:
    """F, F_z, F_c, F_zc in u=z−m, avoiding cancellation at the multiple root."""
    z = acb_poly([m, 1])
    w, b = z, acb_poly([0])
    for _ in range(n):
        b, w = 2 * w * b + 1, w * w + c
    F = w - z
    return F, F.derivative(), b, b.derivative()


def _evaluate(poly: acb_poly, u: acb) -> acb:
    """Taylor shift before ball evaluation keeps cancellation near the multiple root local."""
    if bool(u.rad() > arb(2) ** -64):
        mid = acb(u.mid())
        return poly(acb_poly([mid, 1]))(u - mid)
    return poly(u)


def _count(q: int, c: acb, m: acb, r: arb, n: int | None = None) -> acb:
    n = q if n is None else n

    F, Fp, *_ = _polynomials(n, c, m)
    return _circle(r, lambda u: u * _evaluate(Fp, u) / _evaluate(F, u))


def _s_and_ds(q: int, c: acb, m: acb, r: arb, lam0: acb) -> tuple[acb, acb]:
    F, Fp, Fc, Fpc = _polynomials(q, c, m)

    def s(u):
        return u * (2 * (m + u) / lam0).log(analytic=True) * _evaluate(Fp, u) / _evaluate(F, u)

    def ds(u):
        f, fp = _evaluate(F, u), _evaluate(Fp, u)
        return u * (2 * (m + u) / lam0).log(analytic=True) * (_evaluate(Fpc, u) / f - (fp / f) * (_evaluate(Fc, u) / f))

    return _circle(r, s), _circle(r, ds)


def _is_int(x: acb, k: int) -> bool:
    return x.real.is_finite() and bool(abs(x.real - k).upper() < 0.25) and bool(abs(x.imag).upper() < 0.25)


@dataclass(frozen=True, slots=True)
class RootSide:
    p: int
    q: int
    boxes: int
    accepted: bool
    reason: str


def certify_root_side(p: int, q: int, c_root: acb, c1: acb, m: acb, r: arb, r_in: arb,
                      min_dt: Fraction = Fraction(1, 2**20)) -> RootSide:
    """Cover c_root + t·v, t ∈ [0, 1], v = c₁ − mid(c_root), by boxes satisfying 1–3 of the module docstring.

    Each box carries the root's own enclosure radius, so it covers the segment from the true root too. On
    the last box, which contains c₁, also |ρ_C| < 1 and |λ| > 1 throughout: the box is convex, so c₁ lies in
    the same component as the segment, and the attracting cycle at c₁ is C, not α."""
    old = ctx.prec
    ctx.prec = WORK_PREC
    try:
        lam0 = lambda_ball(p, q, WORK_PREC)
        # log(2z/λ₀) must be holomorphic on all of D̄′, not only on the annulus used by the contour bound:
        # |2z/λ₀ − 1| ≤ |2m/λ₀ − 1| + 2r′ < 1 keeps 2z/λ₀ in the disk |w − 1| < 1, inside Re w > 0.
        if not bool((abs(2 * m / lam0 - 1) + 2 * r_in).upper() < 1):
            return RootSide(p, q, 0, False, "log(2z/λ₀) not holomorphic on D̄′")
        v = c1 - acb(c_root.mid())
        divisors = [d for d in range(1, q) if q % d == 0]
        nbox, t, dt = 0, Fraction(0), Fraction(1, 8)
        while t < 1:
            t1 = min(t + dt, Fraction(1))
            mid = acb(c_root.mid()) + v * _arb_exact((t + t1) / 2)
            half = abs(v).upper() * _arb_exact((t1 - t) / 2) + c_root.rad() + _arb_exact(Fraction(1, 2**80))
            box = acb(arb(mid.real.mid(), half + mid.real.rad()), arb(mid.imag.mid(), half + mid.imag.rad()))
            ok, why = _box_ok(q, box, v, m, r, r_in, lam0, divisors, last=t1 == 1)
            if ok:
                nbox, t, dt = nbox + 1, t1, min(dt * 2, Fraction(1, 4))
            else:
                dt /= 2
                if dt < min_dt:
                    return RootSide(p, q, nbox, False, f"{why} near t = {float(t):.4f}")
        return RootSide(p, q, nbox, True, "accepted")
    finally:
        ctx.prec = old


def _box_ok(q, box, v, m, r, r_in, lam0, divisors, last: bool = False) -> tuple[bool, str]:
    # f(m + δ) − m = (f(m) − m) + (2m + δ)δ, so on D̄′ (|δ| ≤ r′): |f − m| ≤ |f(m) − m| + r′(|2m| + r′)
    if not bool((abs(m * m + box - m) + r_in * (abs(2 * m) + r_in)).upper() < r):
        return False, "f(D′) ⊄ D"
    if not (_is_int(_count(q, box, m, r), q + 1) and _is_int(_count(q, box, m, r_in), q + 1)):
        return False, "count"
    if not all(_is_int(_count(q, box, m, r, d), 1) for d in divisors):
        return False, "divisor count"
    root = (1 - 4 * box).sqrt()
    if not bool(root.real.lower() > 0):
        return False, "branch"
    lam, dlam = 1 - root, 2 / root
    if not bool((lam.conjugate() * dlam * v).real.lower() > 0):
        return False, "|λ| not increasing"
    s, ds = _s_and_ds(q, box, m, r_in, lam0)
    if not (s.real.is_finite() and ds.real.is_finite()):
        return False, "S"
    rho = lam0 ** (q + 1) * s.exp() / lam
    drho = rho * (ds - dlam / lam)
    if not bool((rho.conjugate() * drho * v).real.upper() < 0):
        return False, "|ρ| not decreasing"
    if last and not (bool(abs(rho).upper() < 1) and bool(abs(lam).lower() > 1)):
        return False, "last box: not |ρ| < 1 < |λ| throughout"
    return True, ""
