"""Horn maps of the parabolic germs f(z) = e^{2πip/q}z + z² (the root of the cardioid's p/q-satellite) and the one-sided
limits of κ at rationals (V50).  g = f^q = z(1 + h(z)), h = a z^q + …, has q attracting and q repelling petals, and f
permutes them (a step of 1/q in Fatou coordinates).

Formal Fatou coordinate of g:  Φ(z) = Σ_{k=−q, k≠0}^{K} s_k z^k + β log z,  Φ∘g = Φ + 1  (β = résit(g); the equation at
order z^{k+q} fixes s_k with the factor k·a, the one at z^q fixes β with the factor a).  Horn maps from the repelling petal
R₀ that follows the base attracting petal P₀ counterclockwise: Ψ_rep(Z) on R₀, pushed by g into some attracting petal P_j,
carried back to P₀ by f^m (m ≡ −j·p⁻¹ mod q), E(Z) = Φ_att(f^m(g^·(Ψ_rep(Z)))) − m/q.  Upper map: Im Z ≫ 0 (expansion in
e^{2πinZ}); lower map: Im Z ≪ 0 (expansion in e^{−2πinZ}).  q = 1 is the cusp z + z² (κ₀), q = 2 the ½-root (κ_½)."""
from __future__ import annotations
from functools import lru_cache
import mpmath as mp


def _mul(a, b, L):
    return [mp.fsum(a[i] * b[k - i] for i in range(k + 1)) for k in range(L + 1)]


def _log1p(h, L):
    out, term = [mp.mpc(0)] * (L + 1), [mp.mpc(1)] + [mp.mpc(0)] * L
    for j in range(1, L + 1):
        term = _mul(term, h, L)
        if all(t == 0 for t in term):
            break
        out = [o + (-1) ** (j + 1) * t / j for o, t in zip(out, term)]
    return out


def _exp(x, L):
    """exp of a series with zero constant term (E' = x'E)."""
    e = [mp.mpc(1)] + [mp.mpc(0)] * L
    for k in range(1, L + 1):
        e[k] = mp.fsum(j * x[j] * e[k - j] for j in range(1, k + 1)) / k
    return e


def poly_germ(coeffs, digits: int = 120) -> tuple:
    """Hashable germ spec (for the functions below) from coefficients [f₁, f₂, …] of f(z) = Σ f_k z^k, f₁ a primitive
    q-th root of unity; stored as decimal strings, re-read at the working precision."""
    with mp.workdps(digits):
        return tuple((mp.nstr(mp.re(c), digits), mp.nstr(mp.im(c), digits)) for c in coeffs)


def two_cycle_germ(c) -> tuple:
    """First return map f² of z² + c at a point z₀ of its 2-cycle, in w = z − z₀:
    F(w) = 4z₀z₁·w + (2z₁ + 4z₀²)·w² + 4z₀·w³ + w⁴ (z₁ = f(z₀); multiplier 4z₀z₁ = 4(c + 1))."""
    with mp.workdps(130):
        c = mp.mpc(c)
        z0 = (-1 + mp.sqrt(-3 - 4 * c)) / 2                   # roots of z² + z + c + 1
        z1 = z0 * z0 + c
        return poly_germ([4 * z0 * z1, 2 * z1 + 4 * z0 ** 2, 4 * z0, mp.mpc(1)])


def _poly_coeffs(p: int, q: int, poly):
    lam = mp.expjpi(2 * mp.mpf(p) / q)
    if poly is None:
        return [lam, mp.mpc(1)]                               # e^{2πip/q}z + z²
    cs = [mp.mpc(*x) for x in poly]
    assert abs(cs[0] - lam) < mp.mpf(10) ** -20, "multiplier must be e^{2πip/q}"
    return cs


def _orbit_map(p: int, q: int, poly=None):
    cs = _poly_coeffs(p, q, poly)
    def f(z):
        out = mp.mpc(0)
        for c in reversed(cs):
            out = (out + c) * z
        return out
    def g(z):
        for _ in range(q):
            z = f(z)
        return z
    return cs[0], f, g


@lru_cache(maxsize=None)
def fatou_series_pq(p: int, q: int, K: int = 30, poly=None):
    """({k: s_k}, β, g) at the current precision; g = f^q as a callable (f = e^{2πip/q}z + z², or the polynomial germ
    `poly` from poly_germ)."""
    _, f, g = _orbit_map(p, q, poly)
    cs = _poly_coeffs(p, q, poly)
    L = K + q + 2
    F = [mp.mpc(0), mp.mpc(1)] + [mp.mpc(0)] * (L)
    for _ in range(q):                                        # coefficients of f^q up to z^{L+1}: F ← Σ c_k F^k
        new, pw = [mp.mpc(0)] * (L + 2), F
        for c in cs:
            new = [a + c * b for a, b in zip(new, pw)]
            pw = _mul(pw, F, L + 1)
        F = new
    h = [mp.mpc(0)] + F[2:L + 2]                              # g = z(1 + h)
    lg = _log1p(h, L)
    P = {k: [e - (1 if i == 0 else 0) for i, e in enumerate(_exp([k * c for c in lg], L))] for k in range(-q, K + 1) if k}
    s, beta = {}, mp.mpc(0)
    for m in range(K + q + 1):                                # order z^m of Φ∘g − Φ − 1
        known = mp.fsum(c * P[k][m - k] for k, c in s.items() if 0 <= m - k <= L) + beta * lg[m]
        rhs = (1 if m == 0 else 0) - known
        if m == q:
            beta = rhs / lg[q]
        else:
            s[m - q] = rhs / P[m - q][q]                      # [z^q]((1+h)^k − 1) = k·a
    return s, beta, g


def _geometry(p: int, q: int, K: int, n: int, poly=None):
    """Petal geometry at the current precision: (s, g, f, θ₀, Φ, Φ'), with θ₀ the base attracting direction and Φ on a
    branch of log continuous on the sector from P₀ to R₀ (centred between them)."""
    s, beta, g = fatou_series_pq(p, q, K, poly)
    _, f, _ = _orbit_map(p, q, poly)
    a = -1 / (q * s[-q])                                      # s_{−q} = −1/(q a)
    th0 = (mp.pi - mp.arg(a)) / q
    mid = th0 + mp.pi / (2 * q)
    rot = mp.expj(-mid)
    Phi = lambda z: mp.fsum(c * z ** k for k, c in s.items()) + beta * (mp.log(z * rot) + 1j * mid)
    dPhi = lambda z: mp.fsum(k * c * z ** (k - 1) for k, c in s.items()) + beta / z
    return s, g, f, th0, Phi, dPhi


def _start(Z, n, q, s, th0, Phi, dPhi, tol):
    """Ψ_rep(Z) at depth n: the solution of Φ(z) = Z − n on the repelling petal R₀ (direction θ₀ + π/q)."""
    Zp = Z - n
    r = (s[-q] / Zp) ** (mp.mpf(1) / q)
    th_rep = th0 + mp.pi / q
    z = min((r * mp.expjpi(2 * mp.mpf(j) / q) for j in range(q)), key=lambda w: abs(mp.arg(w * mp.expj(-th_rep))))
    for _ in range(100):
        d = (Phi(z) - Zp) / dPhi(z); z -= d
        if abs(d) < tol * abs(z): break
    return z


def _land(z, n, q, g, th0, r_n):
    """Index j of the attracting petal reached after 2n steps of g (None if the orbit escapes); r_n = |s_{−q}/n|^{1/q}
    is the radius at depth n, which sets the scale of the germ."""
    for _ in range(2 * n):
        z = g(z)
        if abs(z) > 1000 * r_n:
            return None, z
    if abs(z) > 4 * r_n:
        return None, z
    return int(mp.nint((mp.arg(z * mp.expj(-th0)) % (2 * mp.pi)) / (2 * mp.pi / q))) % q, z


@lru_cache(maxsize=None)
def escape_band(p: int, q: int, n: int = 600, K: int = 30, lo: float = -30, hi: float = 30, step: float = 0.25,
                poly=None):
    """Heights (h₋, h₊) between which the repelling cylinder of R₀ is not mapped into the petals (orbits escape or
    switch petals).  Below h₋ orbits reach P₀ (lower horn map), above h₊ the counterclockwise petal (upper map).
    Determined by scanning Im Z (one sample per height) — the centre of the cylinder is not Im Φ on the ray, since the
    terms z^{−k}, k < q, of Φ have depth-dependent imaginary parts."""
    with mp.workdps(40):
        s, g, f, th0, Phi, dPhi = _geometry(p, q, K, n, poly)
        tol = mp.mpf(10) ** -35
        r_n = abs(s[-q] / n) ** (mp.mpf(1) / q)
        pattern = []
        h = lo
        while h <= hi:
            z0 = _start(mp.mpf("0.5") + 1j * h, n, q, s, th0, Phi, dPhi, tol)
            pattern.append((h, _land(z0, n, q, g, th0, r_n)[0]))
            h += step
    bad = [h for h, j in pattern if j is None]
    if bad:
        return min(bad), max(bad)
    switch = [h for (h, j), (_, k) in zip(pattern, pattern[1:]) if j != k]
    return switch[0], switch[0] + step


def horn_coeffs_pq(p: int, q: int, M: int = 2, upper: bool = True, dh: float = 1.5, N: int = 64, n: int = 1500,
                   K: int = 30, dps: int = 60, h: float | None = None, poly=None):
    """Fourier coefficients (a_0..a_M) of the upper (b_0..b_M of the lower) horn map, sampled at distance dh beyond the
    escape band (escape_band) or at the absolute height h, and the aliasing level.  `dps` is the target accuracy; the
    working precision adds the digits lost to the amplification e^{2πM|h|} of the M-th mode."""
    if h is None:
        h_lo, h_hi = escape_band(p, q, poly=poly)
        h = h_hi + dh if upper else h_lo - dh
    dps += int(2 * float(mp.pi) * M * abs(h) / 2.302585) + 10
    with mp.workdps(dps):
        s, g, f, th0, Phi, dPhi = _geometry(p, q, K, n, poly)
        pinv = pow(p, -1, q) if q > 1 else 0
        tol = mp.mpf(10) ** (-dps + 5)
        r_n = abs(s[-q] / n) ** (mp.mpf(1) / q)

        def E(Z):
            j, z = _land(_start(Z, n, q, s, th0, Phi, dPhi, tol), n, q, g, th0, r_n)
            if j is None:
                raise ValueError(f"orbit of Z = {Z} escaped")
            m = (-j * pinv) % q
            for _ in range(m):
                z = f(z)
            return Phi(z) - 2 * n - mp.mpf(m) / q, j
        Zs = [mp.mpf(k) / N + 1j * mp.mpf(h) for k in range(N)]
        vals = [E(Z) for Z in Zs]
        if len({j for _, j in vals}) != 1:
            raise ValueError(f"samples land in different petals {sorted({j for _, j in vals})}: increase dh")
        D = [v - Z for (v, _), Z in zip(vals, Zs)]
        c = lambda k: mp.fsum(D[j] * mp.expjpi(-2 * mp.mpf(k) * j / N) for j in range(N)) / N
        sgn = 1 if upper else -1
        coeffs = [c(sgn * k) * mp.exp(2 * mp.pi * k * h * sgn) for k in range(M + 1)]
        return coeffs, abs(c(N // 2))


def horn_index_pq(p: int, q: int, upper: bool = True, **kw):
    """ι − ½ of the horn germ: a₂/(2πi a₁²) for the upper map, −b₂/(2πi b₁²) for the lower map (W = e^{−2πiZ})."""
    dps = kw.get("dps", 60)
    a, _ = horn_coeffs_pq(p, q, upper=upper, **kw)
    with mp.workdps(dps):
        return (1 if upper else -1) * a[2] / (2j * mp.pi * a[1] ** 2)


def root(p: int, q: int, upper: bool = True, poly=None) -> tuple:
    """Germ tag for renorm.py: the upper/lower parabolic renormalization of e^{2πip/q}z + z² (the p/q-root of the
    cardioid) or of the polynomial germ `poly` (poly_germ, two_cycle_germ)."""
    return ("root", p, q, upper, poly)


def root_germ_coeffs(p: int, q: int, upper: bool, poly, M: int, dps: int) -> list:
    """[0, a_1, …, a_M] of the root germ in the form W·exp(2πi Σ a_n W^n) = W + W² + …: the lower map, written in
    W = e^{−2πiZ}, has a_n = −b_n; then Z is translated so that a₁ = 1/(2πi) (a_n ↦ a_n (2πi a₁)^{−n})."""
    a, _ = horn_coeffs_pq(p, q, M, upper=upper, dh=1.0, n=3000, K=45, dps=dps + 4 * M, poly=poly)  # high modes: close
    # to the band and deep (the error of the truncated Fatou series is amplified by e^{2πn·dh} in the n-th mode)
    with mp.workdps(dps + 4 * M):
        a = [x if upper else -x for x in a]
        sc = 2j * mp.pi * a[1]
        return [mp.mpc(0)] + [a[k] / sc ** k for k in range(1, M + 1)]
