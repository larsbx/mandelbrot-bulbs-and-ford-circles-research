"""V52: the root principle at a primitive (cusp) root.  The real period-3 component P₃ (centre ≈ −1.75488, root
c = −7/4): κ(1/q) of its satellites, computed exactly as (ι(F^q) − ½)/q (Theorem 2.1 / Remark 2.2), where F is the
first-return germ f_c³ at the cycle point at the satellite root c(λ₀), λ₀ = e^{2πi/q}, and ι(F^q) comes from the
general-germ normal form (germ.normal_form_index, O(q³)); against the horn map of the normalized cusp germ of f³ at
c = −7/4.  Writes data/primitive_root.txt."""
import sys, time
import mpmath as mp
from bulbford.germ import normal_form_index
from bulbford.horn import horn_coeffs
from bulbford.extrapolate import power_fit

N_PER, C_CENTRE = 3, mp.mpf("-1.7548776662466927600495088963585286918946")


def iterate_series(c, z, L):
    """Coefficients of f_c^3(z + w) − z in w (degree ≤ L)."""
    S = [z, mp.mpf(1)] + [0] * (L - 1)
    for _ in range(N_PER):
        T = [mp.fsum(S[i] * S[k - i] for i in range(k + 1)) for k in range(L + 1)]
        T[0] += c
        S = T
    return [S[0] - z] + S[1:]


def root_cycle(lam, steps=40):
    """(c, z) with f_c³(z) = z, (f_c³)'(z) = lam, continued from (centre, 0) at λ = 0 (2×2 Newton in mpmath)."""
    c, z = mp.mpc(C_CENTRE), mp.mpc(0)
    for t in range(1, steps + 1):
        target = lam * t / steps
        for _ in range(50):
            S = iterate_series(c, z, 2)
            h = mp.mpf(10) ** (-mp.mp.dps // 2)
            Sc = iterate_series(c + h, z, 2)
            f1, f2 = S[0], S[1] - target
            a, e = S[1] - 1, 2 * S[2]                            # ∂/∂z of (f³ − z) and of (f³)'
            b, g = (Sc[0] - S[0]) / h, (Sc[1] - S[1]) / h        # ∂/∂c (forward difference, refined by Newton)
            det = a * g - b * e
            dz, dc = (f1 * g - b * f2) / det, (a * f2 - e * f1) / det
            z, c = z - dz, c - dc
            if abs(dz) + abs(dc) < mp.mpf(10) ** (-mp.mp.dps + 8):
                break
    return c, z


def kappa_primitive(q, dps=60):
    with mp.workdps(dps):
        lam = mp.expjpi(mp.mpf(2) / q)
        c, z = root_cycle(lam)
        g = iterate_series(c, z, 2 * q + 2)                      # the polynomial F has degree 8; higher terms vanish
        return (normal_form_index(g, q) - mp.mpf(1) / 2) / q


def cusp_constant(dps=50):
    with mp.workdps(dps):
        c = mp.mpf(-7) / 4
        z = mp.findroot(lambda z: iterate_series(c, z, 1)[1] - 1, mp.mpf("-0.03"))
        S = iterate_series(c, z, 8)
        germ = tuple(mp.mpc(S[k] / S[2] ** (k - 1)) for k in range(2, 9))
        a, _ = horn_coeffs(2, h=0.25, dps=dps, germ=germ)
        return germ, a[2] / (2j * mp.pi * a[1] ** 2)


if __name__ == "__main__":
    mp.mp.dps = 40
    germ, K = cusp_constant()
    with open("data/primitive_root.txt", "w") as out:
        def emit(s):
            out.write(s + "\n"); out.flush(); print(s, flush=True)
        emit(f"cusp germ at c = -7/4: g3 = {mp.nstr(mp.re(germ[1]), 15)} (= -2/49?), horn kappa = {mp.nstr(K, 16)}")
        rows = []
        for q in map(int, sys.argv[1:] or (32, 48, 64, 96, 128, 192)):
            t = time.time()
            k = kappa_primitive(q)
            rows.append((q, k))
            emit(f"q = {q}: kappa(1/q) = {mp.nstr(k, 20)}   ({time.time() - t:.0f} s)")
        for J in range(2, len(rows) - 1):
            f = power_fit(rows, J)
            emit(f"fit J = {J}: c0 = {mp.nstr(f[0], 14)}, |c0 - horn| = {mp.nstr(abs(f[0] - K), 2)}, c1*2pi/-i = {mp.nstr(f[1] * 2 * mp.pi / (-1j), 8)}")
