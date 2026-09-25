"""Horn map of f(z) = z + z² and κ₀ = a₂/(2πi a₁²) (V38).  Run: python3 -m bulbford.horn [h …]
w = −1/z: F(w) = w²/(w−1) = w + 1 + 1/(w−1).  Formal Fatou coordinate Φ(w) = w − log w + Σ_{k≥1} c_k w^{−k}
(exact rationals from Φ∘F = Φ + 1).  Φ_att = Φ(F^n w) − n (principal log), Ψ_rep(Z) = F^n(Φ⁻¹(Z − n)) (log with
arg ∈ (0, 2π) on the repelling side).  E(Z) = Φ_att(Ψ_rep(Z)) = Z + a₀ + Σ_n a_n e^{±2πinZ} at the upper/lower end;
in W = e^{±2πiZ} the germ is W ↦ e^{±2πia₀}W(1 + A W + B W² + …) and ι = B/A²."""
import sys
from fractions import Fraction as Fr
from math import comb
import numpy as np

def fatou_coeffs(K: int) -> list:
    """c_0..c_K (c_0 = 0) of Φ(w) = w − log w + Σ c_k w^{−k}:  Σ_k c_k[(u−u²)^k − u^k] = −Σ_{m≥2}(1 − 1/m)u^m."""
    c = [Fr(0)] * (K + 1)
    for m in range(2, K + 2):
        s = Fr(m - 1, m) + sum(c[k] * comb(k, m - k) * (-1) ** (m - k) for k in range(1, m - 1) if m - k <= k)
        c[m - 1] = s / (m - 1)
    return c


K = 24
c = fatou_coeffs(K)
cf = np.array([float(x) for x in c])

def phi(w, branch):                              # asymptotic Fatou coordinate
    lg = np.log(w) if branch == "att" else np.log(np.abs(w)) + 1j * np.mod(np.angle(w), 2 * np.pi)
    return w - lg + sum(cf[k] * w ** (-k) for k in range(1, K + 1))

def dphi(w):
    return 1 - 1 / w - sum(k * cf[k] * w ** (-k - 1) for k in range(1, K + 1))

F = lambda w: w + 1 + 1 / (w - 1)

def phi_att(w, R=4e4, maxit=400000):
    n = np.zeros(w.shape); w = w.copy()
    for _ in range(maxit):
        act = (np.abs(w) < R) | (w.real < 0.5 * np.abs(w))
        if not act.any(): break
        w[act] = F(w[act]); n[act] += 1
    return phi(w, "att") - n

def psi_rep(Z, n=40000):
    Zp = Z - n; w = Zp + np.log(-Zp) + 0j        # initial guess on the far left
    for _ in range(60):
        w = w - (phi(w, "rep") - Zp) / dphi(w)
    for _ in range(n):
        w = F(w)
    return w

if __name__ == "__main__":
    N = 64
    for h in map(float, sys.argv[1:] or ["0.6", "0.8", "1.0"]):
        for end in (+1, -1):
            x = np.arange(N) / N; Z = x + 1j * end * h
            D = phi_att(psi_rep(Z)) - Z                    # = a₀ + Σ a_n e^{±2πinZ}
            co = np.fft.fft(D) / N                         # co[n] ↔ e^{2πinx}
            if end > 0:  a1, a2 = co[1] * np.exp(2 * np.pi * h), co[2] * np.exp(4 * np.pi * h); s = 2j * np.pi
            else:        a1, a2 = co[-1] * np.exp(2 * np.pi * h), co[-2] * np.exp(4 * np.pi * h); s = -2j * np.pi
            A, B = s * a1, s * a2 + (s * a1) ** 2 / 2
            print(f"h={h} {'upper' if end > 0 else 'lower'}: a0={co[0]:.10f}  a1={a1:.8f}  a2={a2:.8f}  ι−½ = {B / A**2 - 0.5:.10f}  |a3·e^-6πh|={abs(co[3 if end > 0 else -3]):.1e}", flush=True)


def horn_coeffs(M, h=0.25, N=64, R=400, n=500, dps=40, K=30):
    """Upper horn-map Fourier coefficients a_0..a_M of E(Z) = Z + a₀ + Σ a_n e^{2πinZ} (mpmath), plus aliasing |ĉ_{N/2}|."""
    import mpmath as mp
    with mp.workdps(dps):
        cm = [mp.mpf(x.numerator) / x.denominator for x in fatou_coeffs(K)]
        def ph(w, rep):
            lg = mp.log(w) if not rep else mp.log(abs(w)) + 1j * (mp.arg(w) % (2 * mp.pi))
            return w - lg + mp.fsum(cm[k] * w ** (-k) for k in range(1, K + 1))
        dph = lambda w: 1 - 1 / w - mp.fsum(k * cm[k] * w ** (-k - 1) for k in range(1, K + 1))
        Fm = lambda w: w + 1 + 1 / (w - 1)
        def E(Z):
            Zp = Z - n; w = Zp + mp.log(-Zp)
            for _ in range(80):
                s = (ph(w, True) - Zp) / dph(w); w -= s
                if abs(s) < mp.mpf(10) ** (-dps + 5): break
            for _ in range(n): w = Fm(w)
            m = 0
            while abs(w) < R or w.real < abs(w) / 2:
                w = Fm(w); m += 1
            return ph(w, False) - m
        Zs = [mp.mpf(j) / N + 1j * mp.mpf(h) for j in range(N)]
        D = [E(Z) - Z for Z in Zs]
        coef = lambda k: mp.fsum(D[j] * mp.expjpi(-2 * mp.mpf(k) * j / N) for j in range(N)) / N
        return [coef(k) * mp.exp(2 * mp.pi * k * h) for k in range(M + 1)], abs(coef(N // 2))


def kappa0_mp(h=0.25, N=48, R=400, n=500, dps=40, K=K):
    """(a₁, a₂, κ₀ = a₂/(2πi a₁²), aliasing) for the upper horn map."""
    import mpmath as mp
    a, al = horn_coeffs(2, h=h, N=N, R=R, n=n, dps=dps, K=K)
    with mp.workdps(dps):
        return a[1], a[2], a[2] / (2j * mp.pi * a[1] ** 2), al
