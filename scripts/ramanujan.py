"""Ramanujan inversion: m·b_m = Σ_{q'} β(q') c_{q'}(m),  m²·a_m = Σ α(q') c_{q'}(m),  m²·g_m = Σ γ(q') c_{q'}(m)
(least squares, q' ≤ Q, 1 ≤ m ≤ M).  Jump model predicts β(q') = J(q')/π with J the Im κ̂ jump (V23: ≈ −0.05/q')."""
import sys
import numpy as np
from math import gcd
from sympy import mobius, divisors
q = int(sys.argv[1]) if len(sys.argv) > 1 else 2003
Q, M = int(sys.argv[2]) if len(sys.argv) > 2 else 24, 120
d = np.load(f"data/fourier_q{q}.npy"); ms, a, b, g = d[:M, 0], d[:M, 1], d[:M, 2], d[:M, 3]
def cq(qq, m): return sum(dd * int(mobius(qq // dd)) for dd in divisors(gcd(qq, m)))
C = np.array([[cq(qq, int(m)) for qq in range(2, Q + 1)] for m in ms], float)
sol = {}
for name, y in (("β (m·b_m)", ms * b), ("α (m²·a_m)", ms ** 2 * a), ("γ (m²·g_m)", ms ** 2 * g)):
    x, res, rk, sv = np.linalg.lstsq(C, y, rcond=None); fit = C @ x
    sol[name] = x
    print(f"{name}: relative residual {np.linalg.norm(y - fit) / np.linalg.norm(y):.3f}  (cond {sv[0]/sv[-1]:.1f})")
J = {2: -0.034, 3: -0.0210, 4: -0.0146, 5: -0.0104, 6: -0.0092, 7: -0.0046, 8: -0.0053}     # V23 (q=2003 intercept fits; ½ from V21)
print("\n q' |   β(q')   J(q')/π  ratio |   α(q')  α·q'² |   γ(q')  γ·q'²")
for i, qq in enumerate(range(2, Q + 1)):
    jj = J.get(qq); jp = f"{jj/np.pi:+.4f}  {sol['β (m·b_m)'][i]/(jj/np.pi):5.2f}" if jj else "        —        "
    print(f" {qq:>2} | {sol['β (m·b_m)'][i]:+.4f}  {jp} | {sol['α (m²·a_m)'][i]:+.4f} {sol['α (m²·a_m)'][i]*qq*qq:+.3f} | {sol['γ (m²·g_m)'][i]:+.4f} {sol['γ (m²·g_m)'][i]*qq*qq:+.3f}")
