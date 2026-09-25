"""Mixed jump + cusp inversion for the even functions Ĝ and Re κ̂:
 a_m = (2/(π m)) Σ_{q'} Jr(q') S_{q'}(m) + Σ_{q'≤Qc} A(q') c_{q'}(m),
 S_{q'}(m) = Σ_{0<p'<q'/2, p'⟂q'} sin(2πm p'/q')   (jump +Jr at p'/q', −Jr at −p'/q'),
 c_{q'}(m) Ramanujan sums (cusp area, flat in m)."""
import sys
import numpy as np
from math import gcd, pi
from sympy import mobius, divisors
q = int(sys.argv[1]) if len(sys.argv) > 1 else 2003
Q, Qc, M = 24, 6, 120
d = np.load(f"data/fourier_q{q}.npy"); ms = d[:M, 0].astype(int); a, g = d[:M, 1], d[:M, 3]
S = np.array([[sum(np.sin(2 * pi * m * pp / qq) for pp in range(1, (qq + 1) // 2) if gcd(pp, qq) == 1 and 2 * pp != qq) for qq in range(2, Q + 1)] for m in ms])
C = np.array([[sum(dd * int(mobius(qq // dd)) for dd in divisors(gcd(qq, m))) for qq in range(2, Qc + 1)] for m in ms], float)
X = np.c_[2 / (pi * ms[:, None]) * S, C]
tips = {3: -0.0020, 4: -0.0014, 5: -0.0017}       # V32/V30 Ĝ jumps (⅓, ¼, ⅖)
for name, y in (("Ĝ", g), ("Re κ̂", a)):
    x, res, rk, sv = np.linalg.lstsq(X, y, rcond=None)
    print(f"{name}: relative residual {np.linalg.norm(y - X @ x) / np.linalg.norm(y):.3f}, cond {sv[0]/sv[-1]:.0f}")
    print("   q' | jump Jr(q')   | cusp area A(q')  | local tip jump (V32)")
    for i, qq in enumerate(range(2, Q + 1)):
        if qq <= 9:
            A = f"{x[Q - 1 + qq - 2]:+.5f}" if qq <= Qc else "   —   "
            print(f"   {qq:>2} | {x[i]:+.5f}      | {A}        | {tips.get(qq, ''):}")
