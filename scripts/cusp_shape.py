"""Cusp profile at rational x0: fit y = y0 + s·|x−x0|^α on each side (0.0015 < |x−x0| < w), α free (grid + lstsq)."""
import json, sys
import numpy as np
from fractions import Fraction
q = int(sys.argv[1]) if len(sys.argv) > 1 else 2003
w = float(sys.argv[2]) if len(sys.argv) > 2 else 0.03
rows = [r for r in json.load(open(f"data/kappa_q{q}.json")) if r["p"] > 5]
xt = np.array([r["xt"] for r in rows]); G = np.array([r["G"] for r in rows]); k = np.array([complex(*r["kappa"]) for r in rows])
def fit(y, x, x0, sgn):
    m = (sgn * (x - x0) > 0.0015) & (sgn * (x - x0) < w)
    d, yy = np.abs(x[m] - x0), y[m]; best = None
    for a in np.linspace(0.2, 1.5, 131):
        A = np.c_[np.ones_like(d), d ** a]; c, res, *_ = np.linalg.lstsq(A, yy, rcond=None)
        r = np.sqrt(np.mean((A @c - yy) ** 2))
        if best is None or r < best[0]: best = (r, a, c[0], c[1], m.sum())
    return best
print(f"q={q}, window {w}: rational | quantity | side | α | y0 | s | rms | n")
for r_ in [Fraction(1,3), Fraction(1,4), Fraction(2,5), Fraction(1,5), Fraction(3,7)]:
    x0 = float(r_)
    for name, y, x in (("G", G, np.abs(xt)), ("Re κ", k.real, np.abs(xt))):
        for sgn, lab in ((-1, "−"), (+1, "+")):
            b = fit(y, x, x0, sgn)
            print(f"  {str(r_):>4} | {name:<4} | {lab} | α={b[1]:.2f} | y0={b[2]:+.5f} | s={b[3]:+.4f} | rms={b[0]:.1e} | n={b[4]}")
x = np.abs(xt); b = fit(G, x, 0.5, -1); print(f"   1/2 | G    | − | α={b[1]:.2f} | y0={b[2]:+.5f} | s={b[3]:+.4f} | rms={b[0]:.1e} | n={b[4]}")
b = fit(k.real, x, 0.5, -1); print(f"   1/2 | Re κ | − | α={b[1]:.2f} | y0={b[2]:+.5f} | s={b[3]:+.4f} | rms={b[0]:.1e} | n={b[4]}")
