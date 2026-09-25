"""Cusp/jump laws at rationals p'/q' from a dense sweep: for Re κ and G fit a two-sided line
y = y0 + s_∓ |x − x0| on 0.003 < |x − x0| < w on each side; cusp depth := envelope(mean of both sides at |x−x0|=w) − y0;
jump := y0⁺ − y0⁻ (fitted intercepts). Bounded-p bulbs (p ≤ 5) excluded."""
import json, sys
import numpy as np
from fractions import Fraction

q = int(sys.argv[1]) if len(sys.argv) > 1 else 1009
w = float(sys.argv[2]) if len(sys.argv) > 2 else 0.012
rows = [r for r in json.load(open(f"data/kappa_q{q}.json")) if r["p"] > 5]
xt = np.array([r["xt"] for r in rows]); G = np.array([r["G"] for r in rows]); k = np.array([complex(*r["kappa"]) for r in rows])

def side_fit(y, x, x0, sgn):
    m = (sgn * (x - x0) > 0.0025) & (sgn * (x - x0) < w)
    if m.sum() < 3: return None
    d = np.abs(x[m] - x0); A = np.c_[np.ones(m.sum()), d]
    (y0, s), *_ = np.linalg.lstsq(A, y[m], rcond=None)
    return y0, s, m.sum()

print(f"q={q}, window {w}:  rational | quantity | y0⁻  y0⁺ | jump y0⁺−y0⁻ | cusp depth (env−mean y0) | slopes ∓")
for r_ in [Fraction(1,2), Fraction(1,3), Fraction(1,4), Fraction(1,5), Fraction(2,5), Fraction(1,6), Fraction(1,7), Fraction(2,7), Fraction(3,7), Fraction(3,8)]:
    x0 = float(r_)
    for name, y in (("Re κ", k.real), ("Im κ", k.imag), ("G", G)):
        lo, hi = side_fit(y, xt, x0, -1), side_fit(y, xt, x0, +1)
        if x0 == 0.5:   # x̃ = ½ ≡ −½: the "+" side is the conjugate side; use |x̃| with signed Im
            lo = side_fit(y, np.abs(xt), 0.5, -1)
            yy = np.where(xt < 0, (-1 if name == "Im κ" else 1) * y, y)
            hi = side_fit(yy, -np.abs(xt) + 1.0, 0.5, +1)   # mirror: x = 1 − |x̃| > ½
        if lo is None or hi is None: continue
        env = 0.5 * (lo[0] + lo[1] * w + hi[0] + hi[1] * w); depth = env - 0.5 * (lo[0] + hi[0])
        print(f"  {str(r_):>4} | {name:<4} | {lo[0]:+.5f} {hi[0]:+.5f} | {hi[0]-lo[0]:+.5f} | {depth:+.5f} (·q'² = {depth*r_.denominator**2:+.4f}) | {lo[1]:+.3f} {hi[1]:+.3f}  (n={lo[2]},{hi[2]})")
