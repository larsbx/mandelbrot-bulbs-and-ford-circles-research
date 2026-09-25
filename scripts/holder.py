"""Modulus of continuity of G, Re κ, Im κ on the generic bulbs (p > 5) of a dense sweep:
ω(δ) = max |y(x) − y(x')| over |x − x'| ≤ δ (signed x̃ for Im κ, |x̃| otherwise); log-log slope = Hölder exponent."""
import json, sys
import numpy as np
q = int(sys.argv[1]) if len(sys.argv) > 1 else 1009
rows = [r for r in json.load(open(f"data/kappa_q{q}.json")) if r["p"] > 5]
xt = np.array([r["xt"] for r in rows]); G = np.array([r["G"] for r in rows]); k = np.array([complex(*r["kappa"]) for r in rows])
deltas = np.array([0.0015, 0.002, 0.003, 0.004, 0.006, 0.008, 0.012, 0.016, 0.024, 0.032, 0.048])
def omega(y, x):
    o = np.argsort(x); x, y = x[o], y[o]
    out = []
    for d in deltas:
        m = 0.0
        for i in range(len(x)):
            j = np.searchsorted(x, x[i] + d, side="right")
            if j > i + 1: m = max(m, np.abs(y[i + 1:j] - y[i]).max())
        out.append(m)
    return np.array(out)
print(f"q={q}, generic bulbs n={len(rows)}.   δ:      " + "  ".join(f"{d:.4f}" for d in deltas))
for name, y, x in (("G", G, np.abs(xt)), ("Re κ", k.real, np.abs(xt)), ("Im κ", k.imag, xt)):
    w = omega(y, x); sl = np.polyfit(np.log(deltas[2:]), np.log(w[2:]), 1)[0]
    print(f"  ω_{name:<5}: " + "  ".join(f"{v:.4f}" for v in w) + f"   log-log slope (δ ≥ 0.003) = {sl:.2f}")
# same on the median-smoothed (jump-free) part: |Im κ| to remove the odd jumps' sign
w = omega(np.abs(k.imag), np.abs(xt)); print(f"  ω_|Imκ|: " + "  ".join(f"{v:.4f}" for v in w) + f"   slope = {np.polyfit(np.log(deltas[2:]), np.log(w[2:]), 1)[0]:.2f}")
