"""Dense-q analysis (experiments/data/kappa_q<q>.json): additive model G ≈ F(x*) + H(t), local residuals vs t,
and the local structure of κ, G around rational x*."""
from paths import DATA
import json, sys
import numpy as np
from fractions import Fraction

q = int(sys.argv[1]) if len(sys.argv) > 1 else 1009
rows = json.load(open(DATA / f"kappa_q{q}.json"))
xt = np.array([r["xt"] for r in rows]); xs = np.abs(xt); t = np.array([r["p"] / q for r in rows])
G = np.array([r["G"] for r in rows]); k = np.array([complex(*r["kappa"]) for r in rows]); p = np.array([r["p"] for r in rows])

def binned_fit(y, x, nb):
    b = np.minimum((x * nb).astype(int), nb - 1)
    means = np.array([y[b == i].mean() if (b == i).any() else 0.0 for i in range(nb)])
    return means[b]

def backfit(y, x1, x2, nb1, nb2, it=50):
    f = binned_fit(y, x1, nb1); h = np.zeros_like(y)
    for _ in range(it):
        h = binned_fit(y - f, x2, nb2); h -= h.mean()
        f = binned_fit(y - h, x1, nb1)
    return f, h

var = G.var()
fA = binned_fit(G, xs * 2, 100)                              # x* ∈ (0, ½] → 100 bins
fB, hB = backfit(G, xs * 2, t * 2, 100, 50)                  # t ∈ (0, ½] → 50 bins
print(f"q={q}: n={len(G)}  var(G)={var:.2e}")
print(f"  x*-only (100 bins):        residual var {np.var(G - fA):.2e}  (R²={1 - np.var(G - fA) / var:.4f})")
print(f"  x* + t additive (100+50):  residual var {np.var(G - fB - hB):.2e}  (R²={1 - np.var(G - fB - hB) / var:.4f})")
print(f"  sd of fitted H(t): {hB.std():.2e};  sd of x*-only residual: {(G - fA).std():.2e}")

# local residual: G minus mean of the 8 nearest x* neighbours (excluding self)
order = np.argsort(xs); e = np.zeros_like(G)
for i, idx in enumerate(order):
    nb = order[max(0, i - 4): i] .tolist() + order[i + 1: i + 5].tolist()
    e[idx] = G[idx] - G[nb].mean()
print("\nlocal residual e (G − mean of 8 x*-neighbours) by t-class:")
special = {"t≈0 (p≤5)": p <= 5, "t≈1/2 (|p−q/2|≤2)": np.abs(p - q / 2) <= 2, "t≈1/3 (|p−q/3|≤2)": np.abs(p - q / 3) <= 2,
           "t≈1/4": np.abs(p - q / 4) <= 2, "t≈2/5": np.abs(p - 2 * q / 5) <= 2, "t≈1/5": np.abs(p - q / 5) <= 2}
for name, m in special.items():
    print(f"  {name:<22} n={m.sum():>2}  mean e = {e[m].mean():+.5f}  (values: {np.round(e[m], 4).tolist()})")
print(f"  {'all':<22} n={len(e):>2}  sd e = {e.std():.5f}")
for pp in (1, 2, 3, 4, 5, 6, 7, 8):
    i = np.where(p == pp)[0]
    if len(i): print(f"  p={pp}: x*={xs[i[0]]:.4f} t={t[i[0]]:.4f} G={G[i[0]]:.5f} e={e[i[0]]:+.5f}")

print("\nlocal structure around rational x* (signed x̃; window ±0.012):")
for r_ in (Fraction(1, 2), Fraction(1, 3), Fraction(1, 4), Fraction(2, 5), Fraction(1, 5)):
    x0 = float(r_)
    for sign in (+1, -1):
        m = (np.abs(xt - sign * x0) < 0.012)
        if not m.any(): continue
        idx = np.argsort(xt[m])
        print(f"  x̃ near {'+' if sign > 0 else '-'}{r_}:")
        for j in idx:
            print(f"     x̃={xt[m][j]:+.5f} p={p[m][j]:>4}  κ={k[m][j].real:+.5f}{k[m][j].imag:+.5f}i  G={G[m][j]:.5f}")
