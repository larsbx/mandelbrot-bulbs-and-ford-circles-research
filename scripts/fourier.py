"""Fourier coefficients of κ̂ and Ĝ on the full grid x̃ = j/q (symmetrised by κ(−x̃) = conj κ(x̃)):
Re κ̂ even → a_m (cos), Im κ̂ odd → b_m (sin), Ĝ even → g_m (cos).  Tests C14″: b_m ~ 1/m (jumps), a_m ~ 1/m² (cusps)."""
import json, sys
import numpy as np
q = int(sys.argv[1]) if len(sys.argv) > 1 else 2003
rows = json.load(open(f"data/kappa_q{q}.json"))
x = np.array([r["xt"] for r in rows]); k = np.array([complex(*r["kappa"]) for r in rows]); G = np.array([r["G"] for r in rows])
X = np.r_[x, -x]; K = np.r_[k, np.conj(k)]; GG = np.r_[G, G]                     # full grid j/q, j = 1..q−1
assert len(set(np.round(X * q).astype(int))) == q - 1
M = 120
a = np.array([np.mean(K.real * np.cos(2 * np.pi * m * X)) * 2 for m in range(1, M + 1)])
b = np.array([np.mean(K.imag * np.sin(2 * np.pi * m * X)) * 2 for m in range(1, M + 1)])
g = np.array([np.mean((GG - GG.mean()) * np.cos(2 * np.pi * m * X)) * 2 for m in range(1, M + 1)])
print(f"q={q}: mean Re κ = {K.real.mean():.5f}, mean G = {GG.mean():.5f}")
print("  m |   a_m (Re κ, cos)  m²·a_m |   b_m (Im κ, sin)   m·b_m |   g_m (G, cos)   m²·g_m")
for m in list(range(1, 13)) + [16, 20, 24, 30, 40, 60, 90, 120]:
    print(f" {m:>3} | {a[m-1]:+.5f}  {m*m*a[m-1]:+.4f} | {b[m-1]:+.5f}  {m*b[m-1]:+.4f} | {g[m-1]:+.5f}  {m*m*g[m-1]:+.4f}")
ms = np.arange(1, M + 1)
for name, c in (("a_m", a), ("b_m", b), ("g_m", g)):
    sel = (ms >= 3) & (ms <= 60) & (np.abs(c) > 0)
    al, lc = np.polyfit(np.log(ms[sel]), np.log(np.abs(c[sel])), 1)
    print(f"  |{name}| ~ m^{al:.2f} (fit 3 ≤ m ≤ 60, prefactor {np.exp(lc):.4f})")
np.save(f"data/fourier_q{q}.npy", np.c_[ms, a, b, g])
