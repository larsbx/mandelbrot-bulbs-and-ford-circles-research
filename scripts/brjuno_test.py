"""Is log|a(p/q)|/q (leading parabolic coefficient, V25) a Brjuno-type sum?  Compare with
B_fin(x) = Σ_{k=0}^{n-1} β_{k-1} log(1/α_k)  (Brjuno function at the rational x, infinite last term dropped),
evaluated at x = p/q and at x = x̃ = p̄/q, and with log q / q."""
import json, sys, cmath, math
import numpy as np
from fractions import Fraction

def brjuno_fin(x: Fraction) -> float:
    total, beta, a = 0.0, 1.0, x
    while a != 0:
        total += beta * math.log(1 / float(a)); beta *= float(a); a = (1 / a) % 1
    return total

q = int(sys.argv[1]) if len(sys.argv) > 1 else 59
rows = json.load(open(f"data/leading_coeff_q{q}.json"))
L = np.array([r["loga"] / q for r in rows]); p = np.array([r["p"] for r in rows])
xt = np.array([abs(r["xt"]) for r in rows])
Bt = np.array([brjuno_fin(Fraction(int(pp), q)) for pp in p]); Bx = np.array([brjuno_fin(Fraction(int(round(x * q)), q)) for x in xt])
print(f"q={q}: corr(log|a|/q, B_fin(p/q)/q) = {np.corrcoef(L, Bt/q)[0,1]:.4f};  corr with B_fin(x̃)/q = {np.corrcoef(L, Bx/q)[0,1]:.4f};  corr with x*: {np.corrcoef(L, xt)[0,1]:.4f}")
A = np.c_[np.ones_like(L), Bt / q]; coef, res, *_ = np.linalg.lstsq(A, L, rcond=None)
print(f"  fit log|a|/q = {coef[0]:.4f} + {coef[1]:.4f}·B_fin(p/q)/q ; rms resid {np.sqrt(np.mean((A@coef-L)**2)):.4f} (sd of L {L.std():.4f})")
A2 = np.c_[np.ones_like(L), Bt / q, Bx / q]; coef2, *_ = np.linalg.lstsq(A2, L, rcond=None)
print(f"  fit with both: {coef2[0]:.4f} + {coef2[1]:.4f}·B(p/q)/q + {coef2[2]:.4f}·B(x̃)/q ; rms resid {np.sqrt(np.mean((A2@coef2-L)**2)):.4f}")
for r, l, b, bx in sorted(zip(rows, L, Bt, Bx), key=lambda z: z[0]["p"])[:12]:
    print(f"    p={r['p']:>3} x*={abs(r['xt']):.4f}  log|a|/q={l:.4f}  B(p/q)={b:.3f}  B(x̃)={bx:.3f}")
