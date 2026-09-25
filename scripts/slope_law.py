"""Collect converged (a=16) family slopes s and tips G_tip at rationals from data/tip_map_*.log and
data/canonical_fit.txt; fit s ∝ q'^α and s vs G_tip·q'."""
import re, glob
import numpy as np
from fractions import Fraction
pts = {}
for f in glob.glob("data/tip_map_*.log"):
    for line in open(f):
        m = re.match(r"FIT (\d+/\d+)([+−]) a=\s*16: G_tip=([\d.]+) slope=([\d.]+)", line)
        if m: pts[(Fraction(m[1]), m[2])] = (float(m[3]), float(m[4]))
for line in open("data/canonical_fit.txt"):
    m = re.match(r"(1/3|1/2|0)([+-])\s+16 \|\s+\d+ \| ([\d.]+)\s+([\d.]+)", line)
    if m: pts[(Fraction(m[1]), m[2].replace('-', '−'))] = (float(m[3]), float(m[4]))
rows = sorted(pts.items(), key=lambda kv: (kv[0][0].denominator, kv[0][0]))
print("rational side | q' | G_tip | slope s | s/q' | s/(q'·(G_tip−1))")
for (r, side), (g, s) in rows:
    qd = r.denominator if r != 0 else 1
    print(f"  {str(r):>4}{side} | {qd} | {g:.5f} | {s:.3f} | {s/qd:.3f} | {s/(qd*(g-1)):.2f}")
qs = np.array([r.denominator if r != 0 else 1 for (r, _), _ in rows], float); ss = np.array([s for _, (_, s) in rows]); gs = np.array([g for _, (g, _) in rows])
m = qs >= 2
al, c = np.polyfit(np.log(qs[m]), np.log(ss[m]), 1)
print(f"\nfit s = {np.exp(c):.3f}·q'^{al:.2f} over q' ≥ 2 (n={m.sum()}); rms of log-residual {np.sqrt(np.mean((np.log(ss[m]) - (al*np.log(qs[m])+c))**2)):.3f}")
A = np.c_[qs[m] * (gs[m] - 1)]; k = np.linalg.lstsq(A, ss[m], rcond=None)[0][0]
print(f"fit s = {k:.2f}·q'·(G_tip − 1): rms {np.sqrt(np.mean((A[:,0]*k - ss[m])**2)):.3f}  (sd of s = {ss[m].std():.3f})")
