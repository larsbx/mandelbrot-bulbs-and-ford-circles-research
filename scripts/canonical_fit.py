"""Linear fits G = G_tip + s·δ and Re κ = κ_tip + s'·δ per canonical family, from data/canonical_profile_*.log."""
import re, glob, sys
import numpy as np
rows = []
for f in glob.glob("data/canonical_profile_*.log"):
    for line in open(f):
        m = re.match(r"(\S+)\s+a=\s*(\d+) N=\s*(\d+) q=(\d+)\s+δ=(\S+) κ=([+-][\d.]+)([+-][\d.]+)j G=([\d.]+)", line)
        if m: rows.append((m[1], int(m[2]), int(m[3]), int(m[4]), float(m[5]), float(m[6]), float(m[7]), float(m[8])))
print("family a | n | G_tip  slope_G  rms | Reκ_tip slope_Reκ | Imκ(tip)  | δ-range")
for fam in ("1/3-", "1/3+", "1/2-", "0+"):
    for a in (16, 5):
        R = [r for r in rows if r[0] == fam and r[1] == a and r[4] < 0.01]
        if len(R) < 4: continue
        d = np.array([r[4] for r in R]); G = np.array([r[7] for r in R]); kr = np.array([r[5] for r in R]); ki = np.array([r[6] for r in R])
        A = np.c_[np.ones_like(d), d]
        (g0, sg), *_ = np.linalg.lstsq(A, G, rcond=None); rms = np.sqrt(np.mean((A @ [g0, sg] - G) ** 2))
        (k0, sk), *_ = np.linalg.lstsq(A, kr, rcond=None)
        print(f"{fam:<5} {a:>2} | {len(R)} | {g0:.5f}  {sg:.3f}  {rms:.1e} | {k0:+.5f} {sk:+.3f} | {ki[np.argmin(d)]:+.5f} | [{d.min():.1e}, {d.max():.1e}]")
