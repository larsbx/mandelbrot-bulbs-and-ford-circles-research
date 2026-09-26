"""Extrapolate data/level_diag.tsv in 1/N and compare with the hierarchy (V44, paper Conjecture 6.3):
[0; N, N] → conj κ₂ (control, known from V42), [0; N, N, N] → κ₃.  Prints fits for several J and N-windows."""
import mpmath as mp
from bulbford.extrapolate import power_fit

mp.mp.dps = 40
rows = [l.split("\t") for l in open("data/level_diag.tsv") if l.strip()]
data = {k: sorted((int(r[1]), mp.mpc(r[4], r[5])) for r in rows if int(r[0]) == k) for k in (2, 3)}
lv = [mp.mpc(*l.split("\t")[1:3]) for l in open("data/levels_h100.tsv") if l.strip()]
target = {2: mp.conj(lv[1]), 3: lv[2]}
for k, pts in data.items():
    print(f"k = {k}: {len(pts)} points, N = {pts[0][0]}…{pts[-1][0]}; target {mp.nstr(target[k], 10)}; "
          f"alternatives κ₂ = {mp.nstr(lv[1], 8)}, κ̄₃ = {mp.nstr(mp.conj(lv[2]), 8)}")
    for skip in (0, 2, 4):
        for J in (3, 4, 5, 6):
            use = pts[skip:]
            if len(use) <= J + 1:
                continue
            c0 = power_fit(use, J)[0]
            print(f"  N ≥ {use[0][0]:3d}, J = {J}: {mp.nstr(c0, 10)}   |Δ target| = {mp.nstr(abs(c0 - target[k]), 3)}")
