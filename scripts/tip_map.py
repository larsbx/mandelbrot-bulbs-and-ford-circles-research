"""Tip map and family slopes at rationals x0 = [0; tail]: x̃ = [0; tail, N, a] for a ∈ A, N ∈ NS.
Fits G = G_tip + s·δ per (x0, side, a).  Output: data/tip_map_<label>.json"""
import json, sys
import numpy as np
from bulbford.cf import from_cf, modinv
from bulbford.taylor import kappa_fft
TAILS = {"1/4": (4,), "1/4+": (3, 1), "1/5": (5,), "2/5": (2, 2), "2/5+": (2, 1, 1), "2/7": (3, 2), "1/3": (3,), "1/2": (2,),
         "1/6": (6,), "1/7": (7,), "3/7": (2, 3), "3/8": (2, 1, 2), "2/9": (4, 2), "4/9": (2, 4)}
label = sys.argv[1]; names = [n for n in sys.argv[2:] if not n.startswith("a=")]
A = tuple(int(v[2:]) for v in sys.argv[2:] if v.startswith("a=")) or (2, 3, 5, 16, 64); NS = (24, 32, 48, 64, 96, 128, 192, 256)
rows = []
for name in names:
    tail = TAILS[name]; x0 = from_cf((0,) + tail); x0 = x0[0] / x0[1]
    for a in A:
        for N in NS:
            pb, q = from_cf((0,) + tail + (N, a))
            if q > 16000: continue
            p = modinv(pb, q); t = kappa_fft(p, q); G = abs(t.solve(-1, 2.0)) / 2
            xt = (pb if pb <= q / 2 else pb - q) / q; sd = abs(xt) - x0
            rows.append(dict(x0=name, a=a, N=N, p=p, q=q, xt=xt, sdelta=sd, kappa=[t.coeffs[2].real, t.coeffs[2].imag], G=G))
            print(f"{name:<4} a={a:>2} N={N:>3} q={q:<6} δ={sd:+.2e} κ={t.coeffs[2]:+.5f} G={G:.5f}", flush=True)
json.dump(rows, open(f"data/tip_map_{label}.json", "w"))
for name in names:
    for a in A:
        R = [r for r in rows if r["x0"] == name and r["a"] == a]
        if len(R) < 3: continue
        d = np.array([abs(r["sdelta"]) for r in R]); G = np.array([r["G"] for r in R]); side = "+" if R[0]["sdelta"] > 0 else "−"
        (g0, s), *_ = np.linalg.lstsq(np.c_[np.ones_like(d), d], G, rcond=None)
        print(f"FIT {name.rstrip(chr(43))}{side} a={a:>2}: G_tip={g0:.5f} slope={s:.3f} n={len(R)} δ∈[{d.min():.1e},{d.max():.1e}]", flush=True)
