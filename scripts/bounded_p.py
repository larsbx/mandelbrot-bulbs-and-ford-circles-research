"""Bounded-p excess Ĥ(p) = lim G(bounded-p sequence x̃ = [0; tail, N]) − Ĝ_tip(p'/q', same side)."""
import numpy as np
from bulbford.cf import from_cf, modinv
from bulbford.taylor import kappa_fft
TAILS = {"1/3": (3,), "1/2": (2,), "1/4": (4,), "1/4+": (3, 1), "1/5": (5,), "2/5": (2, 2), "2/5+": (2, 1, 1),
         "2/7": (3, 2), "1/6": (6,), "1/7": (7,), "3/7": (2, 3), "3/8": (2, 1, 2), "2/9": (4, 2), "4/9": (2, 4)}
TIPS = {"1/3": 1.12056, "1/2": 1.10008, "1/4": 1.11829, "1/4+": 1.11696, "1/5": 1.10953, "2/5": 1.14032, "2/5+": 1.14196,
        "2/7": 1.13485, "1/6": 1.09968, "1/7": 1.09048, "3/7": 1.14517, "3/8": 1.14769, "2/9": 1.12328, "4/9": 1.14466}
print("rational | p | G(N=128), G(N=256), G(N=512) | extrap (c/N) | Ĝ_tip | Ĥ(p)")
for name, tail in TAILS.items():
    vals = []
    for N in (128, 256, 512):
        pb, q = from_cf((0,) + tail + (N,)); p = modinv(pb, q); t = kappa_fft(p, q); vals.append((N, p, abs(t.solve(-1, 2.0)) / 2))
    (N1, p, G1), (N2, _, G2), (N3, _, G3) = vals
    ext = 2 * G3 - G2      # Richardson for a c/N approach
    print(f"  {name:>4} | {p:>2} | {G1:.5f}, {G2:.5f}, {G3:.5f} | {ext:.5f} | {TIPS[name]:.5f} | {ext - TIPS[name]:+.4f}", flush=True)
