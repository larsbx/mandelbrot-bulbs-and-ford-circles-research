"""Bulb-side G(p/(pN+r)) for the listed p, N (FFT route); appends p, r, N, q, G to data/G_bounded_p.tsv."""
import sys
from math import gcd
from bulbford.taylor import kappa_fft
ps = [int(x) for x in sys.argv[1].split(",")]; NS = [int(x) for x in sys.argv[2].split(",")]
for N in NS:
    for p in ps:
        for r in ((r for r in range(1, p) if gcd(p, r) == 1) if p > 1 else [1]):
            q = p * N + r; G = abs(kappa_fft(p, q).solve(-1, 2.0)) / 2
            line = f"{p}\t{r}\t{N}\t{q}\t{G:.12f}"; open("data/G_bounded_p.tsv", "a").write(line + "\n"); print(line, flush=True)
