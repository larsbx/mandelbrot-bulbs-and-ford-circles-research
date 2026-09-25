"""Bulb-side κ(p/(pN+1)) for larger p by the normal form, to test C23 (ii) beyond p ≤ 5.
Usage: bounded_p_large.py p N1 N2 …  → data/bounded_p_large.tsv (p, N, q, Re, Im)."""
import sys, time, mpmath as mp
from bulbford.normal_form import normal_form
p = int(sys.argv[1])
for N in map(int, sys.argv[2:]):
    q = p * N + 1; t0 = time.time()
    nf = normal_form(p, q, dps=50)
    with mp.workdps(50): k = (nf.iota - mp.mpf(1) / 2) / q
    line = f"{p}\t{N}\t{q}\t{mp.nstr(k.real, 30)}\t{mp.nstr(k.imag, 30)}\t{time.time()-t0:.0f}s"
    open("data/bounded_p_large.tsv", "a").write(line + "\n"); print(line, flush=True)
