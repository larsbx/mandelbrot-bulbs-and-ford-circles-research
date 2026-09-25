"""κ(1/q) from the O(q²) normal form at fixed working precision; appends q, Re, Im (40 digits), dps, time
to data/p1_nf.tsv.  Usage: p1_nf.py DPS q1 q2 …"""
import sys, time, mpmath as mp
from bulbford.normal_form import normal_form
dps = int(sys.argv[1])
for q in map(int, sys.argv[2:]):
    t0 = time.time(); nf = normal_form(1, q, dps=dps)
    with mp.workdps(dps): k = (nf.iota - mp.mpf(1) / 2) / q
    line = f"{q}\t{mp.nstr(k.real, 40)}\t{mp.nstr(k.imag, 40)}\t{dps}\t{time.time()-t0:.0f}s"
    open("data/p1_nf.tsv", "a").write(line + "\n"); print(line, flush=True)
