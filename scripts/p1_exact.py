"""Rigorous κ(1/q) = (ι_{1/q} − ½)/q in ball arithmetic for the given q; appends to data/p1_exact.tsv
(q, Re mid, Im mid, radius) with 40 significant digits."""
import sys, time
from bulbford.index import index
for q in map(int, sys.argv[1:]):
    t0 = time.time(); r = index(1, q); k = (r - 0.5) / q
    line = f"{q}\t{k.real.str(40, radius=False)}\t{k.imag.str(40, radius=False)}\t{float(k.rad()):.1e}\t{time.time()-t0:.0f}s"
    open("data/p1_exact.tsv", "a").write(line + "\n"); print(line, flush=True)
