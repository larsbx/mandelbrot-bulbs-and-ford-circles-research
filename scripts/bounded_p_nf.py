"""κ(p/q) for fixed p, q = pN + r (all r coprime to p), N ∈ NS, by the normal form (50 digits).
Appends p, r, N, q, Re, Im to data/bounded_p_nf.tsv."""
import sys, time, mpmath as mp
from math import gcd
from bulbford.normal_form import normal_form
ps = [int(x) for x in sys.argv[1].split(",")]; NS = [int(x) for x in sys.argv[2].split(",")]
for N in NS:
    for p in ps:
        for r in (r for r in range(1, p) if gcd(p, r) == 1):
            q = p * N + r; t0 = time.time()
            nf = normal_form(p, q, dps=50)
            with mp.workdps(50): k = (nf.iota - mp.mpf(1) / 2) / q
            line = f"{p}\t{r}\t{N}\t{q}\t{mp.nstr(k.real, 30)}\t{mp.nstr(k.imag, 30)}\t{time.time()-t0:.0f}s"
            open("data/bounded_p_nf.tsv", "a").write(line + "\n"); print(line, flush=True)
