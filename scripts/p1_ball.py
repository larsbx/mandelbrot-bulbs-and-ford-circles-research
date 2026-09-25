"""Certified κ(1/q) (Arb balls, blocked convolution; bulbford/normal_form_ball.py).
Appends q, Re mid, Im mid (45 digits), radius, bits, time to data/p1_ball.tsv."""
import sys, time
from flint import ctx
from bulbford.normal_form_ball import normal_form_ball
for q in map(int, sys.argv[1:]):
    t0 = time.time(); nb = normal_form_ball(1, q); k = nb.kappa
    line = f"{q}\t{k.real.mid().str(45, radius=False)}\t{k.imag.mid().str(45, radius=False)}\t{float(k.rad()):.2e}\t{time.time()-t0:.0f}s"
    open("data/p1_ball.tsv", "a").write(line + "\n"); print(line, flush=True)
