"""Level-2 test, route 1: κ_lim(p, 1) = (ι((e^{−2πi/p}𝒫₀)^p) − ½)/p for large p (general-germ normal form),
then p → ∞.  Horn coefficients a_1..a_30 at 90 digits (h = 1/4, N = 128).  Writes data/level2_bounded_p.tsv."""
import sys, time, pickle, os
import mpmath as mp
from bulbford.horn import horn_coeffs
from bulbford.germ import horn_germ, normal_form_index
mp.mp.dps = 60
cache = "data/horn_coeffs_z2_M30_dps90.pkl"
if os.path.exists(cache):
    a = pickle.load(open(cache, "rb"))
else:
    t0 = time.time(); a, al = horn_coeffs(30, h=0.25, N=128, R=1500, n=1800, dps=90, K=50)
    a = [mp.mpc(x) for x in a]; pickle.dump(a, open(cache, "wb"))
    print(f"horn coefficients: {time.time()-t0:.0f}s, aliasing {mp.nstr(al, 2)}, |a_30| = {mp.nstr(abs(a[30]), 3)}", flush=True)
for p in map(int, sys.argv[1:]):
    t0 = time.time()
    G = horn_germ(a, 2 * p + 1, mp.expjpi(-2 * mp.mpf(1) / p))
    k = (normal_form_index(G, p) - mp.mpf(1) / 2) / p
    line = f"{p}\t{mp.nstr(k.real, 40)}\t{mp.nstr(k.imag, 40)}\t{time.time()-t0:.0f}s"
    open("data/level2_bounded_p.tsv", "a").write(line + "\n"); print(line, flush=True)
