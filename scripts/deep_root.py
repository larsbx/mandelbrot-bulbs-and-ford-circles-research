"""V53: the root principle at a deeper root.  The disc's own 1/3-satellite is attached at c = −1 + e^{2πi/3}/4, where
the germ is the first return f² at the parabolic 2-cycle (three petals, horn_rational.two_cycle_germ).  Its upper and
lower parabolic renormalizations against (a) the disc's satellites approaching that root from above/below and (b) the
satellites at 1/q, (q−1)/q of the period-6 component attached there.  Writes data/deep_root.txt."""
import mpmath as mp
from bulbford.horn_rational import horn_index_pq, two_cycle_germ
from bulbford.cf import from_cf
from bulbford.taylor import taylor
from bulbford.dynamics import DISK2, component, bulb
from bulbford.extrapolate import power_fit

mp.mp.dps = 30
NS = (128, 256, 512, 1024, 2048)
germ = two_cycle_germ(-1 + mp.expjpi(mp.mpf(2) / 3) / 4)
pred = {up: horn_index_pq(1, 3, upper=up, dps=30, poly=germ) for up in (True, False)}
W6 = component("W1/2,1/3", 6, bulb(DISK2, 1, 3).cen.c)
cases = [("disc [0;2,1,N]", True, DISK2, lambda N: from_cf((0, 2, 1, N))),
         ("disc [0;3,N]", False, DISK2, lambda N: from_cf((0, 3, N))),
         ("period-6 1/q", True, W6, lambda N: (1, N)),
         ("period-6 (q-1)/q", False, W6, lambda N: (N - 1, N))]
with open("data/deep_root.txt", "w") as out:
    for label, up, fam, pq in cases:
        rows = [(q, complex(taylor(p, q, fam, r=1.5, N=64).coeffs[2])) for p, q in map(pq, NS)]
        f2, f3 = power_fit(rows, 2)[0], power_fit(rows, 3)[0]
        line = (f"{label}\t{'upper' if up else 'lower'}\t{mp.nstr(pred[up], 12)}\t{mp.nstr(f3, 12)}\t"
                f"{mp.nstr(abs(f3 - pred[up]), 2)}\t{mp.nstr(abs(f3 - f2), 2)}")
        out.write(line + "\n"); out.flush(); print(line, flush=True)
