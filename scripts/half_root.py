"""V49: the ½-root constant from the horn map of f(z) = −z + z² (bulbford/horn_half.py) against the bulbs (V48).
Writes data/half_root.txt."""
import mpmath as mp
from bulbford.horn_half import HALF, kappa_half
from bulbford.renorm import bounded_p_limit, G_limits
from bulbford.dynamics import DISK2, bulb
from bulbford.extrapolate import power_fit

mp.mp.dps = 40
with open("data/half_root.txt", "w") as out:
    def emit(s):
        out.write(s + "\n"); print(s, flush=True)
    for h in (4.5, 5.5):
        a1, a2, k, al = kappa_half(h=h, dps=100, K=40, n=900)
        emit(f"h = {h}: K_1/2 = a2/(2 pi i a1^2) = {mp.nstr(k, 32)}   (aliasing {mp.nstr(al, 2)})")
    for p, r in ((1, 1), (2, 1), (3, 1), (3, 2)):
        emit(f"bounded p = {p}, r = {r}: kappa limit {mp.nstr(bounded_p_limit(p, r, dps=60, M=20, germ=HALF), 13)}")
    for p, r in ((1, 1), (2, 1)):
        emit(f"G limit p = {p}, r = {r}: {G_limits(p, r, germ=HALF)[0]:.10f}")
    rows = [(2 * N + 1, bulb(DISK2, 2, 2 * N + 1).G_ant) for N in (256, 512, 1024, 2048, 4096)]
    emit(f"bulbs: disc G along 2/(2N+1), N = 256..4096, J = 3: {mp.nstr(mp.re(power_fit(rows, 3)[0]), 10)}")
    emit("bulbs (data/disk_root.txt, J = 3): " + "; ".join(l.strip() for l in open("data/disk_root.txt") if "J = 3" in l))
