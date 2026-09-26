"""V51: the root principle at a three-petal root.  Satellites of the period-3 component W_{1/3} (centre ≈ −0.1226 +
0.7449i, parametrised by the multiplier of its 3-cycle, dynamics.component) against the upper/lower parabolic
renormalizations of e^{2πi/3}z + z² at its root.  Writes data/period3_root.txt."""
import mpmath as mp
from bulbford.dynamics import component, bulb
from bulbford.taylor import taylor
from bulbford.extrapolate import power_fit
from bulbford.renorm import G_limits, bounded_p_limit
from bulbford.horn_rational import root, horn_index_pq

mp.mp.dps = 30
W3 = component("W1/3", 3, -0.1225611668766536 + 0.7448617666197442j)
QS = (128, 256, 512, 1024, 2048)


def fits(rows):
    return power_fit(rows, 2)[0], power_fit(rows, 3)[0]


with open("data/period3_root.txt", "w") as out:
    def emit(label, pred, rows):
        f2, f3 = fits(rows)
        line = f"{label}\t{mp.nstr(pred, 12)}\t{mp.nstr(f3, 12)}\t{mp.nstr(abs(f3 - pred), 2)}\t{mp.nstr(abs(f3 - f2), 2)}"
        out.write(line + "\n"); out.flush(); print(line, flush=True)
    for up, side in ((True, "1/q"), (False, "(q-1)/q")):
        num = (lambda q: 1) if up else (lambda q: q - 1)
        emit(f"kappa {side}", horn_index_pq(1, 3, upper=up, dps=30),
             [(q, complex(taylor(num(q), q, W3, r=1.5, N=64).coeffs[2])) for q in QS])
    emit("kappa 2/(2N+1)", bounded_p_limit(2, 1, dps=60, M=20, germ=root(1, 3, True)),
         [(2 * N + 1, complex(taylor(2, 2 * N + 1, W3, r=1.5, N=64).coeffs[2])) for N in QS])
    for up, side in ((True, "1/q"), (False, "(q-1)/q")):
        num = (lambda q: 1) if up else (lambda q: q - 1)
        emit(f"G {side}", G_limits(1, 1, germ=root(1, 3, up))[0],
             [(q, bulb(W3, num(q), q).G_ant) for q in (256, 512, 1024, 2048, 4096)])
