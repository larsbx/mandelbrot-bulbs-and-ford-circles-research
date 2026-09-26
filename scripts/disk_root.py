"""C23 (iv): the period-2 disc near its root c = −3/4 against the cardioid near its ½-root (same parabolic germ).
Disk [0; N, t] (κ from Cauchy coefficients in the 2-cycle multiplier) vs cardioid [0; 1, 1, N − 1, t] (normal form),
N → ∞ with the tail t fixed; also G for t = ().  Fits in 1/q.  Writes data/disk_root.txt."""
import mpmath as mp
from bulbford.taylor import taylor
from bulbford.dynamics import DISK2, MAIN2, bulb
from bulbford.cf import from_cf
from bulbford.normal_form_ball import normal_form_ball
from bulbford.extrapolate import power_fit

mp.mp.dps = 30
NS = (128, 256, 512, 1024, 2048)


def kappa_main(p, q):
    k = normal_form_ball(p, q, prec=128, certified=False).kappa
    return mp.mpc(k.real.mid().str(25, radius=False), k.imag.mid().str(25, radius=False))


def limits(disk, card):
    return [(J, power_fit(disk, J)[0], power_fit(card, J)[0]) for J in (2, 3)]


with open("data/disk_root.txt", "w") as out:
    def emit(s):
        out.write(s + "\n"); out.flush(); print(s, flush=True)
    for tail in ((), (2,), (3,)):
        disk, card = [], []
        for N in NS:
            p, q = from_cf((0, N) + tail)
            P, Q = from_cf((0, 1, 1, N - 1) + tail)
            disk.append((q, complex(taylor(p, q, DISK2, r=1.5, N=64).coeffs[2])))
            card.append((Q, kappa_main(P, Q)))
        for J, d, c in limits(disk, card):
            emit(f"kappa tail {tail}: J = {J}: disk {mp.nstr(d, 11)}  cardioid {mp.nstr(c, 11)}  |Δ| {mp.nstr(abs(d - c), 2)}")
    disk = [(N, bulb(DISK2, 1, N).G_ant) for N in NS[1:] + (4096,)]
    card = [(2 * N - 1, bulb(MAIN2, N, 2 * N - 1).G_ant) for N in NS[1:] + (4096,)]
    for J, d, c in limits(disk, card):
        emit(f"G tail (): J = {J}: disk {mp.nstr(mp.re(d), 11)}  cardioid {mp.nstr(mp.re(c), 11)}  |Δ| {mp.nstr(abs(d - c), 2)}")
