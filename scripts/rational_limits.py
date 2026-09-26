"""V50: one-sided limits of κ at rationals.  For p/q = [0; a₁, …, a_n] and a tail t, the cardioid satellites
[0; a₁, …, a_n, N, t] accumulate at the root of W_{p/q} from above (n even) or below (n odd) as N → ∞.  Prediction: the
bounded-p limit (6) of the upper (above) or lower (below) parabolic renormalization of e^{2πip/q}z + z², with p_t, r_t
read off the tail ([0; t] = r_t/p_t; t = ∅: r_t = 0, the horn index itself) and phase μ_t above, μ̄_t below (V55).  Bulbs: normal form in midpoint arithmetic (128 bits), fits of degree 3 and 4 in 1/q.
Writes data/rational_limits.txt."""
from math import gcd
import mpmath as mp
from bulbford.cf import from_cf
from bulbford.normal_form_ball import normal_form_ball
from bulbford.extrapolate import power_fit
from bulbford.renorm import bounded_p_limit
from bulbford.horn_rational import root, horn_index_pq

mp.mp.dps = 30
NS = (128, 256, 512, 1024, 2048)
CASES = [((3,), ()), ((2, 1), ()), ((4,), ()), ((3, 1), ()), ((2, 1, 1), ()), ((2, 2), ()), ((2, 3), ()), ((2, 2, 1), ()),
         ((2, 1), (2,)), ((3,), (2,)), ((2, 1), (3,)),
         ((2,), (3,)), ((3,), (3,)), ((1,), (3,)), ((3,), (1, 2))]                # V55: lower sides, p_t ≥ 3


def kappa(p, q):
    k = normal_form_ball(p, q, prec=128, certified=False).kappa
    return mp.mpc(k.real.mid().str(25, radius=False), k.imag.mid().str(25, radius=False))


with open("data/rational_limits.txt", "w") as out:
    for pre, tail in CASES:
        P, Q = from_cf((0,) + pre)
        upper = len(pre) % 2 == 0
        rows = [(q, kappa(p, q)) for p, q in (from_cf((0,) + pre + (N,) + tail) for N in NS)]
        f3, f4 = power_fit(rows, 3)[0], power_fit(rows, 4)[0]
        if tail:
            tr, tp = from_cf((0,) + tail)                     # [0; t] = tr/tp, so [0; N, t] = tp/(tp·N + tr)
            pred = bounded_p_limit(tp, tr, dps=60, M=20, germ=root(P, Q, upper))
        else:
            pred = horn_index_pq(P, Q, upper=upper, dps=30)
        line = (f"{P}/{Q}\t{'above' if upper else 'below'}\t{','.join(map(str, tail)) or '-'}\t{mp.nstr(pred, 14)}\t"
                f"{mp.nstr(f4, 14)}\t{mp.nstr(abs(f4 - pred), 2)}\t{mp.nstr(abs(f4 - f3), 2)}")
        out.write(line + "\n"); out.flush(); print(line, flush=True)
