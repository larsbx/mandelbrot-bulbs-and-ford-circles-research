"""Sweep p over 1..q/2 for given q: x*, Taylor r_0..r_K of R_q(u), G_ant, G_cen, series-reconstructed
|u_a|/2 and |u_c|.  Writes experiments/data/taylor_q<q>.json."""
from __future__ import annotations
from paths import DATA
import json, sys, time
import numpy as np
from bulbford.cf import xstar, coprime_numerators, cf
from bulbford.dynamics import MAIN2, bulb
from bulbford.taylor import taylor

K = 16

def row(p, q):
    t = taylor(p, q, r=2.5, N=256)
    b = bulb(MAIN2, p, q)
    ua, uc = t.solve(-1, 2.0), t.solve(0, 1.0)
    return dict(p=p, q=q, xstar=float(xstar(p, q)), cf=cf(p, q),
                r=[[float(z.real), float(z.imag)] for z in t.coeffs[:K]],
                G_ant=b.G_ant, G_cen=b.G_cen, ua=[ua.real, ua.imag], uc=[uc.real, uc.imag],
                G_ant_series=abs(ua) / 2, G_cen_series=abs(uc),
                tail=float(np.abs(t.coeffs[K:K + 8]).max()))

if __name__ == "__main__":
    for q in map(int, sys.argv[1:]):
        t0 = time.time()
        rows = [row(p, q) for p in coprime_numerators(q) if p <= q // 2]
        json.dump(rows, open(DATA / f"taylor_q{q}.json", "w"))
        print(f"q={q}: {len(rows)} bulbs, {time.time() - t0:.1f}s", flush=True)
