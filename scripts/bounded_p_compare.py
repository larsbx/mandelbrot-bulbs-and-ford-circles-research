"""Extrapolate κ(p/(pN+r)) in 1/q (fits with J = 2, 3 over the largest N) and compare with bounded_p_limit."""
import mpmath as mp
from collections import defaultdict
from bulbford.renorm import bounded_p_limit
from bulbford.extrapolate import power_fit
mp.mp.dps = 40
seq = defaultdict(dict)                                   # keyed by q: reruns are de-duplicated
for l in open("data/bounded_p_nf.tsv"):
    p, r, N, q, re, im, *_ = l.split("\t"); seq[(int(p), int(r))][int(q)] = mp.mpc(re, im)
fit = lambda R, J: power_fit(R, J, dps=40)
print(" p r | extrapolated lim κ (J=3, all N)   | J=2 (largest 4) | predicted (horn map)            | |Δ|      | c₁ (coefficient of 1/q)")
for (p, r), D in sorted(seq.items()):
    R = sorted(D.items())
    if len(R) < 4: continue
    c3, c2 = fit(R, min(3, len(R) - 1)), fit(R[-4:], 2)
    pred = bounded_p_limit(p, r)
    print(f" {p} {r} | {mp.nstr(c3[0], 14):>32} | {mp.nstr(abs(c3[0] - c2[0]), 2):>8} | {mp.nstr(pred, 14):>32} | {mp.nstr(abs(c3[0] - pred), 2):>8} | {mp.nstr(c3[1], 6)}")
