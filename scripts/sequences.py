"""κ and Ĝ along sequences x̃ = p̄/q = [0; a1, …, a_k] with a_k = N → ∞ (one-sided limits at rationals),
plus a quadratic-irrational sequence. p = p̄⁻¹ mod q.  Output: data/sequences.json"""
from __future__ import annotations
import json, sys
from bulbford.cf import from_cf, modinv
from bulbford.taylor import kappa_fft

SEQS = {
    "1/3-  [0;3,N]":        lambda N: (0, 3, N),
    "1/3-  [0;3,N,5]":      lambda N: (0, 3, N, 5),
    "1/3+  [0;2,1,N]":      lambda N: (0, 2, 1, N),
    "1/3+  [0;2,1,N,7]":    lambda N: (0, 2, 1, N, 7),
    "1/2-  [0;2,N]":        lambda N: (0, 2, N),
    "1/2-  [0;2,N,3]":      lambda N: (0, 2, N, 3),
    "0+    [0;N]":          lambda N: (0, N),
    "0+    [0;N,4]":        lambda N: (0, N, 4),
    "2-φ   [0;2,1^k]":      lambda N: (0, 2) + (1,) * N,
}
NS = (3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384)

def entry(name, N):
    pb, q = from_cf(SEQS[name](N))
    if q > 2100 or q < 3: return None
    p = modinv(pb, q)
    t = kappa_fft(p, q)
    ua = t.solve(-1, 2.0)
    return dict(seq=name, N=N, p=p, q=q, xt=(pb if pb <= q / 2 else pb - q) / q,
                kappa=[t.coeffs[2].real, t.coeffs[2].imag], r3=[t.coeffs[3].real, t.coeffs[3].imag],
                G=abs(ua) / 2)

if __name__ == "__main__":
    names = sys.argv[1:] or list(SEQS)
    out = [e for name in names for N in NS if (e := entry(name, N))]
    json.dump(out, open("data/sequences.json", "w"))
    for e in out:
        print(f"{e['seq']:<18} N={e['N']:>3} q={e['q']:>5} x̃={e['xt']:+.5f}  κ={e['kappa'][0]:+.5f}{e['kappa'][1]:+.5f}i  r3={e['r3'][0]:+.5f}{e['r3'][1]:+.5f}i  Ĝ={e['G']:.5f}", flush=True)
