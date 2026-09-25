"""z³+c main component: bounded-p and generic limits of κ, G at x̃ → ⅓⁻, ½⁻, 0⁺ (cf. V16 for main2)."""
import json
from bulbford.cf import from_cf, modinv
from bulbford.dynamics import MAIN3, MAIN2
from bulbford.taylor import kappa_fft
SEQ = {"1/3- p=3": lambda N: (0, 3, N), "1/3- generic": lambda N: (0, 3, N, 16),
       "1/2- p=2": lambda N: (0, 2, N), "1/2- generic": lambda N: (0, 2, N, 16),
       "0+ p=1": lambda N: (0, N), "0+ generic": lambda N: (0, N, 16)}
rows = []
for name, f in SEQ.items():
    for N in (32, 64, 128):
        pb, q = from_cf(f(N)); p = modinv(pb, q)
        out = {}
        for fam in (MAIN2, MAIN3):
            t = kappa_fft(p, q, fam); out[fam.name] = (t.coeffs[2], abs(t.solve(-1, 2.0)) / 2)
        rows.append(dict(seq=name, N=N, p=p, q=q, **{k: [v[0].real, v[0].imag, v[1]] for k, v in out.items()}))
        print(f"{name:<13} N={N:>3} p/q={p}/{q:<5} | main2 κ={out['main2'][0]:+.5f} G={out['main2'][1]:.5f} | main3 κ={out['main3'][0]:+.5f} G={out['main3'][1]:.5f} | G3/G2={out['main3'][1]/out['main2'][1]:.4f}", flush=True)
json.dump(rows, open("data/main3_limits.json", "w"))
