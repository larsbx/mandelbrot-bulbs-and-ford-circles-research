"""G and κ for p/q = [0; a1, N, 3] as a1 grows (t = p/q → 0 at fixed tail), main2 and disk2."""
import json, sys
from bulbford.cf import from_cf, modinv
from bulbford.dynamics import MAIN2, DISK2
from bulbford.taylor import kappa_fft
N = int(sys.argv[1]) if len(sys.argv) > 1 else 16
rows = []
for a1 in (1, 2, 3, 4, 5, 7, 9, 12, 16, 24, 32, 48, 64, 100, 150, 200):
    p, q = from_cf((0, a1, N, 3)) if a1 > 1 else from_cf((0, N, 3))
    pb = modinv(p, q); xt = (pb if pb <= q / 2 else pb - q) / q
    out = {}
    for fam in (MAIN2, DISK2):
        t = kappa_fft(p, q, fam); out[fam.name] = (t.coeffs[2], abs(t.solve(-1, 2.0)) / 2)
    rows.append(dict(a1=a1, p=p, q=q, t=p / q, xt=xt, **{k: [v[0].real, v[0].imag, v[1]] for k, v in out.items()}))
    print(f"a1={a1:>3} p/q={p:>4}/{q:<5} t={p/q:.4f} x̃={xt:+.4f} | main2 κ={out['main2'][0]:+.5f} G={out['main2'][1]:.5f} | disk2 κ={out['disk2'][0]:+.5f} G={out['disk2'][1]:.5f} | ratio={out['disk2'][1]/out['main2'][1]:.5f}", flush=True)
json.dump(rows, open(f"data/a1_scan_N{N}.json", "w"))
