"""Second-variable test: p/q = [0; prefix…, N, 3] (so x̃ → 1/3⁻ as N → ∞, prefix fixes t = p/q ≈ [0; prefix]).
If lim G depends smoothly on t only, the limits line up on one curve in t whatever the prefix length."""
from paths import DATA
import json, sys
from bulbford.cf import from_cf, xstar, modinv
from bulbford.taylor import kappa_fft

PREFIXES = [(), (9,), (5,), (4,), (3,), (3, 2), (2,), (2, 2), (2, 3), (2, 5), (1, 1, 2), (1, 2), (1, 3), (1, 5)]
NS = (64, 128, 256)
TAIL = tuple(int(a) for a in sys.argv[1:]) or (3,)

rows = []
for pre in PREFIXES:
    for N in NS:
        p, q = from_cf((0,) + pre + (N,) + TAIL)
        if q > 3000: continue
        pb = modinv(p, q); t = kappa_fft(p, q); ua = t.solve(-1, 2.0)
        r = dict(prefix=pre, N=N, p=p, q=q, t=p / q, xt=(pb if pb <= q / 2 else pb - q) / q,
                 kappa=[t.coeffs[2].real, t.coeffs[2].imag], G=abs(ua) / 2)
        rows.append(r)
        print(f"prefix={str(pre):<10} N={N:>3} p/q={p:>5}/{q:<5} t={r['t']:.4f} x̃={r['xt']:+.5f}  κ={r['kappa'][0]:+.5f}{r['kappa'][1]:+.5f}i  G={r['G']:.5f}", flush=True)
json.dump(rows, open(DATA / f"prefix_test_tail{'_'.join(map(str, TAIL))}.json", "w"))
