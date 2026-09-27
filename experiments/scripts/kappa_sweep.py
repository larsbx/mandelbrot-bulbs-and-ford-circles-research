"""Dense κ(p/q) via FFT for all p ≤ q/2 (double precision, N=64). Writes experiments/data/kappa_q<q>.json."""
from paths import DATA
import json, sys, time
from bulbford.cf import coprime_numerators, xstar, modinv
from bulbford.taylor import kappa_fft

if __name__ == "__main__":
    for q in map(int, sys.argv[1:]):
        t0 = time.time(); rows = []
        for p in coprime_numerators(q):
            if p > q // 2: continue
            t = kappa_fft(p, q); pb = modinv(p, q); ua = t.solve(-1, 2.0)
            rows.append(dict(p=p, q=q, xt=(pb if pb <= q / 2 else pb - q) / q, kappa=[t.coeffs[2].real, t.coeffs[2].imag],
                             r3=[t.coeffs[3].real, t.coeffs[3].imag], G=abs(ua) / 2))
        json.dump(rows, open(DATA / f"kappa_q{q}.json", "w"))
        print(f"q={q}: {len(rows)} rows, {time.time()-t0:.0f}s", flush=True)
