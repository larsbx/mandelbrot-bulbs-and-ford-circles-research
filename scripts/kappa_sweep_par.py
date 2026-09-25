"""Parallel version of kappa_sweep: chunk k of n (p ≡ k mod n) → data/kappa_q<q>_chunk<k>.json;
merge with: python3 scripts/kappa_sweep_par.py <q> merge <n>."""
import json, sys, time, glob
from bulbford.cf import coprime_numerators, modinv
from bulbford.taylor import kappa_fft

q = int(sys.argv[1])
if sys.argv[2] == "merge":
    rows = sorted((r for f in sorted(glob.glob(f"data/kappa_q{q}_chunk*.json")) for r in json.load(open(f))), key=lambda r: r["p"])
    json.dump(rows, open(f"data/kappa_q{q}.json", "w")); print(f"q={q}: merged {len(rows)} rows")
else:
    k, n = int(sys.argv[2]), int(sys.argv[3]); t0 = time.time(); rows = []
    for p in coprime_numerators(q):
        if p > q // 2 or p % n != k: continue
        t = kappa_fft(p, q); pb = modinv(p, q); ua = t.solve(-1, 2.0)
        rows.append(dict(p=p, q=q, xt=(pb if pb <= q / 2 else pb - q) / q, kappa=[t.coeffs[2].real, t.coeffs[2].imag],
                         r3=[t.coeffs[3].real, t.coeffs[3].imag], G=abs(ua) / 2))
    json.dump(rows, open(f"data/kappa_q{q}_chunk{k}.json", "w")); print(f"q={q} chunk {k}/{n}: {len(rows)} rows, {time.time()-t0:.0f}s")
